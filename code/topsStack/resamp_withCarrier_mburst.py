#!/usr/bin/env python3


import os
import argparse
import copy
import json
import datetime
import numpy as np

import isce
import isceobj
import stdproc
from isceobj.Util.Poly2D import Poly2D
from isceobj.Sensor.TOPS import createTOPSSwathSLCProduct
from stdproc.stdproc import crossmul
import s1a_isce_utils as ut


def createParser():
    parser = argparse.ArgumentParser( description='Resampling burst by burst SLCs ')

    # inputs
    parser.add_argument('-m', '--reference', dest='reference', type=str, required=True,
            help='Directory with reference acquisition')

    parser.add_argument('-s', '--secondary', dest='secondary', type=str, required=True,
            help='Directory with secondary acquisition')

    # output
    parser.add_argument('-o', '--coregdir', dest='coreg', type=str, default='coreg_secondary',
            help='Directory with coregistered SLCs and IFGs')

    # additional setups
    parser.add_argument('-a', '--azimuth_misreg', dest='misreg_az', type=str, default=0.0,
            help='File name with the azimuth misregistration')

    parser.add_argument('-r', '--range_misreg', dest='misreg_rng', type=str, default=0.0,
            help='File name with the range misregistration')

    parser.add_argument('--noflat', dest='noflat', action='store_true', default=False,
            help='To turn off flattening. False: flattens the SLC. True: turns off flattening.')

    parser.add_argument('-v', '--overlap', dest='overlap', action='store_true', default=False,
            help='Is this an overlap burst slc. default: False')

    parser.add_argument('-d', '--overlapDir', dest='overlapDir', type=str, default='overlap',
            help='reference overlap directory')

    return parser


def cmdLineParse(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)


def loadCropMetadata(path=os.path.join('geom_reference', 'crop_metadata.json')):
    """Load per-acquisition native burst grids written by pre_data_mburst1."""
    if not os.path.exists(path):
        raise RuntimeError(
            'Missing crop_metadata.json. Run pre_data_mburst1.py from run_00.')
    with open(path, 'r', encoding='utf-8') as f:
        payload = json.load(f)
    if payload.get('version') != 2 or not isinstance(
            payload.get('acquisitions'), dict):
        raise RuntimeError('Unsupported mburst crop metadata: {0}'.format(path))
    print('Using per-acquisition full-burst carrier metadata: {0}'.format(path))
    return payload


def acquisitionDate(path):
    value = os.path.basename(os.path.normpath(path))
    digits = ''.join(ch for ch in value if ch.isdigit())
    if len(digits) >= 8:
        return digits[:8]
    raise RuntimeError('Cannot determine secondary acquisition date: {0}'.format(path))


def getCropMetadata(payload, date, swath, burstNumber):
    iw = 'IW{0}'.format(swath)
    nn = '{0:02d}'.format(int(burstNumber))
    acquisitions = payload.get('acquisitions', {})
    date_info = acquisitions.get(str(date))
    if date_info is None:
        raise RuntimeError('Missing native grid metadata for acquisition {0}'.format(date))
    info = date_info.get(iw, {}).get(nn)
    if info is None:
        raise RuntimeError(
            'Missing native grid metadata for {0} {1} burst {2}'.format(
                date, iw, nn))
    return info


def _parseTime(value):
    parsed = datetime.datetime.fromisoformat(str(value).strip().replace('Z', '+00:00'))
    if parsed.tzinfo is not None:
        parsed = parsed.replace(tzinfo=None)
    return parsed


def estimateNativeCarrierPolynomials(product, burst, cropInfo, offset=0.0,
                                     xstep=500, ystep=50,
                                     azorder=5, rgorder=3):
    """Fit the secondary carrier on its original full-burst coordinates."""
    required = ('xoff', 'yoff', 'width', 'height', 'original_width',
                'original_height', 'original_starting_range',
                'original_sensing_start', 'azimuth_time_interval')
    missing = [key for key in required if cropInfo.get(key) is None]
    if missing:
        raise RuntimeError('Incomplete native carrier metadata: {0}'.format(', '.join(missing)))
    xoff = int(cropInfo['xoff'])
    yoff = int(cropInfo['yoff'])
    width = int(cropInfo['width'])
    height = int(cropInfo['height'])
    originalWidth = int(cropInfo['original_width'])
    originalHeight = int(cropInfo['original_height'])
    originalStartingRange = float(cropInfo['original_starting_range'])
    originalSensingStart = _parseTime(cropInfo['original_sensing_start'])
    dt = float(cropInfo['azimuth_time_interval'])
    if (width, height) != (burst.numberOfSamples, burst.numberOfLines):
        raise RuntimeError(
            'Stale native carrier window for burst {0}: metadata={1}x{2}, product={3}x{4}'.format(
                burst.burstNumber, width, height,
                burst.numberOfSamples, burst.numberOfLines))
    if (xoff < 0 or yoff < 0 or xoff + width > originalWidth or
            yoff + height > originalHeight):
        raise RuntimeError('Native carrier window exceeds original burst {0}'.format(burst.burstNumber))
    if abs(dt - burst.azimuthTimeInterval) > 1.0e-9:
        raise RuntimeError('Azimuth time interval changed for burst {0}'.format(burst.burstNumber))
    expectedRange = originalStartingRange + xoff * burst.rangePixelSize
    if abs(expectedRange - burst.startingRange) > 0.51 * burst.rangePixelSize:
        raise RuntimeError(
            'Cropped starting range is inconsistent with original xoff for burst {0}: {1} vs {2}'.format(
                burst.burstNumber, burst.startingRange, expectedRange))
    originalSensingStop = originalSensingStart + datetime.timedelta(
        seconds=(originalHeight - 1) * dt)
    originalSensingMid = originalSensingStart + 0.5 * (
        originalSensingStop - originalSensingStart)
    velocity = burst.orbit.interpolateOrbit(
        originalSensingMid, method='hermite').getVelocity()
    Vs = np.linalg.norm(velocity)
    Ks = 2 * Vs * burst.azimuthSteeringRate / burst.radarWavelength
    x = np.unique(np.append(
        np.arange(0, burst.numberOfSamples, xstep, dtype=int),
        burst.numberOfSamples - 1))
    y = np.unique(np.append(
        np.arange(0, burst.numberOfLines, ystep, dtype=int),
        burst.numberOfLines - 1))
    xx, yy = np.meshgrid(x, y)
    rng = burst.startingRange + xx * burst.rangePixelSize
    eta = ((yy + yoff) - (originalHeight // 2)) * dt + offset * dt
    f_etac = burst.doppler(rng)
    Ka = burst.azimuthFMRate(rng)
    anchor = (burst.doppler(originalStartingRange) /
              burst.azimuthFMRate(originalStartingRange))
    eta_ref = anchor - (f_etac / Ka)
    Kt = Ks / (1.0 - Ks / Ka)
    carrier = np.pi * Kt * ((eta - eta_ref) ** 2)
    poly = Poly2D()
    poly.initPoly(rangeOrder=rgorder, azimuthOrder=azorder)
    poly.polyfit(xx.flatten() + 1, yy.flatten() + 1, carrier.flatten())
    poly.createPoly2D()
    fit = poly(yy + 1, xx + 1)
    diff = carrier - fit
    maxdiff = float(np.max(np.abs(diff)))
    print('Native-grid carrier burst {0}: xoff={1}, yoff={2}, full={3}x{4}, misfit={5:.6g} rad'.format(
        burst.burstNumber, xoff, yoff, originalWidth, originalHeight, maxdiff))
    if maxdiff > 0.01:
        print('Warning: native-grid azimuth carrier polynomial may not be accurate enough')
    dop = burst.doppler
    dpoly = Poly2D()
    dpoly._meanRange = (dop._mean - burst.startingRange) / burst.rangePixelSize
    dpoly._normRange = dop._norm / burst.rangePixelSize
    coeffs = [2 * np.pi * val * burst.azimuthTimeInterval for val in dop._coeffs]
    dpoly.initPoly(rangeOrder=dop._order, azimuthOrder=0)
    dpoly.setCoeffs([coeffs])
    return poly, dpoly


def resampSecondary(ref, sec, rdict, outname, flatten):
    '''
    Resample burst by burst.
    '''

    azpoly = rdict['azpoly']
    rgpoly = rdict['rgpoly']
    azcarrpoly = rdict['carrPoly']
    dpoly = rdict['doppPoly']

    rngImg = isceobj.createImage()
    rngImg.load(rdict['rangeOff'] + '.xml')
    rngImg.setAccessMode('READ')

    aziImg = isceobj.createImage()
    aziImg.load(rdict['azimuthOff'] + '.xml')
    aziImg.setAccessMode('READ')

    inimg = isceobj.createSlcImage()
    inimg.load(sec.image.filename + '.xml')
    inimg.setAccessMode('READ')
    
    ######################################################################################################################################
    #import sys
    #sys.path.append('/home/jovyan/iscepredata/p')
    #import stdproc0 as stdproc
    ######################################################################################################################################
    rObj = stdproc.createResamp_slc()

    

    rObj.slantRangePixelSpacing = sec.rangePixelSize
    rObj.radarWavelength = sec.radarWavelength
    rObj.azimuthCarrierPoly = azcarrpoly
    rObj.dopplerPoly = dpoly

    rObj.azimuthOffsetsPoly = azpoly
    rObj.rangeOffsetsPoly = rgpoly
    rObj.imageIn = inimg
    

    width = ref.numberOfSamples
    length = ref.numberOfLines
  

    
    
    imgOut = isceobj.createSlcImage()
    imgOut.setWidth(width)
    imgOut.filename = outname
    imgOut.setAccessMode('write')
    
    rObj.outputWidth = width
    rObj.outputLines = length

    rObj.residualRangeImage = rngImg
    rObj.residualAzimuthImage = aziImg
    rObj.flatten = flatten
    rObj.resamp_slc(imageOut=imgOut)

    imgOut.renderHdr()
    imgOut.renderVRT()
    return imgOut


def main(iargs=None):
    '''
    Create coregistered overlap secondarys.
    '''
    inps = cmdLineParse(iargs)
    cropMetadata = loadCropMetadata()
    secondaryDate = acquisitionDate(inps.secondary)
    referenceSwathList = ut.getSwathList(inps.reference)
    secondarySwathList = ut.getSwathList(inps.secondary)
    swathList = list(sorted(set(referenceSwathList + secondarySwathList)))

    for swath in swathList:

        ####Load secondary metadata
        reference = ut.loadProduct( os.path.join(inps.reference , 'IW{0}.xml'.format(swath)))
        secondary = ut.loadProduct( os.path.join(inps.secondary , 'IW{0}.xml'.format(swath)))
        if inps.overlap:
            referenceTop = ut.loadProduct(os.path.join(inps.reference, inps.overlapDir , 'IW{0}_top.xml'.format(swath)))
            referenceBottom = ut.loadProduct(os.path.join(inps.reference, inps.overlapDir , 'IW{0}_bottom.xml'.format(swath)))

        dt = secondary.bursts[0].azimuthTimeInterval
        dr = secondary.bursts[0].rangePixelSize

        if os.path.exists(str(inps.misreg_az)):
             with open(inps.misreg_az, 'r') as f:
                misreg_az = float(f.readline())
        else:
             misreg_az = 0.0

        if os.path.exists(str(inps.misreg_rng)):
             with open(inps.misreg_rng, 'r') as f:
                misreg_rg = float(f.readline())
        else:
             misreg_rg = 0.0


        ###Output directory for coregistered SLCs
        if not inps.overlap:
            outdir = os.path.join(inps.coreg,'IW{0}'.format(swath))
            offdir = os.path.join(inps.coreg,'IW{0}'.format(swath)) 
        else:
            outdir = os.path.join(inps.coreg, inps.overlapDir, 'IW{0}'.format(swath))
            offdir = os.path.join(inps.coreg, inps.overlapDir, 'IW{0}'.format(swath))
        os.makedirs(outdir, exist_ok=True)


        ####Indices w.r.t reference
        burstoffset, minBurst, maxBurst = ut.getCommonBurstLimits(reference, secondary)
        secondaryBurstStart = minBurst +  burstoffset
        secondaryBurstEnd = maxBurst

        relShifts = ut.getRelativeShifts(reference, secondary, minBurst, maxBurst, secondaryBurstStart)
        print('Shifts: ', relShifts)
        if inps.overlap:
            maxBurst = maxBurst - 1 ###For overlaps


        ####Can corporate known misregistration here

        apoly = Poly2D()
        apoly.initPoly(rangeOrder=0,azimuthOrder=0,coeffs=[[0.]])

        rpoly = Poly2D()
        rpoly.initPoly(rangeOrder=0,azimuthOrder=0,coeffs=[[0.]])


        #topCoreg = createTOPSSwathSLCProduct()
        topCoreg = ut.coregSwathSLCProduct()
        topCoreg.configure()

        if inps.overlap:
            botCoreg = ut.coregSwathSLCProduct()
            botCoreg.configure()

        for ii in range(minBurst, maxBurst):
            jj = secondaryBurstStart + ii - minBurst

            if inps.overlap:
                botBurst = referenceBottom.bursts[ii]
                topBurst = referenceTop.bursts[ii]
            else:
                topBurst = reference.bursts[ii]


            secBurst = secondary.bursts[jj]

            #####Top burst processing
            try:
                offset = relShifts[jj]
            except:
                raise Exception('Trying to access shift for secondary burst index {0}, which may not overlap with reference'.format(jj))

            if inps.overlap:
                outname = os.path.join(outdir, 'burst_top_%02d_%02d.slc'%(ii+1,ii+2))

                ####Setup initial polynomials
                ### If no misregs are given, these are zero
                ### If provided, can be used for resampling without running to geo2rdr again for fast results
                rdict = {'azpoly' : apoly,
                         'rgpoly' : rpoly,
                         'rangeOff' : os.path.join(offdir, 'range_top_%02d_%02d.off'%(ii+1,ii+2)),
                         'azimuthOff': os.path.join(offdir, 'azimuth_top_%02d_%02d.off'%(ii+1,ii+2))}


                ###For future - should account for azimuth and range misreg here .. ignoring for now.
                cropInfo = getCropMetadata(
                    cropMetadata, secondaryDate, swath, secBurst.burstNumber)
                azCarrPoly, dpoly = estimateNativeCarrierPolynomials(
                    secondary, secBurst, cropInfo, offset=-1.0 * offset)
                rdict['carrPoly'] = azCarrPoly
                rdict['doppPoly'] = dpoly

                outimg = resampSecondary(topBurst, secBurst, rdict, outname, (not inps.noflat))

                copyBurst = copy.deepcopy(topBurst)
                ut.adjustValidSampleLine(copyBurst)
                copyBurst.image.filename = outimg.filename 
                print('After: ', copyBurst.firstValidLine, copyBurst.numValidLines)
                topCoreg.bursts.append(copyBurst)
                #######################################################


                secBurst = secondary.bursts[jj+1]
                outname = os.path.join(outdir, 'burst_bot_%02d_%02d.slc'%(ii+1,ii+2))

                ####Setup initial polynomials
                ### If no misregs are given, these are zero
                ### If provided, can be used for resampling without running to geo2rdr again for fast results
                rdict = {'azpoly' : apoly,
                         'rgpoly' : rpoly,
                         'rangeOff' : os.path.join(offdir, 'range_bot_%02d_%02d.off'%(ii+1,ii+2)),
                         'azimuthOff': os.path.join(offdir, 'azimuth_bot_%02d_%02d.off'%(ii+1,ii+2))}
                cropInfo = getCropMetadata(
                    cropMetadata, secondaryDate, swath, secBurst.burstNumber)
                azCarrPoly, dpoly = estimateNativeCarrierPolynomials(
                    secondary, secBurst, cropInfo, offset=-1.0 * offset)
                rdict['carrPoly'] = azCarrPoly
                rdict['doppPoly'] = dpoly

                outimg = resampSecondary(botBurst, secBurst, rdict, outname, (not inps.noflat))

                copyBurst = copy.deepcopy(botBurst)
                ut.adjustValidSampleLine(copyBurst)
                copyBurst.image.filename = outimg.filename
                print('After: ', copyBurst.firstValidLine, copyBurst.numValidLines)
                botCoreg.bursts.append(copyBurst)
               #######################################################

            else:
                outname = os.path.join(outdir, 'burst_%02d.slc'%(ii+1))  

                ####Setup initial polynomials
                ### If no misregs are given, these are zero
                ### If provided, can be used for resampling without running to geo2rdr again for fast results
                rdict = {'azpoly' : apoly,
                         'rgpoly' : rpoly,
                         'rangeOff' : os.path.join(offdir, 'range_%02d.off'%(ii+1)),
                         'azimuthOff': os.path.join(offdir, 'azimuth_%02d.off'%(ii+1))}


                ###For future - should account for azimuth and range misreg here .. ignoring for now.
                cropInfo = getCropMetadata(
                    cropMetadata, secondaryDate, swath, secBurst.burstNumber)
                azCarrPoly, dpoly = estimateNativeCarrierPolynomials(
                    secondary, secBurst, cropInfo, offset=-1.0 * offset)
                rdict['carrPoly'] = azCarrPoly
                rdict['doppPoly'] = dpoly

                outimg = resampSecondary(topBurst, secBurst, rdict, outname, (not inps.noflat))
                minAz, maxAz, minRg, maxRg = ut.getValidLines(secBurst, rdict, outname,
                    misreg_az = misreg_az - offset, misreg_rng = misreg_rg)


                copyBurst = copy.deepcopy(topBurst)
                ut.adjustValidSampleLine_V2(copyBurst, secBurst, minAz=minAz, maxAz=maxAz, minRng=minRg, maxRng=maxRg)
                copyBurst.image.filename = outimg.filename
                print('After: ', copyBurst.firstValidLine, copyBurst.numValidLines)
                topCoreg.bursts.append(copyBurst)


        ####################################################### 
        topCoreg.numberOfBursts = len(topCoreg.bursts)
        topCoreg.source = ut.asBaseClass(secondary)

        if inps.overlap:
            botCoreg.numberOfBursts = len(botCoreg.bursts)
            topCoreg.reference = ut.asBaseClass(referenceTop)
            botCoreg.reference = ut.asBaseClass(referenceBottom)
            botCoreg.source = ut.asBaseClass(secondary)
            ut.saveProduct(topCoreg, outdir + '_top.xml')
            ut.saveProduct(botCoreg, outdir + '_bottom.xml')

        else:
            topCoreg.reference = reference
            ut.saveProduct(topCoreg, outdir + '.xml')    


if __name__ == '__main__':
    '''
    Main driver.
    '''
    # Main Driver
    main()
