from coregSwathSLCProduct import coregSwathSLCProduct
import isce
import isceobj
import os
import numpy as np
import json
import datetime
#from isceobj.Sensor.TOPS.coregSwathSLCProduct import coregSwathSLCProduct

class catalog(object):
    def __init__(self):
        pass

    def addItem(self,*args):
        print(' '.join([str(x) for x in args]))
          


def loadProduct(xmlname):
        '''
        Load the product using Product Manager.
        '''

        from iscesys.Component.ProductManager import ProductManager as PM
        pm = PM()
        pm.configure()
        obj = pm.loadProduct(xmlname)

        return obj


def saveProduct( obj, xmlname):
    '''
    Save the product to an XML file using Product Manager.
    '''
    import shelve
    import os
    with shelve.open(os.path.dirname(xmlname) + '/'+ os.path.basename(xmlname)  +'.data') as db:
             db['data'] = obj

    from iscesys.Component.ProductManager import ProductManager as PM

    pm = PM()
    pm.configure()

    pm.dumpProduct(obj, xmlname)

    return None


def getCommonBurstLimits(reference, secondary):
    """Return common burst indices, robust to a cropped reference grid.

    ISCE's native getBurstOffset() compares burst timing/geometry and may fail
    after sensingStart and startingRange are shifted to describe a cropped
    reference.  The burstNumber field is not changed by spatial cropping, so
    use it as a deterministic fallback.  Returned values keep ISCE's original
    convention: (secondary_index - reference_index, reference_start,
    reference_stop_exclusive).
    """
    try:
        return reference.getCommonBurstLimits(secondary)
    except Exception as native_error:
        ref_index = {}
        sec_index = {}
        for index, burst in enumerate(reference.bursts):
            number = int(burst.burstNumber)
            if number in ref_index:
                raise RuntimeError(
                    'Duplicate reference burstNumber: {0}'.format(number))
            ref_index[number] = index
        for index, burst in enumerate(secondary.bursts):
            number = int(burst.burstNumber)
            if number in sec_index:
                raise RuntimeError(
                    'Duplicate secondary burstNumber: {0}'.format(number))
            sec_index[number] = index

        common_numbers = sorted(set(ref_index).intersection(sec_index))
        if not common_numbers:
            raise RuntimeError(
                'Native burst matching failed ({0}); reference and secondary '
                'have no common burstNumber'.format(native_error))

        ref_indices = [ref_index[number] for number in common_numbers]
        sec_indices = [sec_index[number] for number in common_numbers]
        expected_ref = list(range(ref_indices[0], ref_indices[-1] + 1))
        expected_sec = list(range(sec_indices[0], sec_indices[-1] + 1))
        if ref_indices != expected_ref or sec_indices != expected_sec:
            raise RuntimeError(
                'Common burstNumber sequence is not contiguous: {0}'.format(
                    common_numbers))
        if len(ref_indices) != len(sec_indices):
            raise RuntimeError('Reference/secondary common burst count differs')

        burst_offset = sec_indices[0] - ref_indices[0]
        print('Warning: native getCommonBurstLimits failed: {0}'.format(
            native_error))
        print('Using burstNumber fallback: numbers={0}, reference=[{1},{2}), '
              'secondary_start={3}, offset={4}'.format(
                  common_numbers, ref_indices[0], ref_indices[-1] + 1,
                  sec_indices[0], burst_offset))
        return burst_offset, ref_indices[0], ref_indices[-1] + 1


def getRelativeShifts(mFrame, sFrame, minBurst, maxBurst, secondaryBurstStart):
    '''
    Estimate the relative shifts between the start of the bursts.
    '''
    import numpy as np    
    azReferenceOff = {}
    azSecondaryOff = {}
    azRelOff = {}
    tm = mFrame.bursts[minBurst].sensingStart
    dt = mFrame.bursts[minBurst].azimuthTimeInterval
    ts = sFrame.bursts[secondaryBurstStart].sensingStart
    
    for index in range(minBurst, maxBurst):
        burst = mFrame.bursts[index]
        azReferenceOff[index] = int(np.round((burst.sensingStart - tm).total_seconds() / dt))
        
        burst = sFrame.bursts[secondaryBurstStart + index - minBurst]
        azSecondaryOff[secondaryBurstStart + index - minBurst] =  int(np.round((burst.sensingStart - ts).total_seconds() / dt))
        
        azRelOff[secondaryBurstStart + index - minBurst] = azSecondaryOff[secondaryBurstStart + index - minBurst] - azReferenceOff[index]

    
    return azRelOff


def getRelativeShiftsFromReferenceNativeGrid(
        mFrame, sFrame, minBurst, maxBurst, secondaryBurstStart,
        swath, metadataPath=os.path.join('geom_reference',
                                        'crop_metadata.json')):
    """Compute burst shifts using the reference grid before spatial crop.

    Cropping changes each reference burst sensingStart by a different yoff.
    Using those cropped times in getRelativeShifts() incorrectly injects that
    yoff difference into the secondary azimuth-carrier polynomial.  Restore
    the per-burst native sensingStart saved by pre_data_mburst.py instead.
    """
    if not os.path.exists(metadataPath):
        print('Native reference crop metadata not found; using product times')
        return getRelativeShifts(
            mFrame, sFrame, minBurst, maxBurst, secondaryBurstStart)

    with open(metadataPath, 'r', encoding='utf-8') as file:
        payload = json.load(file)
    iw = 'IW{0}'.format(swath)
    # version 3 records the grid before either reference crop. Keep the
    # version-2 fallback so existing geometry can resume at resampling.
    referenceInfo = payload.get('reference_native_swaths', {}).get(iw)
    metadataGrid = 'native'
    if not referenceInfo:
        referenceInfo = payload.get('reference_swaths', {}).get(iw)
        metadataGrid = 'legacy'
    if not referenceInfo:
        raise RuntimeError(
            'Missing native reference timing metadata for {0}'.format(iw))

    dt = float(mFrame.bursts[minBurst].azimuthTimeInterval)
    firstKey = '{0:02d}'.format(minBurst + 1)
    if firstKey not in referenceInfo:
        raise RuntimeError(
            'Missing native reference burst metadata: {0} {1}'.format(
                iw, firstKey))
    nativeReferenceStart = datetime.datetime.fromisoformat(
        str(referenceInfo[firstKey]['original_sensing_start']))
    secondaryStart = sFrame.bursts[secondaryBurstStart].sensingStart

    shifts = {}
    for referenceIndex in range(minBurst, maxBurst):
        key = '{0:02d}'.format(referenceIndex + 1)
        info = referenceInfo.get(key)
        if info is None:
            raise RuntimeError(
                'Missing native reference burst metadata: {0} {1}'.format(
                    iw, key))
        referenceTime = datetime.datetime.fromisoformat(
            str(info['original_sensing_start']))
        referenceOffset = int(np.round(
            (referenceTime - nativeReferenceStart).total_seconds() / dt))

        secondaryIndex = (secondaryBurstStart +
                          referenceIndex - minBurst)
        secondaryTime = sFrame.bursts[secondaryIndex].sensingStart
        secondaryOffset = int(np.round(
            (secondaryTime - secondaryStart).total_seconds() / dt))
        shifts[secondaryIndex] = secondaryOffset - referenceOffset

    print('Native-grid relative burst shifts for {0} ({1} metadata): {2}'.format(
        iw, metadataGrid, shifts))
    return shifts


def adjustValidSampleLine(reference,  minAz=0, maxAz=0, minRng=0, maxRng=0):
    import numpy as np
    import isce
    import isceobj
    # Valid region in the resampled slc based on offsets
    ####Adjust valid samples and first valid sample here
    print ("Adjust valid samples")
    print('Before: ', reference.firstValidSample, reference.numValidSamples)
    print('Offsets : ', minRng, maxRng)
    if (minRng > 0) and (maxRng > 0):
            reference.numValidSamples -= (int(np.ceil(maxRng)) + 8)
            reference.firstValidSample += 4
    elif (minRng < 0) and  (maxRng < 0):
            reference.firstValidSample -= int(np.floor(minRng) - 4)
            reference.numValidSamples += int(np.floor(minRng) - 8)
    elif (minRng < 0) and (maxRng > 0):
            reference.firstValidSample -= int(np.floor(minRng) - 4)
            reference.numValidSamples += int(np.floor(minRng) - 8) - int(np.ceil(maxRng))

    print('After: ', reference.firstValidSample, reference.numValidSamples)

    ###Adjust valid lines and first valid line here
    print ("Adjust valid lines")
    print('Before: ', reference.firstValidLine, reference.numValidLines)
    print('Offsets : ', minAz, maxAz)
    if (minAz > 0) and (maxAz > 0):
            reference.numValidLines -= (int(np.ceil(maxAz)) + 8)
            reference.firstValidLine += 4
    elif (minAz < 0) and  (maxAz < 0):
            reference.firstValidLine -= int(np.floor(minAz) - 4)
            reference.numValidLines += int(np.floor(minAz) - 8)
    elif (minAz < 0) and (maxAz > 0):
            reference.firstValidLine -= int(np.floor(minAz) - 4)
            reference.numValidLines += int(np.floor(minAz) - 8) - int(np.ceil(maxAz))
    print('After:', reference.firstValidLine, reference.numValidLines)


def adjustValidSampleLine_V2(reference, secondary, minAz=0, maxAz=0, minRng=0, maxRng=0): 
    import numpy as np
    import isce
    import isceobj
    ####Adjust valid samples and first valid sample here
    print ("Adjust valid samples")
    print('Before: ', reference.firstValidSample, reference.numValidSamples)
    print('Offsets : ', minRng, maxRng)

    if (minRng > 0) and (maxRng > 0):
        reference.firstValidSample = secondary.firstValidSample - int(np.floor(maxRng)-4)
        lastValidSample = reference.firstValidSample - 8 + secondary.numValidSamples

        if lastValidSample < reference.numberOfSamples:
            reference.numValidSamples = secondary.numValidSamples - 8
        else:
            reference.numValidSamples = reference.numberOfSamples - reference.firstValidSample

    elif (minRng < 0) and (maxRng < 0):
        reference.firstValidSample = secondary.firstValidSample - int(np.floor(minRng) - 4)
        lastValidSample = reference.firstValidSample + secondary.numValidSamples  - 8
        if lastValidSample < reference.numberOfSamples:
            reference.numValidSamples = secondary.numValidSamples - 8
        else:
            reference.numValidSamples = reference.numberOfSamples - reference.firstValidSample
    elif (minRng < 0) and (maxRng > 0):
        reference.firstValidSample = secondary.firstValidSample - int(np.floor(minRng) - 4)
        lastValidSample = reference.firstValidSample + secondary.numValidSamples + int(np.floor(minRng) - 8) - int(np.ceil(maxRng))
        if lastValidSample < reference.numberOfSamples:
            reference.numValidSamples = secondary.numValidSamples + int(np.floor(minRng) - 8) - int(np.ceil(maxRng))
        else:
            reference.numValidSamples = reference.numberOfSamples - reference.firstValidSample

    reference.firstValidSample = np.max([0, reference.firstValidSample])
 
    print('After: ', reference.firstValidSample, reference.numValidSamples)

    ###Adjust valid lines and first valid line here
    print ("Adjust valid lines")
    print('Before: ', reference.firstValidLine, reference.numValidLines)
    print('Offsets : ', minAz, maxAz)

    if (minAz > 0) and (maxAz > 0):

        reference.firstValidLine = secondary.firstValidLine - int(np.floor(maxAz) - 4)
        lastValidLine = reference.firstValidLine - 8  + secondary.numValidLines
        if lastValidLine < reference.numberOfLines:
            reference.numValidLines = secondary.numValidLines - 8
        else:
            reference.numValidLines = reference.numberOfLines - reference.firstValidLine

    elif (minAz < 0) and  (maxAz < 0):
        reference.firstValidLine = secondary.firstValidLine - int(np.floor(minAz) - 4)
        lastValidLine = reference.firstValidLine + secondary.numValidLines  - 8
        if lastValidLine < reference.numberOfLines:
            reference.numValidLines = secondary.numValidLines - 8
        else:
            reference.numValidLines = reference.numberOfLines - reference.firstValidLine

    elif (minAz < 0) and (maxAz > 0):
        reference.firstValidLine = secondary.firstValidLine - int(np.floor(minAz) - 4)
        lastValidLine = reference.firstValidLine + secondary.numValidLines + int(np.floor(minAz) - 8) - int(np.ceil(maxAz))
        if lastValidLine < reference.numberOfLines:
            reference.numValidLines = secondary.numValidLines + int(np.floor(minAz) - 8) - int(np.ceil(maxAz))
        else:
            reference.numValidLines = reference.numberOfLines - reference.firstValidLine

    return reference


def adjustCommonValidRegion(reference,secondary):
    # valid lines between reference and secondary


    reference_lastValidLine = reference.firstValidLine + reference.numValidLines - 1
    reference_lastValidSample = reference.firstValidSample + reference.numValidSamples - 1
    secondary_lastValidLine = secondary.firstValidLine + secondary.numValidLines - 1
    secondary_lastValidSample = secondary.firstValidSample + secondary.numValidSamples - 1

    igram_lastValidLine = min(reference_lastValidLine, secondary_lastValidLine)
    igram_lastValidSample = min(reference_lastValidSample, secondary_lastValidSample)

    reference.firstValidLine = max(reference.firstValidLine, secondary.firstValidLine)
    reference.firstValidSample = max(reference.firstValidSample, secondary.firstValidSample)

    #set to 0 to avoid negative values
    if reference.firstValidLine<0:
        reference.firstValidLine=0
    if reference.firstValidSample<0:
        reference.firstValidSample=0

    reference.numValidLines = igram_lastValidLine - reference.firstValidLine + 1
    reference.numValidSamples = igram_lastValidSample - reference.firstValidSample + 1


def getValidLines(secondary, rdict, inname, misreg_az=0.0, misreg_rng=0.0):
    '''
    Looks at the reference, secondary and azimuth offsets and gets the Interferogram valid lines 
    '''
    import numpy as np
    import isce
    import isceobj

    dimg = isceobj.createSlcImage()
    dimg.load(inname + '.xml')
    shp = (dimg.length, dimg.width)
    az = np.fromfile(rdict['azimuthOff'], dtype=np.float32).reshape(shp)
    az += misreg_az
    aa = np.zeros(az.shape)
    aa[:,:] = az
    aa[aa < -10000.0] = np.nan
    amin = np.nanmin(aa)
    amax = np.nanmax(aa)

    rng = np.fromfile(rdict['rangeOff'], dtype=np.float32).reshape(shp)
    rng += misreg_rng
    rr = np.zeros(rng.shape)
    rr[:,:] = rng
    rr[rr < -10000.0] = np.nan
    rmin = np.nanmin(rr)
    rmax = np.nanmax(rr)

    return amin, amax, rmin, rmax


def adjustValidRegionFromOffsets(reference, secondary, rdict,
                                 kernelMargin=4):
    """Set the valid output rectangle from the actual geo2rdr mapping.

    This works when the reference is cropped but the secondary stays on its
    full native grid.  For each output pixel, geo2rdr provides the secondary
    coordinate as output_index + offset.  Pixels are valid only when that
    coordinate falls inside the secondary valid-data rectangle, with a small
    interpolation-kernel margin.
    """
    height = int(reference.numberOfLines)
    width = int(reference.numberOfSamples)
    azimuth = np.fromfile(
        rdict['azimuthOff'], dtype=np.float32).reshape(height, width)
    rng = np.fromfile(
        rdict['rangeOff'], dtype=np.float32).reshape(height, width)

    rows = np.arange(height, dtype=np.float64)[:, None]
    cols = np.arange(width, dtype=np.float64)[None, :]
    secondaryRows = rows + azimuth
    secondaryCols = cols + rng

    firstLine = int(secondary.firstValidLine) + int(kernelMargin)
    lastLine = (int(secondary.firstValidLine) +
                int(secondary.numValidLines) - 1 - int(kernelMargin))
    firstSample = int(secondary.firstValidSample) + int(kernelMargin)
    lastSample = (int(secondary.firstValidSample) +
                  int(secondary.numValidSamples) - 1 - int(kernelMargin))

    valid = (np.isfinite(azimuth) & np.isfinite(rng) &
             (azimuth > -10000.0) & (rng > -10000.0) &
             (secondaryRows >= firstLine) &
             (secondaryRows <= lastLine) &
             (secondaryCols >= firstSample) &
             (secondaryCols <= lastSample))
    validRows = np.flatnonzero(np.any(valid, axis=1))
    validCols = np.flatnonzero(np.any(valid, axis=0))
    if validRows.size == 0 or validCols.size == 0:
        raise RuntimeError(
            'Resampled burst has no pixels inside the secondary valid region')

    reference.firstValidLine = int(validRows[0])
    reference.numValidLines = int(validRows[-1] - validRows[0] + 1)
    reference.firstValidSample = int(validCols[0])
    reference.numValidSamples = int(validCols[-1] - validCols[0] + 1)
    print('Valid output from geo2rdr offsets: line={0}+{1}, sample={2}+{3}'.format(
        reference.firstValidLine, reference.numValidLines,
        reference.firstValidSample, reference.numValidSamples))
    return reference



def asBaseClass(inobj):
    '''
    Return as TOPSSwathSLCProduct.
    '''
    from isceobj.Sensor.TOPS.TOPSSwathSLCProduct import TOPSSwathSLCProduct
    
    
    def topsproduct(cobj):
        obj = TOPSSwathSLCProduct()
        obj.configure()

        for x in obj.parameter_list:
            val = getattr(cobj, x.attrname)
            setattr(obj, x.attrname, val)

        for x in obj.facility_list:
            attrname = x.public_name
            val = getattr(cobj, x.attrname)
            setattr(obj, x.attrname, val)
        
        return obj


    if isinstance(inobj, coregSwathSLCProduct):
        return topsproduct(inobj)

    elif isinstance(inobj, TOPSSwathSLCProduct):
        return inobj
    else:
        raise Exception('Cannot be converted to TOPSSwathSLCProduct')


def getSwathList(indir):

    swathList = []
    for x in [1,2,3]:
        SW = os.path.join(indir,'IW{0}'.format(x))
        if os.path.exists(SW):
            swathList.append(x)

    return swathList
