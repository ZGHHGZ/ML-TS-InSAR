#########在merged文件夹下运行，用于TIFF转HYP3格式#########
import os
import glob

cor=sorted(glob.glob('./interferograms/*/geo_filt_fine.cor.tif'))
unw=sorted(glob.glob('./interferograms/*/geo_filt_fine.unw.tif'))
con=sorted(glob.glob('./interferograms/*/geo_filt_fine.unw.conncomp.tif'))
intt=sorted(glob.glob('./interferograms/*/geo_filt_fine.int.tif'))



#/home/jovyan/iscepredata/p/merged/interferograms/20200810_20200828/geo_filt_fine.cor.tif

os.system("rm -rf tiff_data")
os.system("mkdir tiff_data")
os.system("mkdir ./tiff_data/int")

for i in range(len(cor)):
    os.system("cp "+str(cor[i])+" ./tiff_data/int/S1_000000_IW0_"+str(cor[i][-39:-22])+"_VV_INT00_0000_corr_clipped.tif")                    
    print("cp "+str(cor[i])+" ./tiff_data/int/S1_000000_IW0_"+str(cor[i][-39:-22])+"_VV_INT00_0000_corr_clipped.tif")
    
for i in range(len(unw)):
    os.system("cp "+str(unw[i])+" ./tiff_data/int/S1_000000_IW0_"+str(unw[i][-39:-22])+"_VV_INT00_0000_unw_phase_clipped.tif")                    
    print("cp "+str(unw[i])+" ./tiff_data/int/S1_000000_IW0_"+str(unw[i][-39:-22])+"_VV_INT00_0000_unw_phase_clipped.tif")
      
hgt='./geom_reference/geo_hgt.rdr.tif'
phi='./geom_reference/geo_los.rdr_phi.tif'
theta='./geom_reference/geo_los.rdr_theta.tif'
os.system("cp "+str(hgt)+" ./tiff_data/int/S1_000000_IW0_"+str(cor[0][-39:-22])+"_VV_INT00_0000_dem_clipped.tif")

######################################################################################################################

import gdal2numpy
import numpy as np

phi = gdal2numpy.GDAL2Numpy(phi)
phi_num = (90 + phi[0]) * np.pi / 180
gdal2numpy.Numpy2GTiff(phi_num, phi[1], phi[2], "./tiff_data/int/S1_000000_IW0_" + str(cor[0][-39:-22]) + "_VV_INT00_0000_lv_phi_clipped.tif")

theta = gdal2numpy.GDAL2Numpy(theta)
theta_num = (90 - theta[0]) * np.pi / 180
gdal2numpy.Numpy2GTiff(theta_num, theta[1], theta[2], "./tiff_data/int/S1_000000_IW0_" + str(cor[0][-39:-22]) + "_VV_INT00_0000_lv_theta_clipped.tif")


#water_mask (nodata_mask)
water=theta_num
water[np.isnan(water)] = 0
water[water != 0] = 1
gdal2numpy.Numpy2GTiff(water, theta[1], theta[2], "./tiff_data/int/S1_000000_IW0_" + str(cor[0][-39:-22]) + "_VV_INT00_0000_water_mask_clipped.tif")

######################################################################################################################


for i in range(len(cor)):
    
    baselines1_path=glob.glob('../baselines/*'+str(cor[i][-39:-31])+'/*'+str(cor[i][-39:-31])+'.txt')[0]
    baselines1 = open(baselines1_path).readlines()
    baselines2_path=glob.glob('../baselines/*'+str(cor[i][-30:-22])+'/*'+str(cor[i][-30:-22])+'.txt')[0]
    baselines2 = open(baselines2_path).readlines()
    
    baseline1=baselines1[1][17:-1]
    baseline2=baselines2[1][17:-1]
    
    if baseline1=="nan":
        baseline=str(float(baseline2))
    elif baseline2=="nan":
        baseline=str(-float(baseline1))
    else:
        baseline=str(float(baseline2)-float(baseline1))
    
    
    ref='S1_000000_IW0_'+str(cor[i][-39:-31])+'T000000_VV_0000-BURST'
    sec='S1_000000_IW0_'+str(cor[i][-30:-22])+'T000000_VV_0000-BURST'
    lines=[]
    lines.append("Reference Granule:"+str(ref)+"\n")
    lines.append("Secondary Granule:"+str(sec)+"\n")
    lines.append("Reference Pass Direction:ASCENDING\n")
    lines.append("Reference Orbit Number:0\n")
    lines.append("Secondary Pass Direction:ASCENDING\n")
    lines.append("Secondary Orbit Number:0\n")
    lines.append("Baseline:"+baseline+"\n")
    lines.append("UTC time:0\n")
    lines.append("Heading:0\n")
    lines.append("Spacecraft height:693000.0\n")
    lines.append("Earth radius at nadir:6337286.638938101\n")
    lines.append("Slant range near:0\n")
    lines.append("Slant range center:0\n")
    lines.append("Slant range far:0\n")
    lines.append("Range looks:0\n")
    lines.append("Azimuth looks:0\n")
    lines.append("INSAR phase filter:yes\n")
    lines.append("Phase filter parameter:0.5\n")
    lines.append("Range bandpass filter:False\n")
    lines.append("Azimuth bandpass filter:False\n")
    lines.append("DEM source:GLO_30\n")
    lines.append("DEM resolution (m):30\n")
    lines.append("Unwrapping type:snaphu_mcf\n")
    lines.append("Speckle filter:True\n")
    lines.append("Water mask:no\n")
    lines.append("Radar n lines:0\n")
    lines.append("Radar n samples:0\n")
    lines.append("Radar first valid line:0\n")
    lines.append("Radar n valid lines:0\n")
    lines.append("Radar first valid sample:0\n")
    lines.append("Radar n valid samples:0\n")
    lines.append("Multilook azimuth time interval:0\n")
    lines.append("Multilook range pixel size:0\n")
    lines.append("Radar sensing stop:0\n")
    txt_path='./tiff_data/int/S1_000000_IW0_'+str(cor[i][-39:-22])+'_VV_INT00_0000.txt'
    print(txt_path)
    open(txt_path,'w').writelines(lines)
#########################################################################################################################
