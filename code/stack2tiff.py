import glob
import os
import multiprocessing
from osgeo import gdal
cor=glob.glob(str(os.getcwd())+'/interferograms/*/filt_fine.cor')
unw=glob.glob('./interferograms/*/filt_fine.unw')
con=glob.glob('./interferograms/*/filt_fine.unw.conncomp')
intt=glob.glob('./interferograms/*/filt_fine.int')
lat='./geom_reference/lat.rdr'
lon='./geom_reference/lon.rdr'
hgt='./geom_reference/hgt.rdr'
los='./geom_reference/los.rdr'

unw_outFile='./interferograms/*/geo_filt_fine.unw'   
cmd2='rm -rf ' + str(unw_outFile)
cmd3='rm -rf ' + str(unw_outFile) + '.vrt'
cmd4='rm -rf ' + str(unw_outFile) + '.xml'
print('rm -rf ' + str(unw_outFile))
os.system(cmd2)
os.system(cmd3)
os.system(cmd4)   

unw_outFile='./interferograms/*/geo_filt_fine.cor'   
cmd2='rm -rf ' + str(unw_outFile)
cmd3='rm -rf ' + str(unw_outFile) + '.vrt'
cmd4='rm -rf ' + str(unw_outFile) + '.xml'
print('rm -rf ' + str(unw_outFile))
os.system(cmd2)
os.system(cmd3)
os.system(cmd4) 

unw_outFile='./interferograms/*/geo_filt_fine.unw.conncomp'   
cmd2='rm -rf ' + str(unw_outFile)
cmd3='rm -rf ' + str(unw_outFile) + '.vrt'
cmd4='rm -rf ' + str(unw_outFile) + '.xml'
print('rm -rf ' + str(unw_outFile))
os.system(cmd2)
os.system(cmd3)
os.system(cmd4)    

unw_outFile='./interferograms/*/geo_filt_fine.int'   
cmd2='rm -rf ' + str(unw_outFile)
cmd3='rm -rf ' + str(unw_outFile) + '.vrt'
cmd4='rm -rf ' + str(unw_outFile) + '.xml'
print('rm -rf ' + str(unw_outFile))
os.system(cmd2)
os.system(cmd3)
os.system(cmd4)    

unw_outFile='./geom_reference/geo_hgt.rdr'   
cmd2='rm -rf ' + str(unw_outFile)
cmd3='rm -rf ' + str(unw_outFile) + '.vrt'
cmd4='rm -rf ' + str(unw_outFile) + '.xml'
print('rm -rf ' + str(unw_outFile))
os.system(cmd2)
os.system(cmd3)
os.system(cmd4)   

unw_outFile='./geom_reference/geo_los.rdr'   
cmd2='rm -rf ' + str(unw_outFile)
cmd3='rm -rf ' + str(unw_outFile) + '.vrt'
cmd4='rm -rf ' + str(unw_outFile) + '.xml'
print('rm -rf ' + str(unw_outFile))
os.system(cmd2)
os.system(cmd3)
os.system(cmd4)  

########################################################################################################

os.system('cp ../../code/topsStack/geocodeGdal.py geocodeGdal.py')

########################################################################################################
def prepare_lat_lon(latFile,lonFile,file):
    width, length =  getSize(latFile)
    widthFile , lengthFile = getSize(file)

    print("size of lat and lon files (width, length) ", width, length)
    print("size of input file to be geocoded (width, length): ", widthFile , lengthFile)

    xOff = 0
    yOff = 0

    cmd = 'gdal_translate -of VRT -srcwin ' + str(xOff) + ' ' + str(yOff) \
           +' '+ str(width - xOff) +' '+ str(length - yOff) +' -outsize ' + str(widthFile) + \
           ' '+ str(lengthFile)  + ' -a_nodata 0 ' + latFile +'.vrt' + ' tempLAT.vrt'

    os.system(cmd)

    cmd = 'gdal_translate -of VRT -srcwin ' + str(xOff) + ' ' + str(yOff) \
          +' '+ str(int(width-xOff)) +' '+ str(int(length-yOff)) +' -outsize ' + str(widthFile) +\
           ' '+ str(lengthFile)  + ' -a_nodata 0 ' + lonFile +'.vrt' + ' tempLON.vrt'

    os.system(cmd)

    return 'tempLAT.vrt', 'tempLON.vrt'



def getSize(infile):    

    ds=gdal.Open(infile + ".vrt")
    b=ds.GetRasterBand(1)
    return b.XSize, b.YSize


########################################################################################################

import gdal2numpy
import numpy as np

os.system('gdal_translate ./geom_reference/lat.rdr ./geom_reference/lat.tif -of GTiff')
os.system('gdal_translate ./geom_reference/lon.rdr ./geom_reference/lon.tif -of GTiff')


lat = './geom_reference/lat.tif'
lon = './geom_reference/lon.tif'
lat_num = gdal2numpy.GDAL2Numpy(lat)[0]
lon_num = gdal2numpy.GDAL2Numpy(lon)[0]
lat_list=[lat_num[0][-1],lat_num[-1][0],lat_num[-1][-1],lat_num[0][0]]
lon_list=[lon_num[0][-1],lon_num[-1][0],lon_num[-1][-1],lon_num[0][0]]
#####±180度经度处可能出现问题
lat_min=min(lat_list)
lat_max=max(lat_list)
lon_min=min(lon_list)
lon_max=max(lon_list)
box="'"+str(lat_min-0.001)+" "+str(lat_max+0.001)+" "+str(lon_min-0.001)+" "+str(lon_max+0.001)+"'"
print(box)

dis=15 #输出分辨率

ave_lat=(float(lat_min)+float(lat_max))/2

ave_lon=(float(lon_min)+float(lon_max))/2

from geopy import distance               #####安装geopy

lat_1=(ave_lat, ave_lon)
lat_2=(ave_lat+0.0001, ave_lon)

lon_1=(ave_lat, ave_lon)
lon_2=(ave_lat, ave_lon+0.0001)

distance_lat=distance.geodesic(lat_1, lat_2).m
distance_lon=distance.geodesic(lon_1, lon_2).m
print(distance_lat)
print(distance_lon)
x=dis/distance_lat*0.0001
y=dis/distance_lon*0.0001






#############################################################################################################
hgt='./geom_reference/hgt.rdr'
los='./geom_reference/los.rdr'
latFile, lonFile = prepare_lat_lon('./geom_reference/lat.rdr','./geom_reference/lon.rdr',str(hgt))
os.system("python ./geocodeGdal.py -l ./geom_reference/lat.rdr -L ./geom_reference/lon.rdr -f "+str(hgt)+" --bbox "+box+" -x "+str(x)+" -y "+str(y)+" -t")
latFile, lonFile = prepare_lat_lon('./geom_reference/lat.rdr','./geom_reference/lon.rdr',str(los))
os.system("python ./geocodeGdal.py -l ./geom_reference/lat.rdr -L ./geom_reference/lon.rdr -f "+str(los)+" --bbox "+box+" -x "+str(x)+" -y "+str(y)+" -t")


unw_outFile='./geom_reference/geo_hgt.rdr'   
cmd2='rm -rf ' + str(unw_outFile)
cmd3='rm -rf ' + str(unw_outFile) + '.vrt'
cmd4='rm -rf ' + str(unw_outFile) + '.xml'
print('rm -rf ' + str(unw_outFile))
os.system(cmd2)
os.system(cmd3)
os.system(cmd4)   

unw_outFile='./geom_reference/los_geo.rdr'   
cmd2='rm -rf ' + str(unw_outFile)
cmd3='rm -rf ' + str(unw_outFile) + '.vrt'
cmd4='rm -rf ' + str(unw_outFile) + '.xml'
print('rm -rf ' + str(unw_outFile))
os.system(cmd2)
os.system(cmd3)
os.system(cmd4)   

#############################################################################################################




def unw_tif(unw):
    os.system("python ./geocodeGdal.py -l ./geom_reference/lat.rdr -L ./geom_reference/lon.rdr -f "+str(unw)+" --bbox "+box+" -x "+str(x)+" -y "+str(y)+" -t") 
        
latFile, lonFile = prepare_lat_lon('./geom_reference/lat.rdr','./geom_reference/lon.rdr',str(unw[0]))
pool = multiprocessing.Pool(8)
results = pool.map(unw_tif, unw)    
pool.close()
pool.join()    
    
unw_outFile='./interferograms/*/geo_filt_fine.unw'   
cmd2='rm -rf ' + str(unw_outFile)
cmd3='rm -rf ' + str(unw_outFile) + '.vrt'
cmd4='rm -rf ' + str(unw_outFile) + '.xml'
print('rm -rf ' + str(unw_outFile))
os.system(cmd2)
os.system(cmd3)
os.system(cmd4)    
    
    
    
def cor_tif(cor):    
    os.system("python ./geocodeGdal.py -l ./geom_reference/lat.rdr -L ./geom_reference/lon.rdr -f "+str(cor)+" --bbox "+box+" -x "+str(x)+" -y "+str(y)+" -t")
    
latFile, lonFile = prepare_lat_lon('./geom_reference/lat.rdr','./geom_reference/lon.rdr',str(cor[0]))
pool = multiprocessing.Pool(8)
results = pool.map(cor_tif, cor)  
pool.close()
pool.join()

unw_outFile='./interferograms/*/geo_filt_fine.cor'   
cmd2='rm -rf ' + str(unw_outFile)
cmd3='rm -rf ' + str(unw_outFile) + '.vrt'
cmd4='rm -rf ' + str(unw_outFile) + '.xml'
print('rm -rf ' + str(unw_outFile))
os.system(cmd2)
os.system(cmd3)
os.system(cmd4)    
    
  
    
def con_tif(con):        
    os.system("python ../../code/topsStack/geocodeGdal.py -l ./geom_reference/lat.rdr -L ./geom_reference/lon.rdr -f "+str(con)+" --bbox "+box+" -x "+str(x)+" -y "+str(y)+" -t")

latFile, lonFile = prepare_lat_lon('./geom_reference/lat.rdr','./geom_reference/lon.rdr',str(con[0]))        
pool = multiprocessing.Pool(8)
results = pool.map(con_tif, con)        
pool.close()
pool.join() 
    
unw_outFile='./interferograms/*/geo_filt_fine.unw.conncomp'   
cmd2='rm -rf ' + str(unw_outFile)
cmd3='rm -rf ' + str(unw_outFile) + '.vrt'
cmd4='rm -rf ' + str(unw_outFile) + '.xml'
print('rm -rf ' + str(unw_outFile))
os.system(cmd2)
os.system(cmd3)
os.system(cmd4)    
    
    
    
def intt_tif(intt):        
    os.system("python ./geocodeGdal.py -l ./geom_reference/lat.rdr -L ./geom_reference/lon.rdr -f "+str(intt)+" --bbox "+box+" -x "+str(x)+" -y "+str(y)+" -t")

latFile, lonFile = prepare_lat_lon('./geom_reference/lat.rdr','./geom_reference/lon.rdr',str(intt[0]))
pool = multiprocessing.Pool(8)
results = pool.map(intt_tif, intt)        
pool.close()
pool.join()

unw_outFile='./interferograms/*/geo_filt_fine.int'   
cmd2='rm -rf ' + str(unw_outFile)
cmd3='rm -rf ' + str(unw_outFile) + '.vrt'
cmd4='rm -rf ' + str(unw_outFile) + '.xml'
print('rm -rf ' + str(unw_outFile))
os.system(cmd2)
os.system(cmd3)
os.system(cmd4)    






