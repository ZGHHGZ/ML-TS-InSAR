#####直接在当前文件所在路径下运行，第2行的project为项目文件夹名称,第三行为分辨率
from scipy import stats

project="ltk"
dis=15 #输出分辨率，单位m。
import glob
from numba import jit
import numpy as np
import gdal2numpy
from osgeo import gdal, osr, ogr
import matplotlib.pyplot as plt
def create_empty_tif(output_file, xmin, ymin, xmax, ymax, pixel_size_x, pixel_size_y):
    # 计算图像的宽度和高度
    width = int((xmax - xmin) / pixel_size_x)
    height = int((ymax - ymin) / pixel_size_y)
    # 创建一个新的tif文件
    driver = gdal.GetDriverByName('GTiff')
    dataset = driver.Create(output_file, width, height, 1, gdal.GDT_Byte)
    # 设置地理变换参数
    geotransform = (xmin, pixel_size_x, 0, ymax, 0, -pixel_size_y)
    dataset.SetGeoTransform(geotransform)
    # 设置投影信息（这里使用WGS84坐标系）
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)  # WGS84
    dataset.SetProjection(srs.ExportToWkt())
    # 关闭数据集
    dataset.FlushCache()
    dataset = None
##########################################
lat = '../'+project+'/merged/geom_reference/lat.tif'
lon = '../'+project+'/merged/geom_reference/lon.tif'
lat_num = gdal.Open(lat).ReadAsArray()
lon_num = gdal.Open(lon).ReadAsArray()
lat_list=[lat_num[0][-1],lat_num[-1][0],lat_num[-1][-1],lat_num[0][0]]
lon_list=[lon_num[0][-1],lon_num[-1][0],lon_num[-1][-1],lon_num[0][0]]
#####±180度经度处可能出现问题
lat_min=min(lat_list)
lat_max=max(lat_list)
lon_min=min(lon_list)
lon_max=max(lon_list)
box="'"+str(lat_min-0.001)+" "+str(lat_max+0.001)+" "+str(lon_min-0.001)+" "+str(lon_max+0.001)+"'"
print(box)

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
#######################################################################################
xmin=lon_min-0.001
ymin=lat_min-0.001
xmax=lon_max+0.001
ymax=lat_max+0.001

# 示例：创建一个空的tif图像
create_empty_tif("../"+project+'/sample.tif', xmin, ymin, xmax, ymax,x,y )
###############################################################################

# 输入栅格文件路径
raster_file = "../"+project+'/sample.tif'

# 打开栅格数据集
dataset = gdal.Open(raster_file)

# 获取地理变换参数
transform = dataset.GetGeoTransform()

# 获取栅格数据的宽度和高度
width = dataset.RasterXSize
height = dataset.RasterYSize

# 创建两个新的栅格数据集，用于存储经度和纬度数据
lon_dataset = gdal.GetDriverByName('GTiff').Create("../"+project+'/lon_geo.tif', width, height, 1, gdal.GDT_Float32)
lat_dataset = gdal.GetDriverByName('GTiff').Create("../"+project+'/lat_geo.tif', width, height, 1, gdal.GDT_Float32)

# 设置新数据集的地理变换参数和投影信息
lon_dataset.SetGeoTransform(transform)
lat_dataset.SetGeoTransform(transform)
srs = ogr.osr.SpatialReference()
srs.ImportFromEPSG(4326)  # WGS84坐标系
lon_dataset.SetProjection(srs.ExportToWkt())
lat_dataset.SetProjection(srs.ExportToWkt())

# 遍历每个像素点，计算经纬度并写入对应的数据集
lon_num_geo=np.zeros((height,width))
lat_num_geo=np.zeros((height,width))
print(transform)
for i in range(height):
    for j in range(width):
        lon_num_geo[i,j] = transform[0] + j * transform[1] + i * transform[2]
        lat_num_geo[i,j] = transform[3]+ j * transform[4] + i * transform[5]
lon_dataset.GetRasterBand(1).WriteArray(np.array(lon_num_geo), xoff=0, yoff=0)
lat_dataset.GetRasterBand(1).WriteArray(np.array(lat_num_geo), xoff=0, yoff=0)

# 关闭数据集

dataset = None
lon_dataset = None
lat_dataset = None
#######################################################################################################
@jit(nopython=True)
def look_up_table():
    index=[]
    for i in range(height):
        for j in range(width):
            distance=(lon_num_geo[i,j]-lon_num)**2+(lat_num_geo[i,j]-lat_num)**2
        ###获取最小值的索引
            if np.min(distance)<y*y+x*x:
                index.append(np.argmin(distance))
            else:
                index.append(-1)
        print("line: "+str(i))
    return index
index=look_up_table()



url="../"+project+"/merged/interferograms/*/filt_fine.int"

data=glob.glob(url)
print(data)
sample=gdal2numpy.GDAL2Numpy("../"+project+'/sample.tif')
for i in range(len(data)):
    wrap = gdal.Open(data[i]).ReadAsArray()
    wrap=np.angle(wrap)
    ###扁平化
    wrap=wrap.flatten()
    wrap_geo=wrap[index]
    wrap_geo=wrap_geo.reshape(height,width)
    ##找众数
    mode_value = wrap_geo[len(wrap_geo)//2,0]
    # 将众数替换为nan
    wrap_geo[wrap_geo == mode_value] = np.nan
    gdal2numpy.Numpy2GTiff(wrap_geo,sample[1],sample[2],"../"+project+"/merged/interferograms/"+data[i].split("/")[-2]+"filt_fine_int.tif")
    print("../"+project+"/merged/interferograms/"+url.split("/")[-2]+"filt_fine_int.tif")
    #plt.imshow(wrap_geo.reshape(height,width),cmap='jet')
    #plt.colorbar()
    #plt.show()
