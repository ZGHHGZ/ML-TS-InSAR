####在code文件夹同级新建项目文件夹，并在项目文件夹内新建data文件夹，data文件夹内放入stack.txt文件###
####在data文件夹下运行此脚本##################################################################
#*****************************************************************************#
#*****************************************************************************#
########################      data download      ##############################
#*****************************************************************************#
#*****************************************************************************#

import os
import numpy
from pathlib import Path
from hyp3lib.get_orb import downloadSentinelOrbitFile
from burst import download_bursts, get_burst_params, get_isce2_burst_bbox, get_region_of_interest
from hyp3_isce2.dem import download_dem_for_isce2
from s1_orbits import fetch_for_scene
import warnings
import multiprocessing
warnings.filterwarnings("ignore")
orbit_dir = Path('orbits')
aux_cal_dir = Path('aux_cal')
dem_dir = Path('dem')
data=numpy.loadtxt('stack.txt',dtype=str,delimiter=',')

#影像数据下载（BURST）
def imgae_download(data):
    scene_name = data[0]
    swath_number = data[1]
    burst_number = data[2]
    params = get_burst_params(scene_name ,swath_number, int(burst_number))
    if not os.path.exists(f'./{params.granule}.SAFE'):
        metadata = download_bursts([params])
        if not os.path.exists("./" + str(params.granule) + '.SAFE/preview/map-overlay.kml'):
            footprint = get_isce2_burst_bbox(params)
            footprint=footprint.buffer(0.2).bounds
            #print("create————————————kml")
            cord_list=str(footprint)[1:-1].split(', ')
            print(cord_list)
            cord_kml=str(cord_list[0])+','+str(cord_list[1])+' '+str(cord_list[2])+','+str(cord_list[1])+' '+str(cord_list[2])+','+str(cord_list[3])+' '+str(cord_list[0])+','+str(cord_list[3])
            path = "./" + str(params.granule) + '.SAFE/preview'
            if not os.path.exists(path):
                os.mkdir(path)
            with open('../../code/map-overlay.txt', 'r') as f:
                kml=f.read()
                f.close()
            kml = kml.replace('          <coordinates>', '          <coordinates>' + cord_kml)  # 在第一个</coordinates>前插入
            with open(path+'/map-overlay.kml', 'w') as f:
                f.write(kml)
            #print(str(params.granule)+"------kml创建成功")
            print("image: " + str(params.granule))
    else:
        print("skip: " + str(params.granule))
################################################################################
pool = multiprocessing.Pool(64)
results = pool.map(imgae_download, data)
pool.close()
pool.join()
################################################################################

#轨道数据下载    
download_aux_cal(aux_cal_dir)

orbit_dir.mkdir(exist_ok=True, parents=True)
def orbit_download(data):
    scene_name = data[0]
    fetch_for_scene(scene_name, str(orbit_dir))
    print('orbit: ',data[0])
################################################################################

pool = multiprocessing.Pool(64)
results = pool.map(orbit_download, data)
pool.close()
pool.join()
################################################################################

#DEM数据下载    
scene_name = data[0][0]
swath_number = data[0][1]
burst_number = data[0][2]
params = get_burst_params(scene_name ,swath_number, int(burst_number))    
dem_roi = get_isce2_burst_bbox(params).bounds
print(dem_roi)
dem_path = download_dem_for_isce2(dem_roi,dem_name='glo_30', dem_dir=dem_dir, buffer=0.2)


#*****************************************************************************#
#*****************************************************************************#
########################      isce process      ###############################
#*****************************************************************************#
#*****************************************************************************#
import os
import numpy
import isce 
import shutil
from pathlib import Path
from burst import download_bursts, get_burst_params, get_isce2_burst_bbox, get_region_of_interest
from topsStack import stackSentinel
from osgeo import gdal
import glob

os.chdir('../')
path= "./run_files"
if os.path.exists(path):
    print("run_files were existed")


else:
    os.chdir('./data')
    data=numpy.loadtxt('stack.txt',dtype=str,delimiter=',')
    ref_params = get_burst_params(data[0][0],data[0][1],int(data[0][2]))
    sec_params = get_burst_params(data[1][0],data[1][1],int(data[1][2]))
    ref_metadata = download_bursts(([ref_params]))
    is_ascending = ref_metadata.orbit_direction == 'ascending'
    ref_footprint = get_isce2_burst_bbox(ref_params)
    sec_footprint = get_isce2_burst_bbox(sec_params)
    insar_roi = get_region_of_interest(ref_footprint, sec_footprint, is_ascending=is_ascending)
    ###################################################################################################
    roi=" -b '" + str(insar_roi[1]) + " " + str(insar_roi[3]) + " " + str(insar_roi[0]) + " " + str(insar_roi[2]) + "'"
    print("burst的选择区域："+roi)
    os.chdir('../')
    
    data = numpy.loadtxt('./data/stack.txt',dtype=str,delimiter=',')
    IW_num = str(data[0][1][2])
    #os.system("python ../code/topsStack/stackSentinel.py -s ./data -d ./data/dem/full_res.dem.wgs84  -a ./data/aux_cal/ -o ./data/orbits -C geometry  -W slc -z 1 -r 1 -n "+IW_num+" "+roi)
    os.system("python ../code/topsStack/stackSentinel.py -s ./data -d ./data/dem/full_res.dem.wgs84  -a ./data/aux_cal/ -o ./data/orbits -C geometry -c 3 -z 1 -r 4 -n "+IW_num+" "+roi)
    line = open('./run_files/run_01_unpack_topo_reference').readline()
    line = line.replace("\n", "")
    line=str(line)+'0'
    print(line)
    open('./run_files/run_00_unpack_topo_reference','w').writelines(str(line))
    lines = open('./configs/config_reference').readlines()
    open('./configs/config_reference0','w').writelines(lines[0:13])
    open('./configs/config_reference','w').writelines(lines[0:4] + lines[16:])
    open('./configs/run_00_unpack_topo_reference','w').writelines(str(line))
    ######删除主影像合并到merged的重复命令
    line_merged = open('./run_files/run_07_merge_reference_secondary_slc').readlines()
    open('./run_files/run_07_merge_reference_secondary_slc','w').writelines(line_merged[1:])


run=[]
run.append("run_00_unpack_topo_reference")    #0
run.append('run_01_unpack_topo_reference')    #1
run.append('run_02_unpack_secondary_slc')    #2
run.append('run_03_average_baseline')    #3
run.append('run_04_fullBurst_geo2rdr')    #4
run.append('run_05_fullBurst_resample')    #5
run.append('run_06_extract_stack_valid_region')    #6
run.append('run_07_merge_reference_secondary_slc')    #7
run.append('run_08_generate_burst_igram')    #8
run.append('run_09_merge_burst_igram')    #9
run.append('run_10_filter_coherence')    #10
run.append('run_11_unwrap')    #11

def swap_burst_vrts():
    """
    Swap the VRTs generated by topsApp for the reference and secondary bursts
    To convince topsApp to process a burst pair, we need to swap the VRTs it generates for the
    reference and secondary bursts with custom VRTs that point to the actual burst rasters.
    """
    ref_vrt_list = glob.glob('./reference/**/*.vrt')
    sec_vrt_list = glob.glob('./secondarys/**/**/*.vrt')
    path=os.getcwd()
    print(str(path))
    print(ref_vrt_list)
    print(sec_vrt_list)
    if len(ref_vrt_list)>0:
        for vrt_path in (ref_vrt_list):
            print(vrt_path)
            print(str(path)+str(vrt_path)[1:])
            vrt = gdal.Open(str(path)+str(vrt_path)[1:])
            base = gdal.Open(vrt.GetFileList()[1])
            del vrt
            gdal.Translate(str(path)+str(vrt_path)[1:], base, format='VRT')
            print("修改VRT")
            del base
    if len(sec_vrt_list)>0:     
        for i in range(len(sec_vrt_list)):
            print(str(path)+str(sec_vrt_list[i])[1:])
            vrt = gdal.Open(str(path)+str(sec_vrt_list[i])[1:])
            base = gdal.Open(vrt.GetFileList()[1])
            print(base)
            del vrt
            gdal.Translate(str(path)+str(sec_vrt_list[i])[1:], base, format='VRT')
            print("修改VRT")
            del base
    
    
#要运行的步骤
     
for i in range(len(run)):
    runstep=run[i]
    if runstep=='run_01_unpack_topo_reference':
        os.system('rm -rf ./geom_reference')
    with open('./run_files/'+str(runstep),"r") as f:
        a=f.readlines()
        print(str(a[0])[0:6])
    if str(a[0])[0:6]!="python":
        for i in range(len(a)):
            a[i]="python ../code/topsStack/"+a[i]
        with open('./run_files/'+str(runstep),"w") as f:
            for i in range(len(a)):
                f.write(str(a[i]))
                
                
    os.system('python ../code/topsStack/run.py -i ./run_files/'+str(runstep)+' -p 8')
    if runstep == run[0] or runstep == run[2]:
        swap_burst_vrts()
    ##############################################################################
    #更改干涉处理中主影像数据路径
    if runstep == run[0]:
        import glob
        pairlist=glob.glob('./configs/config_generate_igram*')
        #print(pairlist)
        for i in range(len(pairlist)):
            lines = open(pairlist[i]).readlines()
            if lines[5][-10:-1]=="reference":
                lines[5]=lines[5][:-10]+'coreg_secondarys/'+pairlist[i][-17:-9]+lines[5][-1]
                open(pairlist[i],'w').writelines(lines)
 
    ##############################################################################
    #↓裁剪处理↓，若处理整个burst，后续直至结尾代码需注释
    ##############################################################################
        
    if runstep == run[1]:
        data = numpy.loadtxt('./data/stack.txt',dtype=str,delimiter=',')
        IW = str(data[0][1])
        #研究区域裁剪设置#################################################################################################
        #↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
        os.system("python ../code/crop_rdr.py -b '37.643509 37.664975 120.274337 120.311787' -lat ./geom_reference/" +IW+ "/lat_01.rdr -lon ./geom_reference/" +IW+ "/lon_01.rdr > ./geom_reference/crop.txt")
        lines = open('./geom_reference/crop.txt').readlines()
        line=str(lines[0][:-1])
        ################################################################################################################
        path= "./geom_reference/III"
        if os.path.exists(path):
            os.system('rm -rf ./geom_reference/III')
        os.system('mv ./geom_reference/'+IW+' ./geom_reference/III')
        os.system('mkdir ./geom_reference/'+IW)   

        os.system(line+' ./geom_reference/III/hgt_01.rdr ./geom_reference/'+IW+'/hgt_01.rdr')
        os.system(line+' ./geom_reference/III/incLocal_01.rdr ./geom_reference/'+IW+'/incLocal_01.rdr')
        os.system(line+' ./geom_reference/III/lat_01.rdr ./geom_reference/'+IW+'/lat_01.rdr')
        os.system(line+' ./geom_reference/III/lon_01.rdr ./geom_reference/'+IW+'/lon_01.rdr')
        os.system(line+' ./geom_reference/III/los_01.rdr ./geom_reference/'+IW+'/los_01.rdr')
        os.system(line+' ./geom_reference/III/shadowMask_01.rdr ./geom_reference/'+IW+'/shadowMask_01.rdr')
        os.system('python ../code/gdal2isce_xml.py -i ./geom_reference/'+IW+'/hgt_01.rdr')
        os.system('python ../code/gdal2isce_xml.py -i ./geom_reference/'+IW+'/incLocal_01.rdr')
        os.system('python ../code/gdal2isce_xml.py -i ./geom_reference/'+IW+'/lat_01.rdr')
        os.system('python ../code/gdal2isce_xml.py -i ./geom_reference/'+IW+'/lon_01.rdr')
        os.system('python ../code/gdal2isce_xml.py -i ./geom_reference/'+IW+'/los_01.rdr')
        os.system('python ../code/gdal2isce_xml.py -i ./geom_reference/'+IW+'/shadowMask_01.rdr')
    ################修改xml
    if runstep == run[3]:
        data = numpy.loadtxt('./data/stack.txt',dtype=str,delimiter=',')
        IW = str(data[0][1])
        size = numpy.loadtxt('./geom_reference/crop.txt',dtype=str,delimiter=' ')
        width=size[4]
        height=size[5]
        lines = open('reference/'+IW+'.xml').readlines()
        for i in range(len(lines)):
            if lines[i][:-1]=='                <property name="firstvalidline">':
                print(lines[i+1])
                lines[i+1]='                    <value>0</value>\n'
            if lines[i][:-1]=='                <property name="firstvalidsample">':
                print(lines[i+1])
                lines[i+1]='                    <value>0</value>\n'
            if lines[i][:-1]=='                <property name="numberoflines">':
                print(lines[i+1])
                lines[i+1]='                    <value>'+height+'</value>\n'
            if lines[i][:-1]=='                <property name="numberofsamples">':
                print(lines[i+1])
                lines[i+1]='                    <value>'+width+'</value>\n'
            if lines[i][:-1]=='                <property name="numberofvalidsamples">':
                print(lines[i+1])
                lines[i+1]='                    <value>'+width+'</value>\n'
        open('reference/'+IW+'.xml','w').writelines(lines)
