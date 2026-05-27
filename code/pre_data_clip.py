import glob
import re
import shutil
import multiprocessing
import subprocess
import threading
import time
from pathlib import Path

import numpy
import rasterio
import os
import numpy as np
import requests
from rasterio.windows import Window
from burst2safe.local2safe import local2safe
import dem_stitcher
import asf_search as asf
from osgeo import gdal
from lxml import etree
import warnings
import os
import shutil
from s1_orbits import fetch_for_scene
warnings.filterwarnings("ignore")
from hyp3_isce2.s1_auxcal import download_aux_cal
from datetime import datetime, timedelta
import xml.etree.ElementTree as ET

def read_par():
    #####读取路径下的pre_parameter.txt文件，获取各参数、
    lat1 = ''
    lat2 = ''
    lon1 = ''
    lon2 = ''
    mp = ''
    unw_mp = ''
    pre_parameter_path = './pre_parameter.txt'
    with open(pre_parameter_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('运行模式：'):
                mode = line.split('：')[1].strip()
            elif line.startswith('起始经度：'):
                lon1 = line.split('：')[1].strip()
            elif line.startswith('结束经度：'):
                lon2 = line.split('：')[1].strip()
            elif line.startswith('起始纬度：'):
                lat1 = line.split('：')[1].strip()
            elif line.startswith('结束纬度：'):
                lat2 = line.split('：')[1].strip()
            elif line.startswith('方位向多视数：'):
                num_azimuth = line.split('：')[1].strip()
            elif line.startswith('距离向多视数：'):
                num_range = line.split('：')[1].strip()
            elif line.startswith('相位滤波强度：'):
                filter_strength = line.split('：')[1].strip()
            elif line.startswith('影像连接数量：'):
                num_connect = line.split('：')[1].strip()
            elif line.startswith('全流程并行运行数：'):
                mp = line.split('：')[1].strip()
            elif line.startswith('解缠步骤并行运行数：'):
                unw_mp = line.split('：')[1].strip()
    return mode, lat1, lat2, lon1, lon2, num_azimuth, num_range, filter_strength, num_connect, mp, unw_mp

def read_list(filename):
    with open(filename) as f:
        data = f.readlines()
    data_url = data[1]
    name = data_url.split(",")[0]
    url_zip = data_url.split(",")[-1]
    url = url_zip[:url_zip.rfind('.')]
    return name, url

def down_aux():
    print("正在下载aux数据...")
    aux_cal_dir = Path('aux_cal')
    #aux数据下载
    download_aux_cal(aux_cal_dir)
    print("aux数据下载完成！")



def orbit_download(url):
    print("正在下载轨道数据...")
    safe_name = url.split("/")[3]
    orbit_dir = Path('orbits')
    orbit_dir.mkdir(exist_ok=True, parents=True)
    eof_file=fetch_for_scene(safe_name, str(orbit_dir))
    print('orbit: ', eof_file)
    ####################################################
    # -------------------------------
    safe_time_str = safe_name.split("_")[5]
    safe_time = datetime.strptime(
        safe_time_str,
        "%Y%m%dT%H%M%S"
    )
    # 前后30分钟
    start_time = safe_time - timedelta(minutes=30)
    end_time = safe_time + timedelta(minutes=30)

    print("SAFE time :", safe_time)
    print("Start time:", start_time)
    print("End time  :", end_time)
    # -------------------------------
    # 解析 EOF
    # -------------------------------
    tree = ET.parse(eof_file)
    print(tree)
    root = tree.getroot()
    # 找到 List_of_OSVs
    osv_container = root.find(".//List_of_OSVs")
    # 所有 OSV
    osv_list = osv_container.findall("OSV")
    # 保存符合条件的 OSV
    selected_osv = []
    for osv in osv_list:
        utc_text = osv.find("UTC").text.strip()
        utc_text = utc_text.replace("UTC=", "")
        osv_time = datetime.strptime(
            utc_text,
            "%Y-%m-%dT%H:%M:%S.%f"
        )
        if start_time <= osv_time <= end_time:
            selected_osv.append(osv)
    # -------------------------------
    # 清空原 OSV
    # -------------------------------
    for osv in osv_list:
        osv_container.remove(osv)
    # -------------------------------
    # 添加裁剪后的 OSV
    # -------------------------------
    for osv in selected_osv:
        osv_container.append(osv)
    # 更新 count
    osv_container.set("count", str(len(selected_osv)))
    # -------------------------------
    # 保存新的 EOF
    # -------------------------------
    tree.write(
        eof_file,
        encoding="UTF-8",
        xml_declaration=True
    )
    ####################################################






#DEM数据下载
def download_dem(name):
    for i in range(9):
        try:
            print("正在下载DEM数据...")
            dem_dir = Path('dem')
            dem_dir.mkdir(exist_ok=True, parents=True)
            data = asf.product_search(name)
            cord_area = data[0].geometry['coordinates'][0]
            #####获取cord_area的最大经纬度范围
            lats = [point[1] for point in cord_area]
            lons = [point[0] for point in cord_area]
            max_lat = max(lats)
            min_lat = min(lats)
            max_lon = max(lons)
            min_lon = min(lons)
            # 向外扩充0.2度
            max_lat += 0.1
            min_lat -= 0.1
            max_lon += 0.1
            min_lon -= 0.1
            X, p = dem_stitcher.stitch_dem([min_lon, min_lat, max_lon, max_lat], dem_name="glo_30")
            import rasterio
            with rasterio.open("./dem/dem.tif", 'w', **p) as ds:
                ds.write(X, 1)
                ds.update_tags(AREA_OR_POINT='Point')
            print("DEM数据下载完成！")
            break
        except:
            print("DEM数据下载失败，正在重试...")
            time.sleep(3)
            pass



####tif dem转isce格式####
def tag_dem_xml_as_ellipsoidal(dem_path: Path) -> str:
    xml_path = str(dem_path) + '.xml'
    assert Path(xml_path).exists()
    root = etree.parse(xml_path).getroot()
    element = etree.Element('property', name='reference')
    etree.SubElement(element, 'value').text = 'WGS84'
    etree.SubElement(element, 'doc').text = 'Geodetic datum'
    root.insert(0, element)
    with open(xml_path, 'wb') as file:
        file.write(etree.tostring(root, pretty_print=True))
    return xml_path

def fix_image_xml(xml_path: str) -> None:
    cmd = ['../code/fixImageXml.py', '-i', xml_path, '--full']
    subprocess.run(cmd, check=True)


def dem2isce():
    dem_name="./dem/dem.tif"
    dem_path = Path('./dem/full_res.dem.wgs84')

    ####dem_path绝对路径
    dem_path=dem_path.resolve()
    gdal.Translate(str(dem_path), dem_name, format='ISCE')
    xml_path = tag_dem_xml_as_ellipsoidal(dem_path)
    fix_image_xml(xml_path)



def donwnload_first_burst(name,url):
    try:
        os.makedirs("./cace/" + str(name))
    except:
        pass
    username = "sarget"
    password = "SARget_SARget_123"
    session = requests.Session()
    session.verify = False
    session.auth = (username, password)
    url_tif=url+".tiff"
    for i in range(99):
        if i>2:
            time.sleep(3)
        response = session.get(url_tif, allow_redirects=False)
        print("第"+str(i+1)+"次跳转",response.url)
        while response.status_code in [301, 302, 303, 307]:
            redirect_url = response.headers.get('Location')
            response = session.get(redirect_url, allow_redirects=False)
        if response.url[0:25] != url[0:25]:
            break

    tif_url = response.url
    print("最终下载链接:", tif_url)

    url_xml=url+".xml"
    print(url_xml)

    for i in range(99):
        if i>2:
            time.sleep(3)
        response = session.get(url_xml, allow_redirects=False)
        print("第"+str(i+1)+"次跳转",response.url)
        while response.status_code in [301, 302, 303, 307]:
            redirect_url = response.headers.get('Location')
            response = session.get(redirect_url, allow_redirects=False)
        if response.url[0:25] != url_xml[0:25]:
            break
    xml_url=response.url
    print("最终下载链接:", xml_url)

    out_tif= "./cace/" + name + "/" + name + ".tiff"
    out_xml = "./cace/" + name + "/" + name + ".xml"
    print("正在下载TIF...")
    r0= requests.get(tif_url, verify=False)
    with open(out_tif, "wb") as f:
        f.write(r0.content)
    print("正在下载XML...")
    r1 = requests.get(xml_url, verify=False)
    with open(out_xml, "wb") as f:
        f.write(r1.content)



    if len(glob.glob("./cace/"+name+"/*.tiff"))>=1 and len(glob.glob("./cace/"+name+"/*.xml"))>=1:
        # 自动解析信息
        parts = url.split("/")
        slc_name = parts[3]
        swath = parts[4]
        polarization = parts[5]
        burst_index = int(parts[6].replace(".zip", ""))

        # 文件路径（解压后在当前目录）
        tiff_path = "./cace/"+name+"/"+name+".tiff"
        xml_path = "./cace/"+name+"/"+name+".xml"

        # 自动生成 slc_dict
        slc_dict = {
            slc_name: {
                swath: {
                    polarization: {
                        burst_index: {
                            "DATA": str(Path(tiff_path).absolute()),
                            "METADATA": str(Path(xml_path).absolute())
                        }
                    }
                }
            }
        }

        # 输出
        print("自动生成的 slc_dict：")
        print(slc_dict)
        safe_path=local2safe(slc_dict,work_dir="./cace/"+name)
        name_path=str(safe_path)+"/name.txt"
        with open(name_path, "w") as f:
            f.write(str(name))
        try:
            shutil.move(safe_path, "./slc")
        except:
            pass



def pre_run(lat1,lat2,lon1,lon2):
    src = r"../code/pre_run"
    dst = r"./pre_run"
    if os.path.exists(dst):
        shutil.rmtree(dst)  # 删除旧目录
    shutil.copytree(src, dst)

    os.system("rm -rf ./geom_reference")
    #获取slc文件夹下的第一个文件夹的文件名
    slc_0= next(os.walk("./slc"))[1][0]
    with open("./pre_run/config_reference0", "r", encoding="utf-8") as f:
        txt = f.read()
    txt = txt.replace("xxx", "./slc/"+slc_0)
    with open("./pre_run/config_reference0", "w", encoding="utf-8") as f:
        f.write(txt)
    os.system("python ../code/topsStack/run.py -i ./pre_run/run_00_unpack_topo_reference")
    os.system("python ../code/topsStack/run.py -i ./pre_run/run_01_unpack_topo_reference")
    ######裁剪地理文件
    if lat1 != '' and lat2 != '' and lon1 != '' and lon2 != '':
        #####获取slc文件夹下以.SAFE结尾的文件夹名称，获取IW号
        slc_list = glob.glob('./slc/*.SAFE')
        slc_list_0 = slc_list[0]
        data = str(numpy.loadtxt(str(slc_list_0) + '/name.txt', dtype=str, delimiter=','))
        IW = data[10:13]
        # 研究区域裁剪设置#################################################################################################
        # ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
        os.system("python ../code/crop_rdr.py -b '" + lat1 + " " + lat2 + " " + lon1 + " " + lon2 + "' -lat ./geom_reference/" + IW + "/lat_01.rdr -lon ./geom_reference/" + IW + "/lon_01.rdr > ./geom_reference/crop.txt")
        lines = open('./geom_reference/crop.txt').readlines()
        line = str(lines[0][:-1])
        ################################################################################################################
        path = "./geom_reference/III"
        if os.path.exists(path):
            os.system('rm -rf ./geom_reference/III')
        os.system('mv ./geom_reference/' + IW + ' ./geom_reference/III')
        os.system('mkdir ./geom_reference/' + IW)
        os.system(line + ' ./geom_reference/III/hgt_01.rdr ./geom_reference/' + IW + '/hgt_01.rdr')
        os.system(line + ' ./geom_reference/III/incLocal_01.rdr ./geom_reference/' + IW + '/incLocal_01.rdr')
        os.system(line + ' ./geom_reference/III/lat_01.rdr ./geom_reference/' + IW + '/lat_01.rdr')
        os.system(line + ' ./geom_reference/III/lon_01.rdr ./geom_reference/' + IW + '/lon_01.rdr')
        os.system(line + ' ./geom_reference/III/los_01.rdr ./geom_reference/' + IW + '/los_01.rdr')
        os.system(line + ' ./geom_reference/III/shadowMask_01.rdr ./geom_reference/' + IW + '/shadowMask_01.rdr')
        os.system('python ../code/gdal2isce_xml.py -i ./geom_reference/' + IW + '/hgt_01.rdr')
        os.system('python ../code/gdal2isce_xml.py -i ./geom_reference/' + IW + '/incLocal_01.rdr')
        os.system('python ../code/gdal2isce_xml.py -i ./geom_reference/' + IW + '/lat_01.rdr')
        os.system('python ../code/gdal2isce_xml.py -i ./geom_reference/' + IW + '/lon_01.rdr')
        os.system('python ../code/gdal2isce_xml.py -i ./geom_reference/' + IW + '/los_01.rdr')
        os.system('python ../code/gdal2isce_xml.py -i ./geom_reference/' + IW + '/shadowMask_01.rdr')
        os.system("rm -rf ./geom_reference/III")



def list_clip_burst(filename):
    with open(filename) as f:
        data = f.readlines()
    data_url = data[1:]
    name_url_list = []
    for line in data_url:
        parts = line.strip().split(",")
        name = parts[0]
        url_zip = parts[-1]
        url = url_zip[:url_zip.rfind('.')]
        name_url_list.append([name, url])

    txt_list = glob.glob("./slc/*/name.txt")
    exist_names = set()
    for txt in txt_list:
        with open(txt, "r") as f:
            name = f.readline().strip()
            if name:
                exist_names.add(name)
    name_url_list = [row for row in name_url_list if row[0] not in exist_names]
    print(len(name_url_list))
    return name_url_list



def mul_download_clip_burst(list_file, mp):

    for i in range(99):
        if i > 2:
            time.sleep(3)
        try:
            os.system("rm -rf ./cace")
            os.system("mkdir ./cace")
            os.system("mkdir ./orbits")
            os.system("mkdir ./slc")
            name_url_list = list_clip_burst(list_file)
            if len(name_url_list) == 0:
                break
            pool = multiprocessing.Pool(processes=int(mp))
            pool.map(download_clip_burst, name_url_list)
            pool.close()
            pool.join()
        except:
            pass



def download_clip_burst(name_url):
    name, url = name_url
    print(name_url)
    try:
        os.makedirs("./cace/" + str(name))
    except:
        pass
    username = "sarget"
    password = "SARget_SARget_123"
    session = requests.Session()
    session.verify = False
    session.auth = (username, password)
    url_tif=url+".tiff"
    for i in range(99):
        if i>2:
            time.sleep(1)
        response = session.get(url_tif, allow_redirects=False)
        print("第"+str(i+1)+"次跳转",response.url)
        while response.status_code in [301, 302, 303, 307]:
            redirect_url = response.headers.get('Location')
            response = session.get(redirect_url, allow_redirects=False)
        if response.url[0:25] != url[0:25]:
            break

    tif_url = response.url
    print("最终下载链接:", tif_url)

    url_xml=url+".xml"
    print(url_xml)

    for i in range(99):
        if i>2:
            time.sleep(1)
        response = session.get(url_xml, allow_redirects=False)
        print("第"+str(i+1)+"次跳转",response.url)
        while response.status_code in [301, 302, 303, 307]:
            redirect_url = response.headers.get('Location')
            response = session.get(redirect_url, allow_redirects=False)
        if response.url[0:25] != url_xml[0:25]:
            break
    xml_url=response.url
    print("最终下载链接:", xml_url)

    out_tif= "./cace/" + name + "/" + name + ".tiff"
    out_xml = "./cace/" + name + "/" + name + ".xml"
    with open("./geom_reference/crop.txt", "r") as f:
        line = f.readline().strip()
        match = re.search(r"-srcwin\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)", line)
        xoff, yoff, xsize, ysize = map(int, match.groups())
        print(xoff, yoff, xsize, ysize)
    row_buffer = 25
    col_buffer = 100
    print("正在下载TIFF: "+str(name))
    with rasterio.open(tif_url) as src:
        # 创建一个全尺寸数组（填充为 NoData）
        nodata_value = 0
        full_data = np.full((src.height, src.width), nodata_value, dtype="complex64")
        row_start_oringin = int(yoff )
        col_start_oringin = int(xoff)
        height_oringin = int(ysize)
        width_oringin = int(xsize)
        if row_start_oringin <= row_buffer:
            row_start_new=0
        else:
            row_start_new=row_start_oringin-row_buffer
        if col_start_oringin <= col_buffer:
            col_start_new=0
        else:
            col_start_new=col_start_oringin-col_buffer

        if row_start_new==0:
            height_new=height_oringin+row_buffer+row_start_oringin
        else:
            height_new=height_oringin+row_buffer*2
        height_new=min(height_new,src.height-row_start_new)
        if col_start_new==0:
            width_new=width_oringin+col_buffer+col_start_oringin
        else:
            width_new=width_oringin+col_buffer*2
        width_new=min(width_new,src.width-col_start_new)
        window = Window(col_start_new, row_start_new, width_new, height_new)
        # 读取窗口数据
        data = src.read(1, window=window)
        # 写入对应位置
        full_data[
            row_start_new:row_start_new + height_new,
            col_start_new:col_start_new + width_new
        ] = data
        profile = src.profile.copy()
        profile.update(
            compress='lzw'
        )
    # 写出 GeoTIFF
    with rasterio.open(out_tif, 'w', **profile) as dst:
            dst.write(full_data, 1)
    print("正在下载XML: "+str(name))
    r1 = requests.get(xml_url, verify=False)
    with open(out_xml, "wb") as f:
        f.write(r1.content)

    if len(glob.glob("./cace/" + name + "/*.tiff")) >= 1 and len(glob.glob("./cace/" + name + "/*.xml")) >= 1:
        # 自动解析信息
        parts = url.split("/")
        slc_name = parts[3]
        swath = parts[4]
        polarization = parts[5]
        burst_index = int(parts[6].replace(".zip", ""))

        # 文件路径（解压后在当前目录）
        tiff_path = "./cace/" + name + "/" + name + ".tiff"
        xml_path = "./cace/" + name + "/" + name + ".xml"

        # 自动生成 slc_dict
        slc_dict = {
            slc_name: {
                swath: {
                    polarization: {
                        burst_index: {
                            "DATA": str(Path(tiff_path).absolute()),
                            "METADATA": str(Path(xml_path).absolute())
                        }
                    }
                }
            }
        }

        # 输出看看
        print("自动生成的 slc_dict：")
        print(slc_dict)
        safe_path=local2safe(slc_dict, work_dir="./cace/" + name)
        input_path = glob.glob("./cace/" + name + "/*/" + "/measurement/*.tiff")
        in_file = input_path[0]
        print("压缩: "+in_file)
        name_path=str(safe_path)+"/name.txt"
        with open(name_path, "w") as f:
            f.write(str(name))
        temp_file = in_file + ".tmp"
        with rasterio.open(in_file) as src:
            profile = src.profile
            profile.update(compress='lzw')
            with rasterio.open(temp_file, 'w', **profile) as dst:
                dst.write(src.read())
        # 用临时文件替换原文件（实现覆盖）
        os.replace(temp_file, in_file)
        orbit_download(url)
        try:
            shutil.move(safe_path, "./slc")
        except:
            pass





def make_run_list():
    path = "./run_files"
    if os.path.exists(path):
        print("run_files were existed")
    else:
        os.system(
            "python ../code/topsStack/stackSentinel.py -s ./slc -d ./dem/full_res.dem.wgs84  -a ./aux_cal/ -o ./orbits -C geometry -c " + num_connect + " -z " + num_azimuth + " -r " + num_range + " -f " + filter_strength)
        line = open('./run_files/run_01_unpack_topo_reference').readline()
        line = line.replace("\n", "")
        line = str(line) + '0'
        open('./run_files/run_00_unpack_topo_reference', 'w').writelines(str(line))
        lines = open('./configs/config_reference').readlines()
        open('./configs/config_reference0', 'w').writelines(lines[0:12])
        open('./configs/config_reference', 'w').writelines(lines[0:4] + lines[15:])
        ######删除主影像合并到merged的重复命令
        line_merged = open('./run_files/run_07_merge_reference_secondary_slc').readlines()
        open('./run_files/run_07_merge_reference_secondary_slc', 'w').writelines(line_merged[1:])

    run = []
    run.append("run_00_unpack_topo_reference")  # 0
    run.append('run_01_unpack_topo_reference')  # 1
    run.append('run_02_unpack_secondary_slc')  # 2
    run.append('run_03_average_baseline')  # 3
    run.append('run_04_fullBurst_geo2rdr')  # 4
    run.append('run_05_fullBurst_resample')  # 5
    run.append('run_06_extract_stack_valid_region')  # 6
    run.append('run_07_merge_reference_secondary_slc')  # 7
    run.append('run_08_generate_burst_igram')  # 8
    run.append('run_09_merge_burst_igram')  # 9
    run.append('run_10_filter_coherence')  # 10
    run.append('run_11_unwrap')  # 11
    return run


def run_run_run(mode,mp,unw_mp):
    run=make_run_list()
    # 要运行的步骤
    ##############################################################################
    # 更改干涉处理中主影像数据路径
    pairlist = glob.glob('./configs/config_generate_igram*')
    # print(pairlist)
    for i in range(len(pairlist)):
        lines = open(pairlist[i]).readlines()
        if lines[5][-10:-1] == "reference":
            lines[5] = lines[5][:-10] + 'coreg_secondarys/' + pairlist[i][-17:-9] + lines[5][-1]
            open(pairlist[i], 'w').writelines(lines)
    ##############################################################################

    for i in range(2,len(run)):
        ##########
        print("正在进行步骤：" + str(run[i]))
        runstep = run[i]
        if runstep == 'run_01_unpack_topo_reference':
            os.system('rm -rf ./geom_reference')
        with open('./run_files/' + str(runstep), "r") as f:
            a = f.readlines()
            print(str(a[0])[0:6])
        if str(a[0])[0:6] != "python":
            for i in range(len(a)):
                a[i] = "python ../code/topsStack/" + a[i]
            with open('./run_files/' + str(runstep), "w") as f:
                for i in range(len(a)):
                    f.write(str(a[i]))

        if unw_mp == "":
            unw_mp = mp
        if i == 8:
            flow_mp = unw_mp
        else:
            flow_mp = mp

        if mode == "SLC":
            if i <= 7:
                os.system('python ../code/topsStack/run.py -i ./run_files/' + str(runstep) + ' -p ' + str(flow_mp))
        else:
            os.system('python ../code/topsStack/run.py -i ./run_files/' + str(runstep) + ' -p ' + str(flow_mp))

            ################修改xml
            if runstep == run[3]:
                #####获取slc文件夹下以.SAFE结尾的文件夹名称，获取IW号
                slc_list = glob.glob('./slc/*.SAFE')
                slc_list_0 = slc_list[0]
                data = str(numpy.loadtxt(str(slc_list_0) + '/name.txt', dtype=str, delimiter=','))
                IW = data[10:13]
                size = numpy.loadtxt('./geom_reference/crop.txt', dtype=str, delimiter=' ')
                width = size[4]
                height = size[5]
                lines = open('reference/' + IW + '.xml').readlines()
                for i in range(len(lines)):
                    if lines[i][:-1] == '                <property name="firstvalidline">':
                        print(lines[i + 1])
                        lines[i + 1] = '                    <value>0</value>\n'
                    if lines[i][:-1] == '                <property name="firstvalidsample">':
                        print(lines[i + 1])
                        lines[i + 1] = '                    <value>0</value>\n'
                    if lines[i][:-1] == '                <property name="numberoflines">':
                        print(lines[i + 1])
                        lines[i + 1] = '                    <value>' + height + '</value>\n'
                    if lines[i][:-1] == '                <property name="numberofsamples">':
                        print(lines[i + 1])
                        lines[i + 1] = '                    <value>' + width + '</value>\n'
                    if lines[i][:-1] == '                <property name="numberofvalidsamples">':
                        print(lines[i + 1])
                        lines[i + 1] = '                    <value>' + width + '</value>\n'
                open('reference/' + IW + '.xml', 'w').writelines(lines)
        ####################################################################################

if __name__ == "__main__":
    list_file = "download_list.sarget"
    name, url = read_list(list_file)
    mode, lat1, lat2, lon1, lon2, num_azimuth, num_range, filter_strength, num_connect, mp, unw_mp=read_par()
    os.system("rm -rf ./cace")
    os.system("mkdir ./cace")
    try:
        os.system("mkdir ./slc")
    except:
        pass
    try:
        os.system("mkdir ./dem")
    except:
        pass
    try:
        os.system("mkdir ./orbits")
    except:
        pass

    crop_file_path=Path("./geom_reference/crop.txt")
    if crop_file_path.exists():
        print("主影影像地理编码已完成，自动跳过次步骤。。。")
    else:
        donwnload_first_burst(name,url)
        down_aux()
        orbit_download(url)
        download_dem(name)
        dem2isce()
        pre_run(lat1,lat2,lon1,lon2)
        os.system("rm -rf ./slc")
        os.system("rm -rf ./orbits")

    mul_download_clip_burst(list_file,mp)
    run_run_run(mode,mp,unw_mp)