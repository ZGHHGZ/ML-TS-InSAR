###### 用于 SARget 下载数据后，进行 isce 处理的预处理脚本（全 burst 模式，无裁剪）
###### 在 SARget 下载路径下运行，不执行任何地理范围裁剪

# *****************************************************************************#
# *****************************************************************************#
########################      isce process      ###############################
# *****************************************************************************#
# *****************************************************************************#
import subprocess
from lxml import etree
import os
import numpy
from pathlib import Path
from osgeo import gdal
import glob
import time


def print_author_info():
    print("==============================================")
    print("Author: Guanghui Zhang")
    print("Institution: Shandong University of Science and Technology (SDUST)")
    print("Email: 2336164866@qq.com")
    print("Create Time: 2024")
    print("Description: Efficient Sentinel-1 Coregistration, Interferometry and Phase Unwrapping Pipeline (Full Burst, No Crop)")
    print("==============================================")
    time.sleep(3)


print_author_info()

##### 读取路径下的 pre_parameter.txt 文件，获取各参数 #####
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

#### tif dem 转 isce 格式 ####
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


dem_name = "./dem/dem.tif"
dem_path = Path('./dem/full_res.dem.wgs84')

#### dem_path 绝对路径 ####
dem_path = dem_path.resolve()
gdal.Translate(str(dem_path), dem_name, format='ISCE')
xml_path = tag_dem_xml_as_ellipsoidal(dem_path)
fix_image_xml(xml_path)


path = "./run_files"
if os.path.exists(path):
    print("run_files were existed")
else:
    # 全 burst 处理：使用 geometry 配准，不做地理裁剪
    os.system("python ../code/topsStack/stackSentinel.py -s ./slc -d ./dem/full_res.dem.wgs84  -a ./aux_cal/ -o ./orbits -C geometry -c "+num_connect+" -z "+num_azimuth+" -r "+num_range+" -f "+filter_strength)
    line = open('./run_files/run_01_unpack_topo_reference').readline()
    line = line.replace("\n", "")
    line = str(line) + '0'
    print(line)
    open('./run_files/run_00_unpack_topo_reference', 'w').writelines(str(line))
    lines = open('./configs/config_reference').readlines()
    open('./configs/config_reference0', 'w').writelines(lines[0:12])
    open('./configs/config_reference', 'w').writelines(lines[0:4] + lines[15:])
    ###### 删除主影像合并到 merged 的重复命令
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


# 要运行的步骤
for i in range(len(run)):
    ##########
    print("正在进行步骤：" + str(run[i]))
    runstep = run[i]
    if runstep == 'run_01_unpack_topo_reference':
        os.system('rm -rf ./geom_reference')
    with open('./run_files/' + str(runstep), "r") as f:
        a = f.readlines()
        print(str(a[0])[0:6])
    if str(a[0])[0:6] != "python":
        for j in range(len(a)):
            a[j] = "python ../code/topsStack/" + a[j]
        with open('./run_files/' + str(runstep), "w") as f:
            for j in range(len(a)):
                f.write(str(a[j]))

    if unw_mp == "":
        unw_mp = mp
    if i == 8:
        flow_mp = unw_mp
    else:
        flow_mp = mp

    if mode == "SLC":
        if i <= 7:
            os.system('python ../code/topsStack/run.py -i ./run_files/' + str(runstep) + ' -p '+str(flow_mp))
    else:
        os.system('python ../code/topsStack/run.py -i ./run_files/' + str(runstep) + ' -p '+str(flow_mp))

    ##############################################################################
    # 更改干涉处理中主影像数据路径
    if runstep == run[0]:
        pairlist = glob.glob('./configs/config_generate_igram*')
        for j in range(len(pairlist)):
            lines = open(pairlist[j]).readlines()
            if lines[5][-10:-1] == "reference":
                lines[5] = lines[5][:-10] + 'coreg_secondarys/' + pairlist[j][-17:-9] + lines[5][-1]
                open(pairlist[j], 'w').writelines(lines)

print("全部步骤执行完毕（全 burst 无裁剪模式）。")
