######用于SARget下载数据后，进行isce处理的预处理脚本，在SARget下载路径下运行
######多burst版本：对 geom_reference 下所有 IW* 子目录（每个 burst）均执行地理范围裁剪，
######同步裁剪 reference/IW*/burst_*.slc.vrt 到对应窗口，并更新 reference/IWn.xml 中每个 burst
######的尺寸、sensingStart/sensingStop/burstStartUTC/burstStopUTC/startingRange 信息，
######保证后续 geo2rdr、resample、merge 等步骤在裁剪后的参考网格上运行。

# *****************************************************************************#
# *****************************************************************************#
########################      isce process      ###############################
# *****************************************************************************#
# *****************************************************************************#
import subprocess
import sys
from lxml import etree
import os
import re
import shutil
import copy
import json
import math
import numpy
from pathlib import Path
from osgeo import gdal
import glob
import time
from datetime import datetime, timedelta


def print_author_info():
    print("==============================================")
    print("Author: Guanghui Zhang")
    print("Institution: Shandong University of Science and Technology (SDUST)")
    print("Email: 2336164866@qq.com")
    print("Create Time: 2026")
    print("Description: Efficient Sentinel-1 Coregistration, Interferometry and Phase Unwrapping Pipeline (Multi-Burst Crop)")
    print("==============================================")
    time.sleep(3)


print_author_info()
#####读取路径下的pre_parameter.txt文件，获取各参数、
lat1=''
lat2=''
lon1=''
lon2=''
mp=''
unw_mp=''
pre_parameter_path = './pre_parameter.txt'
with open(pre_parameter_path, 'r') as f:
    lines = f.readlines()
    for line in lines:
        if line.startswith('运行模式：'):
            mode = line.split('：')[1].strip()
        elif line.startswith('起始经度：'):
            lon1= line.split('：')[1].strip()
        elif line.startswith('结束经度：'):
            lon2= line.split('：')[1].strip()
        elif line.startswith('起始纬度：'):
            lat1= line.split('：')[1].strip()
        elif line.startswith('结束纬度：'):
            lat2= line.split('：')[1].strip()
        elif line.startswith('方位向多视数：'):
            num_azimuth= line.split('：')[1].strip()
        elif  line.startswith('距离向多视数：'):
            num_range= line.split('：')[1].strip()
        elif line.startswith('相位滤波强度：'):
            filter_strength= line.split('：')[1].strip()
        elif line.startswith('影像连接数量：'):
            num_connect= line.split('：')[1].strip()
        elif line.startswith('全流程并行运行数：'):
            mp= line.split('：')[1].strip()
        elif line.startswith('解缠步骤并行运行数：'):
            unw_mp= line.split('：')[1].strip()

# 只有四个边界值都提供时才启用经纬度二次裁剪；
# 任意一项缺失时保持原有的“非零在线数据窗”处理。
HAS_GEO_BBOX = all(value != '' for value in (lat1, lat2, lon1, lon2))

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

def get_iw_list() -> list:
    """返回 geom_reference 下所有 IW* 子目录名（多burst），如 ['IW1','IW2','IW3']。"""
    iw_dirs = sorted(glob.glob('./geom_reference/IW*'))
    iw_list = []
    for p in iw_dirs:
        if os.path.isdir(p):
            name = os.path.basename(p.rstrip('/'))
            if re.match(r'^IW\d+$', name):
                iw_list.append(name)
    return iw_list


def get_bursts_in_iw(iw_dir: str) -> list:
    """返回某 IW 目录下所有 burst 序号（零填充两位），依据 lat_*.rdr 推断。
    例：存在 lat_01.rdr / lat_02.rdr 则返回 ['01','02']。"""
    lat_files = sorted(glob.glob(os.path.join(iw_dir, 'lat_*.rdr')))
    bursts = []
    for lf in lat_files:
        m = re.search(r'lat_(\d+)\.rdr$', os.path.basename(lf))
        if m:
            nn = m.group(1).zfill(2)
            if nn not in bursts:
                bursts.append(nn)
    return bursts


def find_nonzero_slc_window(vrt_path: str, block_rows: int = 64) -> tuple:
    """Find the online-downloaded nonzero rectangle in a full burst VRT."""
    dataset = gdal.Open(vrt_path, gdal.GA_ReadOnly)
    if dataset is None:
        raise RuntimeError('无法打开 reference burst VRT: ' + vrt_path)
    width, height = dataset.RasterXSize, dataset.RasterYSize
    min_x, min_y, max_x, max_y = width, height, -1, -1
    band = dataset.GetRasterBand(1)
    for y0 in range(0, height, block_rows):
        count = min(block_rows, height - y0)
        data = band.ReadAsArray(0, y0, width, count)
        if data is None:
            dataset = None
            raise RuntimeError('读取 reference burst VRT 失败: ' + vrt_path)
        rows, cols = numpy.nonzero(numpy.abs(data) > 0)
        if rows.size:
            min_x = min(min_x, int(cols.min()))
            max_x = max(max_x, int(cols.max()))
            min_y = min(min_y, y0 + int(rows.min()))
            max_y = max(max_y, y0 + int(rows.max()))
    dataset = None
    if max_x < min_x or max_y < min_y:
        raise RuntimeError('reference burst 没有非零在线数据: ' + vrt_path)
    return min_x, min_y, max_x - min_x + 1, max_y - min_y + 1, width, height


def build_product_nonzero_crop_info(product_root: str, IW: str,
                                    allowed_bursts=None) -> dict:
    """Build this acquisition's own online-data windows from burst VRTs."""
    allowed = ({str(value).zfill(2) for value in allowed_bursts}
               if allowed_bursts is not None else None)
    iw_dir = os.path.join(product_root, IW)
    windows = {}
    for vrt_path in sorted(glob.glob(os.path.join(iw_dir, 'burst_*.slc.vrt'))):
        match = re.search(r'burst_(\d+)\.slc\.vrt$', os.path.basename(vrt_path))
        if not match:
            continue
        nn = match.group(1).zfill(2)
        if allowed is not None and nn not in allowed:
            continue
        xoff, yoff, width, height, full_width, full_height = (
            find_nonzero_slc_window(vrt_path))
        windows[nn] = {
            'xoff': xoff, 'yoff': yoff,
            'width': width, 'height': height,
            'geometry_width': full_width,
            'geometry_length': full_height,
            'footprint': [],
        }
    return windows


def prepare_reference_crop_before_topo() -> None:
    """Crop reference radar grids before topo so topo only covers downloaded data."""
    global REFERENCE_PRETOPO_CROP
    iw_list = get_product_iw_list('./reference')
    if not iw_list:
        raise RuntimeError('run_00 后 reference 中未发现 IW*.xml')
    for IW in iw_list:
        iw_dir = os.path.join('./reference', IW)
        windows = {}
        for vrt_path in sorted(glob.glob(os.path.join(iw_dir, 'burst_*.slc.vrt'))):
            match = re.search(r'burst_(\d+)\.slc\.vrt$', os.path.basename(vrt_path))
            if not match:
                continue
            nn = match.group(1).zfill(2)
            xoff, yoff, width, height, full_width, full_height = (
                find_nonzero_slc_window(vrt_path))
            windows[nn] = {
                'xoff': xoff, 'yoff': yoff,
                'width': width, 'height': height,
                'geometry_width': full_width,
                'geometry_length': full_height,
                'footprint': [],
            }
            print('>>> topo 前检测 ' + IW + ' burst ' + nn +
                  ' 在线窗口: xoff=' + str(xoff) + ' yoff=' + str(yoff) +
                  ' width=' + str(width) + ' height=' + str(height) +
                  ' / full=' + str(full_width) + 'x' + str(full_height))
        if not windows:
            prune_product_bursts('./reference', IW, {})
            continue
        ref_xml = os.path.join('./reference', IW + '.xml')
        add_original_burst_metadata(ref_xml, windows)
        REF_CROP_INFO[IW] = windows
    if not REF_CROP_INFO:
        raise RuntimeError('所有 reference IW 都没有在线下载的非零数据')
    reference_date = get_product_acquisition_date('./reference')
    ACQUISITION_CROP_INFO[reference_date] = copy.deepcopy(REF_CROP_INFO)
    for IW, crop_info in REF_CROP_INFO.items():
        print('>>> topo 前裁剪 reference ' + IW + ' SLC VRT/XML...')
        crop_reference_burst_slc(IW, crop_info)
        prune_product_bursts('./reference', IW, crop_info)
        update_iw_xml_geometry('./reference/' + IW + '.xml', crop_info)
    REFERENCE_PRETOPO_CROP = True
    print('>>> reference 已在 topo 前裁剪；run_01 topo 将只计算在线数据区域')
    # 写标记文件，供 run_01 判定磁盘 reference 是否已完成无黑边裁剪
    # （避免 crop_state.json 残留旧状态导致裁剪被跳过、topo 误处理全幅）。
    try:
        open('./reference/.pretopo_cropped', 'w').close()
    except Exception as e:
        print('警告：写入 pre-topo 裁剪标记失败: ' + str(e))


def calculate_iw_crop_info(IW: str, lat1, lat2, lon1, lon2) -> dict:
    """计算单个 IW 的 burst 裁剪窗口，但暂不裁剪文件。

    必须先收集全部 IW 的窗口，才能在绝对方位时间上对齐不同 IW。
    """
    iw_dir = './geom_reference/' + IW
    bursts = get_bursts_in_iw(iw_dir)
    if not bursts:
        print("警告：" + iw_dir + " 下未检测到 lat_*.rdr，跳过该 IW 裁剪")
        return {}
    print(">>> 正在计算 IW [" + IW + "] 的裁剪窗口，检测到 burst: " + ", ".join(bursts))

    # 先计算每个 burst 覆盖研究区所需的窗口。方位向窗口必须分别保留，
    # 因为相邻 burst 的绝对成像时间不同；距离向则在下一步统一为所有
    # burst 窗口的并集，避免 merge 后出现左右阶梯状 NoData 缺口。
    raw_crop_info = {}
    for nn in bursts:
        # 1) 用该 burst 自身的 lat/lon 计算其在自身雷达坐标系下的裁剪窗口
        crop_txt = './geom_reference/crop_' + IW + '_' + nn + '.txt'
        crop_kml = './geom_reference/crop_extent_' + IW + '_' + nn + '.kml'
        command = [
            sys.executable, '../code/crop_rdr.py',
            '-b', lat1 + ' ' + lat2 + ' ' + lon1 + ' ' + lon2,
            '-lat', iw_dir + '/lat_' + nn + '.rdr',
            '-lon', iw_dir + '/lon_' + nn + '.rdr',
            '-kml', crop_kml,
        ]
        result = subprocess.run(command, text=True, capture_output=True)
        if result.returncode != 0:
            no_coverage = ('不在 lat/lon 数据覆盖范围内' in result.stderr or
                           'BBox 不在覆盖范围内' in result.stderr)
            if no_coverage:
                print("    跳过 " + IW + " burst " + nn +
                      "：裁剪范围与该 burst 无交集")
                continue
            raise RuntimeError(
                "计算 " + IW + " burst " + nn + " 裁剪窗口失败:\n" +
                (result.stderr.strip() or result.stdout.strip()))
        if result.stderr:
            # 保留 crop_rdr.py 的 KML 写出信息和 GDAL 警告。
            print(result.stderr, end='' if result.stderr.endswith('\n') else '\n')
        with open(crop_txt, 'w', encoding='utf-8') as f:
            f.write(result.stdout)
        # 2) 记录该 burst 的原始裁剪窗口。
        #    crop_rdr.py 输出: gdal_translate -srcwin <xoff> <yoff> <width> <height> ...
        size = numpy.loadtxt(crop_txt, dtype=str, delimiter=' ')
        size = numpy.atleast_1d(size)
        if size.size < 6:
            raise RuntimeError(IW + " burst " + nn +
                               " 的裁剪窗口输出格式无效: " + result.stdout.strip())
        lat_dataset = gdal.Open(iw_dir + "/lat_" + nn + ".rdr", gdal.GA_ReadOnly)
        if lat_dataset is None:
            raise RuntimeError("无法读取 " + iw_dir + "/lat_" + nn + ".rdr")
        source_width = int(lat_dataset.RasterXSize)
        source_height = int(lat_dataset.RasterYSize)
        lat_dataset = None

        raw_crop_info[nn] = {
            'xoff': str(size[2]),
            'yoff': str(size[3]),
            'width': str(size[4]),
            'height': str(size[5]),
            'source_width': source_width,
            'source_height': source_height,
        }
        print("    burst " + nn + " raw crop: xoff=" + raw_crop_info[nn]['xoff'] +
              " yoff=" + raw_crop_info[nn]['yoff'] + " width=" + raw_crop_info[nn]['width'] +
              " height=" + raw_crop_info[nn]['height'])

    if not raw_crop_info:
        print(">>> 跳过 IW [" + IW + "]：其中没有 burst 覆盖裁剪范围")
        return {}

    # 所有 burst 使用相同的距离向起点和宽度。使用并集而不是交集，保证
    # 每个 burst 根据经纬度计算出的目标范围都不会被再次截掉。
    common_xoff = min(int(info['xoff']) for info in raw_crop_info.values())
    common_xend = max(int(info['xoff']) + int(info['width'])
                      for info in raw_crop_info.values())
    common_width = common_xend - common_xoff
    if common_width <= 0:
        raise RuntimeError(IW + " 的公共距离向裁剪窗口无效")
    for nn, info in raw_crop_info.items():
        if common_xoff < 0 or common_xend > int(info['source_width']):
            raise RuntimeError(
                IW + " burst " + nn + " 无法容纳公共距离向窗口: xoff=" +
                str(common_xoff) + " width=" + str(common_width) +
                " source_width=" + str(info['source_width']))
    print("    公共距离向窗口: xoff=" + str(common_xoff) +
          " width=" + str(common_width))

    crop_info = {}
    for nn in sorted(raw_crop_info):
        raw = raw_crop_info[nn]
        crop_info[nn] = {
            'xoff': str(common_xoff),
            'yoff': raw['yoff'],
            'width': str(common_width),
            'height': raw['height'],
        }

    return crop_info


def apply_iw_geom_crop(IW: str, crop_info: dict) -> None:
    """按已经完成 burst/IW 对齐的窗口裁剪单个 IW 的全部 geometry。"""
    iw_dir = './geom_reference/' + IW
    tmp = './geom_reference/_III_' + IW
    geom_bases = ['hgt', 'incLocal', 'lat', 'lon', 'los', 'shadowMask']

    if os.path.exists(tmp):
        shutil.rmtree(tmp)
    shutil.move(iw_dir, tmp)
    os.makedirs(iw_dir)

    completed = False
    try:
        for nn in sorted(crop_info):
            info = crop_info[nn]
            crop_line = ('gdal_translate -srcwin ' + info['xoff'] + ' ' +
                         info['yoff'] + ' ' + info['width'] + ' ' +
                         info['height'] + ' -of envi -co INTERLEAVE=BIP')
            for b in geom_bases:
                src = tmp + "/" + b + "_" + nn + ".rdr"
                dst = iw_dir + "/" + b + "_" + nn + ".rdr"
                if os.system(crop_line + ' ' + src + ' ' + dst) != 0:
                    raise RuntimeError("裁剪 geometry 失败: " + src)
                if os.system('python ../code/gdal2isce_xml.py -i ' + dst) != 0:
                    raise RuntimeError("生成 ISCE XML 失败: " + dst)
            print("    burst " + nn + " aligned crop: xoff=" + info['xoff'] +
                  " yoff=" + info['yoff'] + " width=" + info['width'] +
                  " height=" + info['height'])
        completed = True
    finally:
        if completed:
            shutil.rmtree(tmp)


def _set_prop(parent, name, value):
    """在 parent 下所有 name=name 的 property 中写入 value。"""
    for prop in parent.findall("property[@name='%s']" % name):
        val = prop.find('value')
        if val is not None:
            val.text = str(value)


def _get_prop(parent, name, default=None):
    """读取 parent 下第一个 name=name 的 property 的 value 文本。"""
    for prop in parent.findall("property[@name='%s']" % name):
        val = prop.find('value')
        if val is not None:
            return val.text
    return default


def get_product_iw_list(product_root: str) -> list:
    """返回 reference/secondary 产品目录中实际存在的 IW 目录名。"""
    iw_list = []
    for path in sorted(glob.glob(os.path.join(product_root, 'IW*'))):
        name = os.path.basename(path.rstrip('/'))
        if os.path.isdir(path) and re.match(r'^IW\d+$', name):
            iw_list.append(name)
    return iw_list


def prune_product_bursts(product_root: str, IW: str, crop_info: dict) -> None:
    """从生成的 reference/secondary 产品中移除未覆盖的 burst。

    只操作当前流程生成的中间产品。保留 burst 的原始编号，保证 reference 与
    secondary 的 burstNumber、文件名及相位参考信息保持一致。
    """
    keep = set(crop_info)
    iw_dir = os.path.join(product_root, IW)
    xml_path = os.path.join(product_root, IW + '.xml')

    if not keep:
        if os.path.isdir(iw_dir):
            shutil.rmtree(iw_dir)
        if os.path.exists(xml_path):
            os.remove(xml_path)
        print("    已停用 " + product_root + "/" + IW +
              "（该 IW 没有 burst 覆盖裁剪范围）")
        return

    # 删除未保留 burst 的 SLC/VRT/XML/HDR 等同名前缀文件，防止后续 merge
    # 通过文件 glob 再次把它们加入处理列表。
    removed_numbers = set()
    if os.path.isdir(iw_dir):
        for path in glob.glob(os.path.join(iw_dir, 'burst_*')):
            match = re.match(r'^burst_(\d+)', os.path.basename(path))
            if match and match.group(1).zfill(2) not in keep:
                removed_numbers.add(match.group(1).zfill(2))
                if os.path.isdir(path):
                    shutil.rmtree(path)
                else:
                    os.remove(path)

    if not os.path.exists(xml_path):
        if removed_numbers:
            print("    已移除 " + product_root + "/" + IW + " 未覆盖 burst 文件: " +
                  ", ".join(sorted(removed_numbers)))
        return

    tree = etree.parse(xml_path)
    root = tree.getroot()
    bursts_comp = root.find(".//component[@name='bursts']")
    if bursts_comp is None:
        raise RuntimeError(xml_path + " 中未找到 bursts 组件，无法跳过未覆盖 burst")

    kept_names = []
    removed_names = []
    for comp in list(bursts_comp):
        name = comp.get('name') or ''
        if not (isinstance(comp.tag, str) and comp.tag == 'component'
                and name.startswith('burst')):
            continue
        nn = name[5:].zfill(2)
        if nn in keep:
            kept_names.append(name)
        else:
            removed_names.append(name)
            bursts_comp.remove(comp)

    missing = sorted(keep - {name[5:].zfill(2) for name in kept_names})
    if missing:
        raise RuntimeError(xml_path + " 缺少需保留的 burst: " + ", ".join(missing))

    for prop in bursts_comp.findall("property[@name='name']"):
        value = prop.find('value')
        if value is not None:
            value.text = str(kept_names)
    _set_prop(root, 'numberofbursts', str(len(kept_names)))
    tree.write(xml_path, pretty_print=True)

    removed = sorted(set(removed_names) | {'burst' + str(int(x)) for x in removed_numbers})
    if removed:
        print("    已从 " + xml_path + " 跳过未覆盖 burst: " + ", ".join(removed))


def add_original_burst_metadata(xml_path: str, crop_info: dict) -> None:
    """把裁剪前的 TOPS 相位参考网格保存到 crop_info。

    ISCE2 的方位载频以完整 burst 的中心行和原始 startingRange 为参考。
    numberOfLines/startingRange 在裁剪后会改变，因此必须在改写 IW*.xml 前保存
    原始值，供 resamp_withCarrier.py 在裁剪网格上重建原始载频。
    """
    if not os.path.exists(xml_path):
        raise FileNotFoundError("未找到裁剪前的 burst 元数据: " + xml_path)

    tree = etree.parse(xml_path)
    bursts_comp = tree.getroot().find(".//component[@name='bursts']")
    if bursts_comp is None:
        raise RuntimeError(xml_path + " 中未找到 bursts 组件")

    found = set()
    for comp in bursts_comp:
        name = comp.get('name') or ''
        if not (isinstance(comp.tag, str) and comp.tag == 'component'
                and name.startswith('burst')):
            continue
        nn = name[5:].zfill(2)
        if nn not in crop_info:
            continue
        info = crop_info[nn]
        info['original_width'] = int(_get_prop(comp, 'numberofsamples', '0') or '0')
        info['original_height'] = int(_get_prop(comp, 'numberoflines', '0') or '0')
        info['original_starting_range'] = float(_get_prop(comp, 'startingrange', '0') or '0')
        info['original_sensing_start'] = _get_prop(comp, 'sensingstart')
        info['original_sensing_stop'] = _get_prop(comp, 'sensingstop')
        info['azimuth_time_interval'] = float(
            _get_prop(comp, 'azimuthtimeinterval', '0') or '0')
        if info['original_width'] <= 0 or info['original_height'] <= 0:
            raise RuntimeError(xml_path + " 中 " + name + " 的原始尺寸无效")
        if info['azimuth_time_interval'] <= 0:
            raise RuntimeError(xml_path + " 中 " + name + " 的方位时间间隔无效")
        found.add(nn)

    missing = sorted(set(crop_info) - found)
    if missing:
        raise RuntimeError(xml_path + " 缺少 burst 元数据: " + ", ".join(missing))


def align_iw_azimuth_windows(crop_info_by_iw: dict) -> None:
    """把不同 IW 的总体裁剪范围对齐到相同的绝对方位时间区间。

    burst 的 yoff 不能跨 IW 直接比较，因为各 IW/burst 的 sensingStart 不同。
    这里先求每个原始窗口的绝对时间，再取所有 IW 的时间并集，并将该并集
    映射回每个 burst 的本地行号。相邻 burst 只承担与自身原始时段相交的部分。
    """
    fmt = '%Y-%m-%d %H:%M:%S.%f'
    raw_intervals = []
    for IW, bursts in crop_info_by_iw.items():
        for nn, info in bursts.items():
            original_start = datetime.strptime(info['original_sensing_start'], fmt)
            dt = float(info['azimuth_time_interval'])
            crop_start = original_start + timedelta(seconds=int(info['yoff']) * dt)
            crop_stop = crop_start + timedelta(seconds=(int(info['height']) - 1) * dt)
            raw_intervals.append((crop_start, crop_stop, IW, nn))

    if not raw_intervals:
        return
    common_start = min(item[0] for item in raw_intervals)
    common_stop = max(item[1] for item in raw_intervals)
    print(">>> 各 IW 公共绝对方位时间: " + common_start.strftime(fmt) +
          " ~ " + common_stop.strftime(fmt))

    for IW, bursts in crop_info_by_iw.items():
        new_intervals = []
        for nn, info in bursts.items():
            original_start = datetime.strptime(info['original_sensing_start'], fmt)
            dt = float(info['azimuth_time_interval'])
            original_height = int(info['original_height'])
            original_stop = original_start + timedelta(
                seconds=(original_height - 1) * dt)

            target_start = max(common_start, original_start)
            target_stop = min(common_stop, original_stop)
            if target_stop < target_start:
                raise RuntimeError(IW + " burst " + nn +
                                   " 与公共方位时间范围没有交集")

            # 与 ISCE merge 的时间到行号转换一致，使用最近整数行。
            start_float = (target_start - original_start).total_seconds() / dt
            stop_float = (target_stop - original_start).total_seconds() / dt
            new_yoff = max(0, min(original_height - 1,
                                  int(math.floor(start_float + 0.5))))
            new_last = max(new_yoff, min(original_height - 1,
                                         int(math.floor(stop_float + 0.5))))
            old_yoff = info['yoff']
            old_height = info['height']
            info['yoff'] = str(new_yoff)
            info['height'] = str(new_last - new_yoff + 1)

            actual_start = original_start + timedelta(seconds=new_yoff * dt)
            actual_stop = original_start + timedelta(seconds=new_last * dt)
            new_intervals.append((actual_start, actual_stop, dt))
            print("    " + IW + " burst " + nn + ": yoff " + old_yoff +
                  " -> " + info['yoff'] + ", height " + old_height +
                  " -> " + info['height'])

        iw_start = min(item[0] for item in new_intervals)
        iw_stop = max(item[1] for item in new_intervals)
        dt = min(item[2] for item in new_intervals)
        if (abs((iw_start - common_start).total_seconds()) > 0.51 * dt or
                abs((iw_stop - common_stop).total_seconds()) > 0.51 * dt):
            raise RuntimeError(
                IW + " 的原始数据不能完整覆盖公共方位时间范围；无法生成无缺口矩形")


def _kml_coordinates(points: list) -> str:
    """将经纬度点列转换为闭合 KML coordinates。"""
    if not points:
        return ''
    closed = points + [points[0]]
    return ' '.join('%.8f,%.8f,0' % (lon, lat) for lon, lat in closed)


def _read_crop_kml_footprint(path: str) -> list:
    """读取 crop_rdr.py 已写出的独立 KML 中的裁剪窗口足迹。"""
    if not os.path.exists(path):
        raise FileNotFoundError("未找到独立裁剪 KML: " + path)
    tree = etree.parse(path)
    ns = '{http://www.opengis.net/kml/2.2}'
    placemarks = tree.getroot().findall('.//' + ns + 'Placemark')
    if len(placemarks) < 2:
        raise RuntimeError(path + " 中未找到裁剪窗口 Placemark")
    coordinates = placemarks[-1].find('.//' + ns + 'coordinates')
    if coordinates is None or not coordinates.text:
        raise RuntimeError(path + " 中未找到裁剪窗口坐标")
    points = []
    for value in coordinates.text.split():
        fields = value.split(',')
        if len(fields) >= 2:
            points.append((float(fields[0]), float(fields[1])))
    if len(points) > 1 and points[-1] == points[0]:
        points.pop()
    if len(points) < 3:
        raise RuntimeError(path + " 中的裁剪窗口坐标无效")
    return points


def write_crop_summary_kml(crop_info_by_iw: dict, lat1, lat2, lon1, lon2,
                           geom_root: str = './geom_reference') -> None:
    """汇总各独立 KML；绝不改写 crop_rdr.py 已生成的单独文件。"""
    kml_ns = 'http://www.opengis.net/kml/2.2'
    ns = '{' + kml_ns + '}'
    colors = ['ff00ff00', 'ffffff00', 'ffff00ff', 'ff00a5ff', 'ffffaa00']

    def add_style(document, style_id, color, width):
        style = etree.SubElement(document, ns + 'Style', id=style_id)
        line = etree.SubElement(style, ns + 'LineStyle')
        etree.SubElement(line, ns + 'color').text = color
        etree.SubElement(line, ns + 'width').text = str(width)
        poly = etree.SubElement(style, ns + 'PolyStyle')
        etree.SubElement(poly, ns + 'fill').text = '0'

    def add_polygon(parent, name, style_id, points):
        placemark = etree.SubElement(parent, ns + 'Placemark')
        etree.SubElement(placemark, ns + 'name').text = name
        etree.SubElement(placemark, ns + 'styleUrl').text = '#' + style_id
        polygon = etree.SubElement(placemark, ns + 'Polygon')
        etree.SubElement(polygon, ns + 'altitudeMode').text = 'clampToGround'
        outer = etree.SubElement(polygon, ns + 'outerBoundaryIs')
        ring = etree.SubElement(outer, ns + 'LinearRing')
        etree.SubElement(ring, ns + 'coordinates').text = _kml_coordinates(points)

    south, north = sorted((float(lat1), float(lat2)))
    west, east = sorted((float(lon1), float(lon2)))
    input_box = [(west, south), (east, south), (east, north), (west, north)]

    footprints = []
    for iw_index, IW in enumerate(sorted(crop_info_by_iw)):
        for nn in sorted(crop_info_by_iw[IW]):
            individual_path = os.path.join(
                geom_root, 'crop_extent_' + IW + '_' + nn + '.kml')
            # 该文件是在每个 burst 自己的原始 lat/lon 上计算窗口时生成的。
            # 后续为对齐 IW 而扩展物理窗口，不应反过来覆盖这个独立范围。
            points = _read_crop_kml_footprint(individual_path)
            style_id = 'crop_' + IW
            footprints.append((IW, nn, style_id, colors[iw_index % len(colors)], points))

    root = etree.Element(ns + 'kml', nsmap={None: kml_ns})
    document = etree.SubElement(root, ns + 'Document')
    etree.SubElement(document, ns + 'name').text = 'all_crop_extents'
    add_style(document, 'inputBox', 'ff0000ff', 3)
    add_polygon(document, '输入范围', 'inputBox', input_box)
    added_styles = set()
    for IW, nn, style_id, color, points in footprints:
        if style_id not in added_styles:
            add_style(document, style_id, color, 3)
            added_styles.add(style_id)
        add_polygon(document, IW + ' burst ' + nn + ' 裁剪窗口实际范围',
                    style_id, points)
    combined_path = os.path.join(geom_root, 'crop_extent.kml')
    etree.ElementTree(root).write(
        combined_path, encoding='UTF-8', xml_declaration=True,
        pretty_print=True)
    print(">>> 已汇总全部 IW/burst 独立裁剪 KML: " + combined_path)


def get_product_acquisition_date(product_root: str) -> str:
    """Read YYYYMMDD from the first retained burst sensingStart."""
    for xml_path in sorted(glob.glob(os.path.join(product_root, 'IW*.xml'))):
        tree = etree.parse(xml_path)
        bursts = tree.getroot().find(".//component[@name='bursts']")
        if bursts is None:
            continue
        for comp in bursts:
            name = comp.get('name') or ''
            if (isinstance(comp.tag, str) and comp.tag == 'component'
                    and name.startswith('burst')):
                value = _get_prop(comp, 'sensingstart')
                match = re.match(r'(\d{4})-(\d{2})-(\d{2})', value or '')
                if match:
                    return ''.join(match.groups())
    match = re.search(r'(\d{8})', os.path.basename(os.path.normpath(product_root)))
    if match:
        return match.group(1)
    raise RuntimeError('无法确定产品日期: ' + product_root)


def save_crop_metadata(crop_info_by_iw: dict, acquisitions: dict,
                       path: str = './geom_reference/crop_metadata.json') -> None:
    """Save per-acquisition full-burst carrier origins for mburst resampling."""
    payload = {
        'version': 2,
        'description': 'Per-acquisition TOPS full-burst grid used after physical crop',
        'reference_swaths': crop_info_by_iw,
        'acquisitions': acquisitions,
    }
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
    print(">>> 已保存各日期 TOPS 原始相位参考信息: " + path)


def _shift_datetime_str(s: str, dt_seconds: float, n_lines: int) -> str:
    """将 ISO 格式日期时间字符串平移 n_lines * dt_seconds，保持原格式。"""
    fmt = '%Y-%m-%d %H:%M:%S.%f'
    d = datetime.strptime(s.strip(), fmt)
    d += timedelta(seconds=n_lines * dt_seconds)
    return d.strftime(fmt)


def update_burst_slc_xml(xml_path: str, width: int, height: int) -> None:
    """更新单个 burst_*.slc.xml / burst_*.int.xml 的 width/length 及其 coordinate。"""
    if not os.path.exists(xml_path):
        return
    tree = etree.parse(xml_path)
    root = tree.getroot()
    _set_prop(root, 'width', str(width))
    _set_prop(root, 'length', str(height))
    _set_prop(root, 'xmax', str(width))
    for coord in root.findall("component[@name='coordinate1']"):
        _set_prop(coord, 'size', str(width))
        _set_prop(coord, 'endingvalue', str(width))
    for coord in root.findall("component[@name='coordinate2']"):
        _set_prop(coord, 'size', str(height))
        _set_prop(coord, 'endingvalue', str(height))
    tree.write(xml_path, pretty_print=True)


def crop_reference_burst_slc(IW: str, crop_info: dict, slc_dir: str = None) -> None:
    """按裁剪窗口修正 <slc_dir>/burst_*.slc.vrt（或实际 .slc 数据）及其 xml。
    这是多 burst 正确合并的关键：只有 geom 裁剪而 SLC 仍按左上角截取，
    会导致 merge 时各 burst 不在同一地理范围、出现黑带。
    slc_dir 默认为 './reference/IW'，可传入 './secondarys/DATE/IW' 以裁剪从影像。"""
    if slc_dir is None:
        slc_dir = './reference/' + IW
    iw_ref_dir = slc_dir
    if not os.path.isdir(iw_ref_dir):
        print("警告：未找到 SLC 目录 " + iw_ref_dir + "，跳过 SLC 裁剪")
        return
    for nn, info in sorted(crop_info.items()):
        xoff = int(info['xoff'])
        yoff = int(info['yoff'])
        width = int(info['width'])
        height = int(info['height'])
        base = os.path.join(iw_ref_dir, 'burst_' + nn)
        vrt_path = base + '.slc.vrt'
        xml_path = base + '.slc.xml'
        slc_path = base + '.slc'

        if os.path.exists(vrt_path):
            try:
                tree = etree.parse(vrt_path)
                root = tree.getroot()
                old_rx = int(root.get('rasterXSize', '0'))
                old_ry = int(root.get('rasterYSize', '0'))
                # 若裁剪窗口与原始画布完全一致，则无需修改 VRT
                if xoff == 0 and yoff == 0 and width == old_rx and height == old_ry:
                    print("    跳过 VRT: " + vrt_path + "（裁剪窗口与原始尺寸一致）")
                    continue

                # 新画布尺寸
                root.set('rasterXSize', str(width))
                root.set('rasterYSize', str(height))

                # 对每个 SimpleSource，计算 desired output 与 available output 的交集，
                # 再映射回 source 坐标，避免读取源文件越界。
                for src in root.findall('.//SimpleSource'):
                    src_rect = src.find('SrcRect')
                    dst_rect = src.find('DstRect')
                    if src_rect is None or dst_rect is None:
                        continue
                    sx0 = int(src_rect.get('xOff', '0'))
                    sy0 = int(src_rect.get('yOff', '0'))
                    sx1 = sx0 + int(src_rect.get('xSize', '0'))
                    sy1 = sy0 + int(src_rect.get('ySize', '0'))
                    dx0 = int(dst_rect.get('xOff', '0'))
                    dy0 = int(dst_rect.get('yOff', '0'))
                    dx1 = dx0 + int(dst_rect.get('xSize', '0'))
                    dy1 = dy0 + int(dst_rect.get('ySize', '0'))

                    # 期望的输出区域（在原始 burst 坐标系中）
                    ox0, oy0 = xoff, yoff
                    ox1, oy1 = xoff + width, yoff + height

                    # 可用输出区域与期望输出区域的交集
                    nx0 = max(dx0, ox0)
                    ny0 = max(dy0, oy0)
                    nx1 = min(dx1, ox1)
                    ny1 = min(dy1, oy1)
                    if nx1 <= nx0 or ny1 <= ny0:
                        # 无交集：将该 source 置空
                        src_rect.set('xSize', '0')
                        src_rect.set('ySize', '0')
                        dst_rect.set('xSize', '0')
                        dst_rect.set('ySize', '0')
                        continue

                    # 输出坐标 -> source 坐标：source = src_origin + (output - dst_origin)
                    nsx0 = sx0 + (nx0 - dx0)
                    nsy0 = sy0 + (ny0 - dy0)
                    nsx1 = sx0 + (nx1 - dx0)
                    nsy1 = sy0 + (ny1 - dy0)

                    # 写回新的 SrcRect / DstRect（新输出坐标系下以 (0,0) 为起点）
                    src_rect.set('xOff', str(nsx0))
                    src_rect.set('yOff', str(nsy0))
                    src_rect.set('xSize', str(nsx1 - nsx0))
                    src_rect.set('ySize', str(nsy1 - nsy0))
                    dst_rect.set('xOff', str(nx0 - ox0))
                    dst_rect.set('yOff', str(ny0 - oy0))
                    dst_rect.set('xSize', str(nx1 - nx0))
                    dst_rect.set('ySize', str(ny1 - ny0))

                tree.write(vrt_path, pretty_print=True)
                print("    已裁剪 VRT: " + vrt_path +
                      " (xoff=" + str(xoff) + " yoff=" + str(yoff) +
                      " w=" + str(width) + " h=" + str(height) + ")")
            except Exception as e:
                print("警告：裁剪 VRT " + vrt_path + " 失败: " + str(e))
        elif os.path.exists(slc_path) and os.path.getsize(slc_path) > 0:
            # 存在实际数据文件时，用 gdal_translate 按窗口重新裁剪
            tmp_slc = slc_path + '.orig'
            os.system('mv ' + slc_path + ' ' + tmp_slc)
            os.system('gdal_translate -srcwin ' + str(xoff) + ' ' + str(yoff) + ' ' +
                      str(width) + ' ' + str(height) +
                      ' -of envi -co INTERLEAVE=BIP ' + tmp_slc + ' ' + slc_path)
            os.system('rm -f ' + tmp_slc + ' ' + tmp_slc + '.hdr')
            print("    已裁剪 SLC 数据: " + slc_path)

        # 同步更新 burst_*.slc.xml 中的尺寸信息
        update_burst_slc_xml(xml_path, width, height)


def update_iw_xml_geometry(xml_path: str, crop_info: dict) -> None:
    """依据各 burst 的裁剪窗口，更新 IWn.xml 中每个 burst 组件的尺寸与几何位置。
    crop_info: {burst序号(两位): {'xoff','yoff','width','height'}}，
    如 {'01':{'xoff':'520','yoff':'19','width':'2298','height':'465'}}"""
    if not os.path.exists(xml_path):
        print("警告：未找到 " + xml_path + "，跳过 xml 更新")
        return
    tree = etree.parse(xml_path)
    root = tree.getroot()
    bursts_comp = root.find(".//component[@name='bursts']")
    if bursts_comp is None:
        print("警告：" + xml_path + " 中未找到 bursts 组件，跳过 xml 更新")
        return
    changed = 0
    for comp in bursts_comp:
        if not (isinstance(comp.tag, str) and comp.tag == 'component'
                and (comp.get('name') or '').startswith('burst')):
            continue
        bnum = comp.get('name')[5:]          # 'burst1' -> '1'
        nn = bnum.zfill(2)
        if nn not in crop_info:
            print("警告：" + xml_path + " 中 " + comp.get('name') + " 无对应裁剪信息，跳过")
            continue
        info = crop_info[nn]
        xoff = int(info['xoff'])
        yoff = int(info['yoff'])
        width = int(info['width'])
        height = int(info['height'])

        # 读取旧尺寸与起始可用像素
        old_w = int(_get_prop(comp, 'numberofsamples', '0') or '0')
        old_h = int(_get_prop(comp, 'numberoflines', '0') or '0')
        fvl_old = int(_get_prop(comp, 'firstvalidline', '0') or '0')
        fvs_old = int(_get_prop(comp, 'firstvalidsample', '0') or '0')
        nvl_old = int(_get_prop(comp, 'numberofvalidlines', '0') or '0')
        nvs_old = int(_get_prop(comp, 'numberofvalidsamples', '0') or '0')

        # 若裁剪窗口与原始画布完全一致，则无需修改 xml
        if xoff == 0 and yoff == 0 and width == old_w and height == old_h:
            print("    跳过 " + comp.get('name') + " xml 更新（裁剪窗口与原始尺寸一致）")
            continue

        # 新文件中的有效区域：旧有效区域平移 (-yoff, -xoff) 后再与新画布求交
        new_fvl = max(0, fvl_old - yoff)
        new_fvs = max(0, fvs_old - xoff)
        old_vl_end = fvl_old + nvl_old - 1
        old_vs_end = fvs_old + nvs_old - 1
        new_vl_end = min(old_vl_end - yoff, height - 1)
        new_vs_end = min(old_vs_end - xoff, width - 1)
        new_nvl = max(0, new_vl_end - new_fvl + 1)
        new_nvs = max(0, new_vs_end - new_fvs + 1)

        # 1) 更新尺寸
        _set_prop(comp, 'firstvalidline', str(new_fvl))
        _set_prop(comp, 'firstvalidsample', str(new_fvs))
        _set_prop(comp, 'numberoflines', str(height))
        _set_prop(comp, 'numberofvalidlines', str(new_nvl))
        _set_prop(comp, 'numberofsamples', str(width))
        _set_prop(comp, 'numberofvalidsamples', str(new_nvs))
        img = comp.find("component[@name='image']")
        if img is not None:
            _set_prop(img, 'width', str(width))
            _set_prop(img, 'length', str(height))
            _set_prop(img, 'xmax', str(width))
            for coord in img.findall("component[@name='coordinate1']"):
                _set_prop(coord, 'size', str(width))
                _set_prop(coord, 'endingvalue', str(width))
            for coord in img.findall("component[@name='coordinate2']"):
                _set_prop(coord, 'size', str(height))
                _set_prop(coord, 'endingvalue', str(height))

        # 2) 更新时间/距离几何：sensingstart/startingrange 描述文件第 0 行/列；
        #    burststartutc/burststoputc 为 annotation 原始 burst 时间，不应修改。
        dt = float(_get_prop(comp, 'azimuthtimeinterval', '0'))
        dr = float(_get_prop(comp, 'rangepixelsize', '0'))
        if dt > 0:
            ss_old = _get_prop(comp, 'sensingstart')
            if ss_old:
                _set_prop(comp, 'sensingstart',
                          _shift_datetime_str(ss_old, dt, yoff))
                _set_prop(comp, 'sensingstop',
                          _shift_datetime_str(ss_old, dt, yoff + height - 1))
        if dr > 0:
            sr_old = _get_prop(comp, 'startingrange')
            if sr_old:
                _set_prop(comp, 'startingrange',
                          str(float(sr_old) + xoff * dr))

        print("    " + comp.get('name') + ": xoff=" + str(xoff) + " yoff=" + str(yoff) +
              " width=" + str(width) + " height=" + str(height))
        changed += 1
    if changed:
        tree.write(xml_path, pretty_print=True)
        print(">>> 已更新 " + xml_path + " 中 " + str(changed) + " 个 burst 的尺寸与几何")

dem_name="./dem/dem.tif"
dem_path = Path('./dem/full_res.dem.wgs84')

####dem_path绝对路径
dem_path=dem_path.resolve()
gdal.Translate(str(dem_path), dem_name, format='ISCE')
xml_path = tag_dem_xml_as_ellipsoidal(dem_path)
fix_image_xml(xml_path)




#os.chdir('../')
path = "./run_files"
if os.path.exists(path):
    print("run_files were existed")
else:
    # 保持 geometry 配准，不启用 NESD。
    os.system("python ../code/topsStack/stackSentinel.py -s ./slc -d ./dem/full_res.dem.wgs84  -a ./aux_cal/ -o ./orbits -C geometry -c "+num_connect+" -z "+num_azimuth+" -r "+num_range+" -f "+filter_strength)
    line = open('./run_files/run_01_unpack_topo_reference').readline()
    line = line.replace("\n", "")
    line = str(line) + '0'
    print(line)
    open('./run_files/run_00_unpack_topo_reference', 'w').writelines(str(line))
    lines = open('./configs/config_reference').readlines()
    open('./configs/config_reference0', 'w').writelines(lines[0:12])
    open('./configs/config_reference', 'w').writelines(lines[0:4] + lines[15:])
    ######删除主影像合并到merged的重复命令（merge 步骤编号随 coregistration 方式变化，动态匹配）
    merge_run_files = sorted(glob.glob('./run_files/run_*_merge_reference_secondary_slc'))
    for merge_run in merge_run_files:
        line_merged = open(merge_run).readlines()
        open(merge_run, 'w').writelines(line_merged[1:])

def configure_mburst_resampler():
    """Make only this workflow call the crop-aware resampler module."""
    changed = 0
    for config_path in sorted(glob.glob('./configs/config_*resample*')):
        text = Path(config_path).read_text(encoding='utf-8')
        updated = re.sub(r'(?m)^resamp_withCarrier\s*:',
                         'resamp_withCarrier_mburst :', text)
        if updated != text:
            Path(config_path).write_text(updated, encoding='utf-8')
            changed += 1
    print('>>> mburst1 已启用原始网格载波重采样配置: ' + str(changed) + ' 个')


configure_mburst_resampler()

# 存储 reference 各 IW 的裁剪窗口，供后续对 secondarys 复用同一套 xoff/yoff
REF_CROP_INFO = {}
# 各日期独立的完整 burst 元数据；不能把 reference 高度/起始距离套给 secondary。
ACQUISITION_CROP_INFO = {}
# True 表示 reference VRT/XML 已在 run_01 topo 前完成物理裁剪。
REFERENCE_PRETOPO_CROP = False
# True 表示提供了完整经纬度框，且已在 topo 后完成二次裁剪。
REFERENCE_POSTTOPO_CROP = False

# 动态发现 run_files 下所有步骤（按名称匹配，不依赖硬编码编号；切换 coregistration 方式时步骤数会变）
run = sorted([os.path.basename(p) for p in glob.glob('./run_files/run_*')])
if len(run) == 0:
    raise RuntimeError("未在 ./run_files 下发现任何 run_* 文件，请检查 stackSentinel.py 是否成功生成 run_files")




# 要运行的步骤

for i in range(len(run)):
    ##########
    print("正在进行步骤：" + str(run[i]))
    runstep = run[i]
    if runstep == 'run_01_unpack_topo_reference':
        os.system('rm -rf ./geom_reference')
        # 无论是否提供经纬度框，都在 topo 前按「已下载 SLC 非零窗口（无黑边区）」裁剪
        # reference，使 topo 只处理有效（无黑边）数据、避免全幅/黑边计算。经纬度框的
        # 精确裁剪交给下方 lat/lon 块在 topo 后统一执行（crop_reference_burst_slc
        # 偏移安全，可在此基础上再裁到框）。
        # 注意：不能用 REF_CROP_INFO/经纬度是否提供来判断是否需裁剪——REF_CROP_INFO
        # 会从 crop_state.json 跨进程恢复，若残留旧状态而磁盘 reference 已被重新解包
        # 成全幅，会错误地跳过裁剪，导致 topo 处理整幅（含黑边）。改用标记文件判断
        # 磁盘 reference 是否已完成无黑边裁剪，缺失或元数据丢失时一律重裁。
        # 每次根据当前磁盘上的 VRT 重建窗口，避免旧标记或旧状态
        # 使本次 topo 误用全幅/过期几何。已裁剪数据再检测时偏移为 0，是幂等的。
        prepare_reference_crop_before_topo()
    with open('./run_files/' + str(runstep), "r") as f:
        a = f.readlines()
        print(str(a[0])[0:6])
    if str(a[0])[0:6] != "python":
        for i in range(len(a)):
            a[i] = "python ../code/topsStack/" + a[i]
        with open('./run_files/' + str(runstep), "w") as f:
            for i in range(len(a)):
                f.write(str(a[i]))

    if unw_mp =="":
        unw_mp=mp
    # generate_burst_igram 步骤使用解缠专用并行数
    if 'generate_burst_igram' in runstep:
        flow_mp=unw_mp
    else:
        flow_mp=mp

    # extractCommonValidRegion.py 会复用已有 stack/IW*.xml。若该目录来自
    # 未裁剪或旧参数的运行，其 firstValidSample/firstValidLine 会把大黑边
    # 和错误有效区带回新的裁剪结果。必须在 run_06 启动前清掉旧 stack。
    if ('extract_stack_valid_region' in runstep and
            HAS_GEO_BBOX):
        if os.path.isdir('./stack'):
            print(">>> 检测到旧 stack，删除后按当前裁剪后的 reference 重新计算公共有效区")
            shutil.rmtree('./stack')

    if mode=="SLC":
        # SLC 模式只跑到 merge_reference_secondary_slc（含）为止
        if 'merge_reference_secondary_slc' in runstep:
            subprocess.run([sys.executable, '../code/topsStack/run.py',
                            '-i', './run_files/' + str(runstep),
                            '-p', str(flow_mp)], check=True)
            break
        else:
            subprocess.run([sys.executable, '../code/topsStack/run.py',
                            '-i', './run_files/' + str(runstep),
                            '-p', str(flow_mp)], check=True)
    else:
        subprocess.run([sys.executable, '../code/topsStack/run.py',
                        '-i', './run_files/' + str(runstep),
                        '-p', str(flow_mp)], check=True)

    ##############################################################################
    # 更改干涉处理中主影像数据路径
    if runstep == 'run_00_unpack_topo_reference':
        pairlist = glob.glob('./configs/config_generate_igram*')
        # print(pairlist)
        for i in range(len(pairlist)):
            lines = open(pairlist[i]).readlines()
            if lines[5][-10:-1] == "reference":
                lines[5] = lines[5][:-10] + 'coreg_secondarys/' + pairlist[i][-17:-9] + lines[5][-1]
                open(pairlist[i], 'w').writelines(lines)

    ##############################################################################
    # ↓↓↓ 多burst裁剪处理：对 geom_reference 下每个 IW（swath），再对每个 burst ↓↓↓
    # （lat_01 / lat_02 / ...）分别按其自身 lat/lon 计算窗口并裁剪其全部 geom 文件，
    # 随后在 run_03 阶段按各 burst 真实尺寸更新 reference/IWn.xml。
    # 若处理整个 burst（不裁剪），将下方 if 整段注释即可
    ##############################################################################
    if HAS_GEO_BBOX:
        # ---------------- 步骤 run_01：裁剪 geom_reference + reference SLC VRT ----------------
        # 物理裁剪数据，使后续所有步骤都基于裁剪后的参考网格。
        # IWn.xml 的更新推迟到 run_03，与 secondarys 一起完成，避免 run_02 解包前
        # reference xml 已被改写而 secondary xml 仍为原始值。
        if runstep == 'run_01_unpack_topo_reference' and not REFERENCE_POSTTOPO_CROP:
            # topo 前的 REF_CROP_INFO 是非零数据窗；从这里起必须换成
            # 以 topo 产生的 lat/lon 为基准计算的二次裁剪窗口。
            REF_CROP_INFO.clear()
            iw_list = get_iw_list()
            if len(iw_list) == 0:
                print("警告：geom_reference 下未检测到 IW* 目录，跳过裁剪")
            # 第一阶段：只计算所有 IW/burst 的原始窗口并保存完整 burst 元数据。
            for IW in iw_list:
                crop_info = calculate_iw_crop_info(IW, lat1, lat2, lon1, lon2)
                if not crop_info:
                    # 整个 IW 不覆盖目标范围：从 reference 与 geometry 中停用，
                    # 使后续 getSwathList/merge 不会再次发现并处理它。
                    prune_product_bursts('./reference', IW, {})
                    skipped_geom = os.path.join('./geom_reference', IW)
                    if os.path.isdir(skipped_geom):
                        shutil.rmtree(skipped_geom)
                    continue
                # 必须在 IW*.xml 被改写前保存完整 burst 的中心行和起始距离。
                # TOPS 载频不是以裁剪后影像中心为参考；丢失这些值会在 burst
                # 交界处产生相位常数/斜坡跳变。
                add_original_burst_metadata('./reference/' + IW + '.xml', crop_info)
                # 存下该 IW 的裁剪窗口，供 secondarys 复用（保证主从影像 burst 几何对齐）
                REF_CROP_INFO[IW] = crop_info

            if not REF_CROP_INFO:
                raise RuntimeError("裁剪范围与 reference 的所有 IW/burst 均无交集")
            # 第二阶段：跨 IW 统一绝对方位时间，然后才物理裁剪 geometry/SLC。
            align_iw_azimuth_windows(REF_CROP_INFO)
            reference_date = get_product_acquisition_date('./reference')
            ACQUISITION_CROP_INFO[reference_date] = copy.deepcopy(REF_CROP_INFO)
            print(">>> 已保存 reference 日期原始载波网格: " + reference_date)
            for IW, crop_info in REF_CROP_INFO.items():
                print(">>> 正在裁剪 IW [" + IW + "] 的 geometry...")
                apply_iw_geom_crop(IW, crop_info)
                print(">>> 正在更新 IW [" + IW + "] 的 reference SLC VRT...")
                crop_reference_burst_slc(IW, crop_info)
                prune_product_bursts('./reference', IW, crop_info)
            if REF_CROP_INFO:
                # 独立 KML 已在各 burst 计算窗口时写出。这里只汇总它们，
                # 不使用跨 IW 对齐后扩展的物理范围覆盖原始裁剪窗口。
                write_crop_summary_kml(REF_CROP_INFO, lat1, lat2, lon1, lon2)
                save_crop_metadata(REF_CROP_INFO, ACQUISITION_CROP_INFO)
            REFERENCE_POSTTOPO_CROP = True

        # ---------------- 步骤 run_02：secondarys 解包完成后，与主影像同步裁剪 SLC VRT ----------------
        # 关键修复：若仅裁剪 reference 而 secondarys 保持原始几何，则主从 burst 的
        # sensingStart / startingRange 错位，geo2rdr 的 getBurstOffset 找不到合法偏移，
        # 从 run_04 起级联崩溃。此处用与 reference 完全相同的 xoff/yoff 裁剪每个
        # secondary 的 SLC VRT，xml 几何平移放到 run_03 统一完成。
        if runstep == 'run_02_unpack_secondary_slc':
            if not REF_CROP_INFO:
                print("警告：未获取到 reference 裁剪窗口，跳过 secondarys 裁剪")
            else:
                sec_dirs = sorted(glob.glob('./secondarys/*/'))
                if len(sec_dirs) == 0:
                    print("警告：未检测到 secondarys/* 目录，跳过 secondarys 裁剪")
                for sec_dir in sec_dirs:
                    sec_date = sec_dir.rstrip('/').split('/')[-1]
                    print(">>> 正在对 secondary [" + sec_date + "] 执行与 reference 同步的多 burst VRT 裁剪...")
                    secondary_phase_info = {}
                    # secondary 解包仍可能生成所有原始 IW；这里仅保留 reference
                    # 实际覆盖目标范围的 IW/burst。
                    for IW in get_product_iw_list(sec_dir):
                        if IW not in REF_CROP_INFO:
                            prune_product_bursts(sec_dir, IW, {})
                            continue
                        sec_iw_xml = sec_dir + IW + '.xml'
                        if not os.path.exists(sec_iw_xml):
                            print("    跳过 " + sec_iw_xml + "（不存在）")
                            continue
                        ref_ci = REF_CROP_INFO[IW]
                        sec_ci = build_product_nonzero_crop_info(
                            sec_dir, IW, allowed_bursts=ref_ci.keys())
                        if not sec_ci:
                            print("    警告：secondary " + sec_date + " " + IW +
                                  " 没有非零在线窗口，跳过")
                            prune_product_bursts(sec_dir, IW, {})
                            continue
                        add_original_burst_metadata(sec_iw_xml, sec_ci)
                        secondary_phase_info[IW] = sec_ci
                        for nn, info in sorted(sec_ci.items()):
                            print("    secondary 自身窗口 " + IW + " burst " + nn +
                                  ": xoff=" + str(info['xoff']) +
                                  " yoff=" + str(info['yoff']) +
                                  " width=" + str(info['width']) +
                                  " height=" + str(info['height']))
                        try:
                            crop_reference_burst_slc(IW, sec_ci, slc_dir=sec_dir + IW)
                            prune_product_bursts(sec_dir, IW, sec_ci)
                        except Exception as e:
                            print("    警告：裁剪 secondary " + sec_date + " " + IW + " VRT 失败: " + str(e))
                    if secondary_phase_info:
                        ACQUISITION_CROP_INFO[sec_date.replace('-', '')] = secondary_phase_info
                        print("    已保存 secondary 原始载波网格: " + sec_date)


        # ---------------- 步骤 run_03：统一更新 reference 与 secondarys 的 IWn.xml 几何 ----------------
        # 放在 baseline 计算之后、geo2rdr 之前，确保 xml 中 numberOfLines/Samples、
        # sensingStart、startingRange 等参数与裁剪后的 VRT 严格一致。
        if runstep == 'run_03_average_baseline':
            if not REF_CROP_INFO:
                print("警告：未获取到 reference 裁剪窗口，跳过 xml 几何更新")
            else:
                save_crop_metadata(REF_CROP_INFO, ACQUISITION_CROP_INFO)
                for IW, ci in REF_CROP_INFO.items():
                    ref_xml = 'reference/' + IW + '.xml'
                    if os.path.exists(ref_xml):
                        if REFERENCE_PRETOPO_CROP and not REFERENCE_POSTTOPO_CROP:
                            print(">>> reference " + IW + ".xml 已在 topo 前更新，跳过重复平移")
                        else:
                            print(">>> 正在更新 reference " + IW + ".xml 几何...")
                            update_iw_xml_geometry(ref_xml, ci)
                    else:
                        print("    跳过 " + ref_xml + "（不存在）")
                    sec_dirs = sorted(glob.glob('./secondarys/*/'))
                    for sec_dir in sec_dirs:
                        sec_date = sec_dir.rstrip('/').split('/')[-1]
                        sec_xml = sec_dir + IW + '.xml'
                        if os.path.exists(sec_xml):
                            date_key = sec_date.replace('-', '')
                            sec_ci = ACQUISITION_CROP_INFO.get(date_key, {}).get(IW)
                            if not sec_ci:
                                raise RuntimeError(
                                    "缺少 secondary " + sec_date + " " + IW +
                                    " 的独立裁剪窗口")
                            print("    正在按自身窗口更新 secondary " + sec_date +
                                  " " + IW + ".xml 几何...")
                            update_iw_xml_geometry(sec_xml, sec_ci)
                        else:
                            print("    跳过 " + sec_xml + "（不存在）")
