#!/usr/bin/env python3
# ============================================================================
#by: ZYF-SDUST
# 本版改用 "在 lat/lon rdr (每个雷达像元的经纬度) 中直接检索落在框内的像元,
# 取其行列极值" 的精确栅格查找, 消除上述偏差。
#
# 用法 :
#   python crop_rdr.py -b 'S N W E' -lat lat_01.rdr -lon lon_01.rdr [-dfac N]
# ============================================================================

import sys
from osgeo import gdal
import argparse
import numpy as np


def cmdLineParse():
    parser = argparse.ArgumentParser(
        description='Generate the gdal_translate -srcwin command to crop RDR data '
                    'based on a lon-lat BBox (exact raster lookup, no plane fit).')
    parser.add_argument('-i', '--input', dest='file', type=str, required=False,
                        help='Input filename (GDAL supported)')
    parser.add_argument('-b', '--bbox', dest='bbox', type=str, required=True,
                        help='Lat/Lon Bounding SNWE: "S N W E"')
    parser.add_argument('-dfac', '--downfac', dest='down_sample', type=str, required=False,
                        default='3', help='Coarse scan step (speed only); final box is exact on full-res.')
    parser.add_argument('-nd', '--nodata', dest='nodata', type=str, required=False,
                        default='0', help='Lon/Lat no-data value')
    parser.add_argument('-lat', '--lat', dest='latfile', type=str, required=False,
                        default='lat_01.rdr', help='Lat filename')
    parser.add_argument('-lon', '--lon', dest='lonfile', type=str, required=False,
                        default='lon_01.rdr', help='Lon filename')
    parser.add_argument('-kml', '--kml', dest='kml', type=str, required=False,
                        default="./geom_reference/crop_extent.kml", help='把裁剪窗口四角经纬度写成 KML (如 crop_extent.kml) 便于可视化核对')
    return parser.parse_args()


def _write_kml(path, S, N, W, E, Lat_crop, Lon_crop):
    """写 KML: 输入框(S,N,W,E,红) + 裁剪窗口实际范围(绿)。"""
    import xml.etree.ElementTree as ET
    mask = np.isfinite(Lat_crop) & np.isfinite(Lon_crop)
    if not mask.any():
        return
    lat_v = Lat_crop[mask]; lon_v = Lon_crop[mask]
    clat_min, clat_max = float(lat_v.min()), float(lat_v.max())
    clon_min, clon_max = float(lon_v.min()), float(lon_v.max())
    # 四角足迹
    ys, xs = np.where(mask)
    y0, y1 = int(ys.min()), int(ys.max())
    x0, x1 = int(xs.min()), int(xs.max())
    foot = [(float(Lon_crop[y0, x0]), float(Lat_crop[y0, x0])),
            (float(Lon_crop[y0, x1]), float(Lat_crop[y0, x1])),
            (float(Lon_crop[y1, x1]), float(Lat_crop[y1, x1])),
            (float(Lon_crop[y1, x0]), float(Lat_crop[y1, x0]))]
    cb = [(clon_min, clat_min), (clon_max, clat_min), (clon_max, clat_max), (clon_min, clat_max)]
    ib = [(W, S), (E, S), (E, N), (W, N)]

    def poly(name, style, coords):
        pts = " ".join("%.8f,%.8f,0" % (p[0], p[1]) for p in coords + [coords[0]])
        return ('    <Placemark>\n      <name>%s</name>\n      <styleUrl>#%s</styleUrl>\n'
                '      <Polygon><altitudeMode>clampToGround</altitudeMode>'
                '<outerBoundaryIs><LinearRing><coordinates>%s</coordinates>'
                '</LinearRing></outerBoundaryIs></Polygon>\n    </Placemark>\n') % (name, style, pts)

    kml = ['<?xml version="1.0" encoding="UTF-8"?>']
    kml.append('<kml xmlns="http://www.opengis.net/kml/2.2">')
    kml.append('  <Document><name>crop_extent</name>')
    kml.append('    <Style id="inputBox"><LineStyle><color>ff0000ff</color><width>3</width></LineStyle>'
               '<PolyStyle><fill>0</fill></PolyStyle></Style>')
    kml.append('    <Style id="cropFoot"><LineStyle><color>ff00ff00</color><width>3</width></LineStyle>'
               '<PolyStyle><fill>0</fill></PolyStyle></Style>')
    kml.append(poly('输入范围 (S,N,W,E)', 'inputBox', ib))
    kml.append(poly('裁剪窗口实际范围', 'cropFoot', foot))
    kml.append('  </Document></kml>')
    kml = "\n".join(kml)
    ET.fromstring(kml)  # 校验
    with open(path, 'w') as f:
        f.write(kml)
    print('[kml] 已写出 %s  (裁剪范围 lat[%.6f,%.6f] lon[%.6f,%.6f])'
          % (path, clat_min, clat_max, clon_min, clon_max), file=sys.stderr)


def _lookup_box(Lat, Lon, S, N, W, E):
    """在给定 Lat/Lon 数组中返回落在 [S,N]x[W,E] 内的像元极值 (0-based 全分辨率)。"""
    mask = (Lat >= S) & (Lat <= N) & (Lon >= W) & (Lon <= E)
    if not mask.any():
        return None
    rows, cols = np.where(mask)
    return (int(rows.min()), int(cols.min()),
            int(rows.max() - rows.min() + 1), int(cols.max() - cols.min() + 1))


if __name__ == '__main__':
    inps = cmdLineParse()
    ds = max(1, int(inps.down_sample))
    nodata = float(inps.nodata)
    bbox = np.fromstring(inps.bbox, dtype=float, sep=' ')
    S, N, W, E = bbox[0], bbox[1], bbox[2], bbox[3]

    LatData = gdal.Open(inps.latfile)
    LonData = gdal.Open(inps.lonfile)
    Lat = LatData.GetRasterBand(1).ReadAsArray().astype(np.float64)
    Lon = LonData.GetRasterBand(1).ReadAsArray().astype(np.float64)
    n_lines_full, n_pixels_full = Lon.shape

    # 处理 nodata (默认 0)
    Lat[Lat == nodata] = np.nan
    Lon[Lon == nodata] = np.nan

    # 1) 粗扫描快速定位大致区域
    Lat_c = Lat[::ds, ::ds]
    Lon_c = Lon[::ds, ::ds]
    box = _lookup_box(Lat_c, Lon_c, S, N, W, E)

    if box is None:
        # 粗采样下无命中 -> 直接全分辨率检索
        box = _lookup_box(Lat, Lon, S, N, W, E)
        if box is None:
            raise RuntimeError('BBox (S,N,W,E)=%s 不在 lat/lon 数据覆盖范围内' % ((S, N, W, E),))
        yoff, xoff, ysize, xsize = box
    else:
        rmin, cmin, _, _ = box
        # 粗采样 0-based 索引 i -> 全分辨率 0-based 起点 i*ds; 向外扩 1 个粗步长以确保覆盖
        y0 = max(0, rmin * ds - ds)
        y1 = min(n_lines_full, (rmin + box[2]) * ds + ds)
        x0 = max(0, cmin * ds - ds)
        x1 = min(n_pixels_full, (cmin + box[3]) * ds + ds)
        sub_box = _lookup_box(Lat[y0:y1, x0:x1], Lon[y0:y1, x0:x1], S, N, W, E)
        if sub_box is None:
            # 退化: 退回全分辨率全局检索
            gbox = _lookup_box(Lat, Lon, S, N, W, E)
            if gbox is None:
                raise RuntimeError('BBox 不在覆盖范围内')
            yoff, xoff, ysize, xsize = gbox
        else:
            yoff = y0 + sub_box[0]
            xoff = x0 + sub_box[1]
            ysize = sub_box[2]
            xsize = sub_box[3]

    # 边界保护
    yoff = max(0, yoff)
    xoff = max(0, xoff)
    if yoff + ysize > n_lines_full:
        ysize = n_lines_full - yoff
    if xoff + xsize > n_pixels_full:
        xsize = n_pixels_full - xoff

    print('gdal_translate -srcwin ' + str(xoff) + ' ' + str(yoff) + ' ' +
          str(xsize) + ' ' + str(ysize) + ' -of envi -co INTERLEAVE=BIP ')

    # 可选: 输出裁剪窗口四角经纬度 KML, 便于核对裁剪范围是否等于输入地理范围
    if inps.kml:
        Lat_crop = Lat[yoff:yoff + ysize, xoff:xoff + xsize]
        Lon_crop = Lon[yoff:yoff + ysize, xoff:xoff + xsize]
        _write_kml(inps.kml, S, N, W, E, Lat_crop, Lon_crop)
