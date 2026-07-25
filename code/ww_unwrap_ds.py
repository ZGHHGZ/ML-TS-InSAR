#!/usr/bin/env python3

import os
import time
import argparse
import subprocess
import numpy as np

from multiprocessing import Pool
from tqdm import tqdm



# =====================================================
# 命令行参数
# =====================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description="Parallel Whirlwind InSAR phase unwrapping"
    )

    parser.add_argument(
        "-p",
        "--process",
        type=int,
        default=8,
        help="并行进程数量"
    )

    parser.add_argument(
        "-n",
        "--nlooks",
        type=int,
        default=5,
        help="Whirlwind nlooks"
    )

    parser.add_argument(
        "-m",
        "--max-ncomps",
        type=int,
        default=255,
        help="Whirlwind 最大连通分量数（1字节 conncomp 上限 255）"
    )

    parser.add_argument(
        "-r",
        "--root",
        type=str,
        default="./merged/interferograms",
        help="干涉图目录"
    )

    parser.add_argument(
        "-w",
        "--whirlwind",
        type=str,
        default="whirlwind",
        help="whirlwind 可执行文件路径（默认从 PATH 查找）"
    )

    return parser.parse_args()



# =====================================================
# 通过 gdalinfo 读取影像尺寸（利用 .aux.xml）
# =====================================================

def get_size(filepath):

    result = subprocess.run(

        ["gdalinfo", filepath],

        stdout=subprocess.PIPE,

        stderr=subprocess.PIPE,

        text=True

    )


    if result.returncode != 0:

        raise RuntimeError(
            "gdalinfo 失败: " + result.stderr.strip()
        )


    for line in result.stdout.splitlines():

        if "Size is" in line:

            parts = line.split()

            width = int(parts[2].rstrip(","))

            length = int(parts[3])

            return length, width


    raise RuntimeError(
        "无法解析影像尺寸: " + filepath
    )



# =====================================================
# ISCE XML 元数据模板（与示例 filt_fine.unw.xml 一致）
# =====================================================

ISCE_XML_TEMPLATE = """<imageFile>
    <property name="ISCE_VERSION">
        <value>Release: 2.6.3, svn-, 20230418. Current: svn-.</value>
    </property>
    <property name="access_mode">
        <value>read</value>
        <doc>Image access mode.</doc>
    </property>
    <property name="byte_order">
        <value>l</value>
        <doc>Endianness of the image.</doc>
    </property>
    <component name="coordinate1">
        <factorymodule>isceobj.Image</factorymodule>
        <factoryname>createCoordinate</factoryname>
        <doc>First coordinate of a 2D image (width).</doc>
        <property name="delta">
            <value>1</value>
            <doc>Coordinate quantization.</doc>
        </property>
        <property name="endingvalue">
            <value>{width}</value>
            <doc>Ending value of the coordinate.</doc>
        </property>
        <property name="family">
            <value>imagecoordinate</value>
            <doc>Instance family name</doc>
        </property>
        <property name="name">
            <value>imagecoordinate_name</value>
            <doc>Instance name</doc>
        </property>
        <property name="size">
            <value>{width}</value>
            <doc>Coordinate size.</doc>
        </property>
        <property name="startingvalue">
            <value>0</value>
            <doc>Starting value of the coordinate.</doc>
        </property>
    </component>
    <component name="coordinate2">
        <factorymodule>isceobj.Image</factorymodule>
        <factoryname>createCoordinate</factoryname>
        <doc>Second coordinate of a 2D image (length).</doc>
        <property name="delta">
            <value>1</value>
            <doc>Coordinate quantization.</doc>
        </property>
        <property name="endingvalue">
            <value>{length}</value>
            <doc>Ending value of the coordinate.</doc>
        </property>
        <property name="family">
            <value>imagecoordinate</value>
            <doc>Instance family name</doc>
        </property>
        <property name="name">
            <value>imagecoordinate_name</value>
            <doc>Instance name</doc>
        </property>
        <property name="size">
            <value>{length}</value>
            <doc>Coordinate size.</doc>
        </property>
        <property name="startingvalue">
            <value>0</value>
            <doc>Starting value of the coordinate.</doc>
        </property>
    </component>
    <property name="data_type">
        <value>{data_type}</value>
        <doc>Image data type.</doc>
    </property>
    <property name="extra_file_name">
        <value>{vrt}</value>
        <doc>For example name of vrt metadata.</doc>
    </property>
    <property name="family">
        <value>image</value>
        <doc>Instance family name</doc>
    </property>
    <property name="file_name">
        <value>{data_file}</value>
        <doc>Name of the image file.</doc>
    </property>
    <property name="length">
        <value>{length}</value>
        <doc>Image length</doc>
    </property>
    <property name="name">
        <value>image_name</value>
        <doc>Instance name</doc>
    </property>
    <property name="number_bands">
        <value>{number_bands}</value>
        <doc>Number of image bands.</doc>
    </property>
    <property name="scheme">
        <value>{scheme}</value>
        <doc>Interleaving scheme of the image.</doc>
    </property>
    <property name="width">
        <value>{width}</value>
        <doc>Image width</doc>
    </property>
    <property name="xmax">
        <value>{width}</value>
        <doc>Maximum range value</doc>
    </property>
    <property name="xmin">
        <value>0</value>
        <doc>Minimum range value</doc>
    </property>
</imageFile>
"""



def write_isce_xml(
        outfile,
        data_file,
        vrt_file,
        data_type,
        number_bands,
        scheme,
        width,
        length
):

    content = ISCE_XML_TEMPLATE.format(

        width=width,
        length=length,
        data_type=data_type,
        vrt=os.path.abspath(vrt_file),
        data_file=os.path.abspath(data_file),
        number_bands=number_bands,
        scheme=scheme

    )


    with open(outfile, "w") as f:

        f.write(content)



# =====================================================
# VRT（与示例 filt_fine.unw.vrt / .conncomp.vrt 一致）
# =====================================================

def write_unw_vrt(path, basename, width, length):

    line_offset = width * 8
    band2_offset = width * 4

    content = (
        '<VRTDataset rasterXSize="{w}" rasterYSize="{h}">\n'
        '    <VRTRasterBand dataType="Float32" band="1" subClass="VRTRawRasterBand">\n'
        '        <SourceFilename relativeToVRT="1">{b}</SourceFilename>\n'
        '        <ByteOrder>LSB</ByteOrder>\n'
        '        <ImageOffset>0</ImageOffset>\n'
        '        <PixelOffset>4</PixelOffset>\n'
        '        <LineOffset>{lo}</LineOffset>\n'
        '    </VRTRasterBand>\n'
        '    <VRTRasterBand dataType="Float32" band="2" subClass="VRTRawRasterBand">\n'
        '        <SourceFilename relativeToVRT="1">{b}</SourceFilename>\n'
        '        <ByteOrder>LSB</ByteOrder>\n'
        '        <ImageOffset>{o2}</ImageOffset>\n'
        '        <PixelOffset>4</PixelOffset>\n'
        '        <LineOffset>{lo}</LineOffset>\n'
        '    </VRTRasterBand>\n'
        '</VRTDataset>\n'
    ).format(
        w=width,
        h=length,
        b=basename,
        lo=line_offset,
        o2=band2_offset
    )


    with open(path, "w") as f:

        f.write(content)



def write_conncomp_vrt(path, basename, width, length):

    content = (
        '<VRTDataset rasterXSize="{w}" rasterYSize="{h}">\n'
        '    <VRTRasterBand dataType="Byte" band="1" subClass="VRTRawRasterBand">\n'
        '        <SourceFilename relativeToVRT="1">{b}</SourceFilename>\n'
        '        <ByteOrder>LSB</ByteOrder>\n'
        '        <ImageOffset>0</ImageOffset>\n'
        '        <PixelOffset>1</PixelOffset>\n'
        '        <LineOffset>{lo}</LineOffset>\n'
        '    </VRTRasterBand>\n'
        '</VRTDataset>\n'
    ).format(
        w=width,
        h=length,
        b=basename,
        lo=width
    )


    with open(path, "w") as f:

        f.write(content)



# =====================================================
# conncomp
# =====================================================

def create_conncomp(filename, length, width):

    conn = np.ones(
        (length, width),
        dtype=np.uint8
    )

    conn.tofile(filename)



# =====================================================
# unw 收尾：保证为 2-band (amplitude + phase) BIL
# =====================================================

def finalize_unw(unw_file, int_file, length, width):

    size = os.path.getsize(unw_file)

    expected_1 = length * width * 4    # 单 band 相位 float32
    expected_2 = length * width * 8    # 2 band float32 (ISCE unw)


    if size == expected_2:

        return   # whirlwind 已直接输出 2-band


    if size == expected_1:

        phase = np.fromfile(
            unw_file,
            dtype=np.float32
        ).reshape(length, width)


        cpx = np.fromfile(
            int_file,
            dtype=np.complex64
        ).reshape(length, width)


        unw = np.zeros(
            (length * 2, width),
            dtype=np.float32
        )


        # ISCE unw: band0 amplitude, band1 phase
        unw[0::2, :] = np.abs(cpx)
        unw[1::2, :] = phase


        tmp2 = unw_file + ".tmp2"

        unw.tofile(tmp2)

        os.replace(tmp2, unw_file)

        return


    raise RuntimeError(
        "unw 文件大小异常: %d 字节 (期望 %d 或 %d)"
        % (size, expected_1, expected_2)
    )



# =====================================================
# 单个干涉对
# =====================================================

def unwrap_one(task):


    ifg, nlooks, whirlwind, max_ncomps = task


    name = os.path.basename(ifg)


    int_file = os.path.join(
        ifg,
        "filt_fine.int"
    )

    cor_file = os.path.join(
        ifg,
        "filt_fine.cor"
    )

    unw_file = os.path.join(
        ifg,
        "filt_fine.unw"
    )

    conn_file = os.path.join(
        ifg,
        "filt_fine.unw.conncomp"
    )


    try:


        t = time.time()


        length, width = get_size(int_file)


        # -------------------------------------
        # 解缠
        # -------------------------------------

        if not os.path.exists(unw_file):


            cmd = [

                whirlwind,

                "--ifg",
                int_file,

                "--cor",
                cor_file,

                "--cols",
                str(width),

                "--nlooks",
                str(nlooks),

                "--max-ncomps",
                str(max_ncomps),

                "--out",
                unw_file

            ]


            proc = subprocess.run(

                cmd,

                stdout=subprocess.DEVNULL,

                stderr=subprocess.PIPE,

                text=True

            )


            if proc.returncode != 0:

                diag = []

                for f in (int_file, cor_file):

                    if os.path.exists(f):

                        diag.append(
                            "%s 存在(%d字节)"
                            % (os.path.basename(f), os.path.getsize(f))
                        )

                    else:

                        diag.append(
                            "%s 缺失!" % os.path.basename(f)
                        )


                raise RuntimeError(
                    "whirlwind 退出码 %d: %s | 输入: %s"
                    % (
                        proc.returncode,
                        proc.stderr.strip() or "(无 stderr 输出)",
                        " ; ".join(diag)
                    )
                )


            finalize_unw(
                unw_file,
                int_file,
                length,
                width
            )


            status = "解缠完成"


        else:

            status = "已有unw"


        # -------------------------------------
        # conncomp
        # -------------------------------------

        if not os.path.exists(conn_file):


            create_conncomp(
                conn_file,
                length,
                width
            )


        # -------------------------------------
        # ISCE XML 元数据
        # -------------------------------------

        if not os.path.exists(
            unw_file + ".xml"
        ):

            write_isce_xml(
                unw_file + ".xml",
                unw_file,
                unw_file + ".vrt",
                "FLOAT",
                2,
                "BIL",
                width,
                length
            )


        if not os.path.exists(
            conn_file + ".xml"
        ):

            write_isce_xml(
                conn_file + ".xml",
                conn_file,
                conn_file + ".vrt",
                "BYTE",
                1,
                "BIL",
                width,
                length
            )


        # -------------------------------------
        # VRT
        # -------------------------------------

        if not os.path.exists(
            unw_file + ".vrt"
        ):

            write_unw_vrt(
                unw_file + ".vrt",
                os.path.basename(unw_file),
                width,
                length
            )


        if not os.path.exists(
            conn_file + ".vrt"
        ):

            write_conncomp_vrt(
                conn_file + ".vrt",
                os.path.basename(conn_file),
                width,
                length
            )


        return (

            "\n" + name,

            status,

            time.time() - t

        )


    except Exception as e:


        return (

            name,

            "失败:" + str(e),

            0

        )





# =====================================================
# 主程序
# =====================================================

if __name__ == "__main__":


    args = parse_args()


    ifgs = []


    for d in sorted(
        os.listdir(args.root)
    ):

        p = os.path.join(
            args.root,
            d
        )


        if (
            os.path.isdir(p)

            and "_"
            in d
        ):

            ifgs.append(
                (
                    p,
                    args.nlooks,
                    args.whirlwind,
                    args.max_ncomps
                )
            )


    print("=" * 60)

    print(
        "Whirlwind 并行相位解缠"
    )

    print(
        "干涉对数量:",
        len(ifgs)
    )

    print(
        "并行数:",
        args.process
    )

    print(
        "nlooks:",
        args.nlooks
    )

    print(
        "max-ncomps:",
        args.max_ncomps
    )

    print(
        "whirlwind:",
        args.whirlwind
    )

    print("=" * 60)


    t0 = time.time()


    with Pool(
        processes=args.process
    ) as pool:

        for r in tqdm(

            pool.imap_unordered(

                unwrap_one,

                ifgs

            ),

            total=len(ifgs),

            desc="相位解缠"

        ):

            name, status, cost = r

            tqdm.write(
                f"{name}: {status}"
                + (f" ({cost:.1f}s)" if cost > 0 else "")
            )


    print()

    print("=" * 60)

    print(
        "所有干涉图解缠完成，解缠耗时： %.2f 分钟"
        %
        ((time.time() - t0) / 60)
    )

    print("=" * 60)
