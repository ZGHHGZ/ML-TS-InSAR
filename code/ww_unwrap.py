#!/usr/bin/env python3

import os
import time
import argparse
import subprocess
import shutil
import numpy as np

from multiprocessing import Pool
from tqdm import tqdm
from xml.etree import ElementTree as ET



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
        "-r",
        "--root",
        type=str,
        default="./merged/interferograms",
        help="干涉图目录"
    )


    return parser.parse_args()





# =====================================================
# 读取 ISCE XML 尺寸
# =====================================================

def read_size(xml):

    tree = ET.parse(xml)

    root = tree.getroot()

    width=None
    length=None


    for p in root.iter("property"):

        name=p.attrib.get("name")

        v=p.find("value")

        if v is None:
            continue


        if name=="width":
            width=int(v.text)

        elif name=="length":
            length=int(v.text)



    if width is None or length is None:

        raise RuntimeError(
            "无法读取尺寸 "+xml
        )


    return length,width





# =====================================================
# 生成 XML
# =====================================================

def make_xml(

        template,

        outfile,

        filename,

        datatype,

        image_type,

        bands,

        scheme

):


    tree=ET.parse(template)

    root=tree.getroot()


    for p in root.iter("property"):

        name=p.attrib.get("name")

        v=p.find("value")


        if v is None:
            continue


        if name=="file_name":

            v.text=os.path.abspath(filename)


        elif name=="data_type":

            v.text=datatype


        elif name=="image_type":

            v.text=image_type


        elif name=="number_bands":

            v.text=str(bands)


        elif name=="scheme":

            v.text=scheme



    tree.write(

        outfile,

        encoding="UTF-8",

        xml_declaration=True

    )






# =====================================================
# 生成unw
# =====================================================

def convert_whirlwind_output(

        phase_file,

        int_file,

        outfile,

        length,

        width

):


    phase=np.fromfile(

        phase_file,

        dtype=np.float32

    ).reshape(

        length,

        width

    )



    cpx=np.fromfile(

        int_file,

        dtype=np.complex64

    ).reshape(

        length,

        width

    )



    unw=np.zeros(

        (

            length*2,

            width

        ),

        dtype=np.float32

    )


    # ISCE unw:
    # band0 amplitude
    # band1 phase

    unw[0::2,:]=np.abs(cpx)

    unw[1::2,:]=phase



    unw.tofile(outfile)





# =====================================================
# conncomp
# =====================================================

def create_conncomp(

        filename,

        length,

        width

):


    conn=np.ones(

        (

            length,

            width

        ),

        dtype=np.uint8

    )


    conn.tofile(filename)






# =====================================================
# 单个干涉对
# =====================================================

def unwrap_one(task):


    ifg,nlooks=task


    name=os.path.basename(ifg)


    int_file=os.path.join(
        ifg,
        "filt_fine.int"
    )


    int_xml=int_file+".xml"


    cor_file=os.path.join(
        ifg,
        "fine.cor"
    )


    unw_file=os.path.join(
        ifg,
        "filt_fine.unw"
    )


    conn_file=os.path.join(
        ifg,
        "filt_fine.unw.conncomp"
    )


    tmp=os.path.join(
        ifg,
        "whirlwind.tmp"
    )



    try:


        t=time.time()


        if not os.path.exists(int_xml):

            return name,"缺少int xml",0



        length,width=read_size(int_xml)



        # -------------------------------------
        # 解缠
        # -------------------------------------

        if not os.path.exists(unw_file):


            cmd=[

                "whirlwind",

                "--ifg",
                int_file,

                "--cor",
                cor_file,

                "--nlooks",
                str(nlooks),

                "--out",
                tmp

            ]


            subprocess.run(

                cmd,

                stdout=subprocess.DEVNULL,

                stderr=subprocess.DEVNULL,

                check=True

            )



            convert_whirlwind_output(

                tmp,

                int_file,

                unw_file,

                length,

                width

            )


            if os.path.exists(tmp):

                os.remove(tmp)


            status="解缠完成"


        else:

            status="已有unw"



        # -------------------------------------
        # unw xml
        # -------------------------------------

        if not os.path.exists(
            unw_file+".xml"
        ):


            make_xml(

                int_xml,

                unw_file+".xml",

                unw_file,

                "float",

                "unw",

                2,

                "BIL"

            )



        # -------------------------------------
        # conncomp
        # -------------------------------------

        if not os.path.exists(conn_file):


            create_conncomp(

                conn_file,

                length,

                width

            )



        if not os.path.exists(
            conn_file+".xml"
        ):


            make_xml(

                int_xml,

                conn_file+".xml",

                conn_file,

                "byte",

                "conncomp",

                1,

                "BIP"

            )



        # -------------------------------------
        # VRT
        # -------------------------------------

        for f in [
            unw_file,
            conn_file
        ]:


            if not os.path.exists(
                f+".vrt"
            ):


                subprocess.run(

                    [

                        "gdal_translate",

                        "-of",
                        "VRT",

                        f,

                        f+".vrt"

                    ],

                    stdout=subprocess.DEVNULL,

                    stderr=subprocess.DEVNULL

                )



        return (

            "\n" + name,

            status,

            time.time()-t

        )



    except Exception as e:


        return (

            name,

            "失败:"+str(e),

            0

        )







# =====================================================
# 主程序
# =====================================================

if __name__=="__main__":


    args=parse_args()



    ifgs=[]


    for d in sorted(
        os.listdir(args.root)
    ):


        p=os.path.join(
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
                    args.nlooks
                )

            )



    print("="*60)

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

    print("="*60)




    t0=time.time()



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


            name,status,cost=r

            tqdm.write(
                f"{name}: {status}"
                + (f" ({cost:.1f}s)" if cost > 0 else "")
            )



    print()

    print("="*60)

    print(
        "所有干涉图解缠完成，解缠耗时： %.2f 分钟"
        %
        ((time.time()-t0)/60)
    )

    print("="*60)