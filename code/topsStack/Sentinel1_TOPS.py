#!/usr/bin/env python3 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Copyright 2014 California Institute of Technology. ALL RIGHTS RESERVED.
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
# United States Government Sponsorship acknowledged. This software is subject to
# U.S. export control laws and regulations and has been classified as 'EAR99 NLR'
# (No [Export] License Required except when exporting to an embargoed country,
# end user, or in support of a prohibited end use). By downloading this software,
# the user agrees to comply with all applicable U.S. export laws and regulations.
# The user has the responsibility to obtain export licenses, or other export
# authority as may be required before exporting this software to any 'EAR99'
# embargoed foreign country or citizen of those countries.
#
# Author: Piyush Agram
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



import isce
import datetime
import isceobj
from isceobj.Planet.Planet import Planet
from isceobj.Orbit.Orbit import StateVector, Orbit
from isceobj.Orbit.OrbitExtender import OrbitExtender
from isceobj.Planet.AstronomicalHandbook import Const
from iscesys.Component.Component import Component
from iscesys.DateTimeUtil.DateTimeUtil import DateTimeUtil as DTUtil
import os
import glob
import numpy as np


def _install_all_valid_line_compat(Sentinel1):
    """兼容裁剪后所有方位行均有效的标准 Sentinel-1 annotation。

    ISCE2 原始 populateBurstSpecificMetadata 仅在 firstValidSample 从非负值
    再变回 -1 时设置 numValidLines。精确裁剪可能恰好只保留有效区，此时数组
    全为非负值，ISCE2 会让 numValidLines 保持 None 并在后续加法中崩溃。

    这里仅在解析期间给内存中的数组临时追加一个 -1 结束标记，使原算法得到
    numValidLines == 实际行数；finally 中立即恢复 XML 文本，因此不会修改 SAFE，
    也不会损失裁剪影像的最后一行。
    """
    if getattr(Sentinel1, '_all_valid_line_compat_installed', False):
        return

    original = Sentinel1.populateBurstSpecificMetadata

    def compatible_populate(self):
        touched = []
        burst_list = self.getxmlelement('swathTiming/burstList')
        for burst in burst_list:
            first = burst.find('firstValidSample')
            last = burst.find('lastValidSample')
            if first is None or last is None:
                continue
            first_values = (first.text or '').split()
            last_values = (last.text or '').split()
            if (first_values and len(first_values) == len(last_values)
                    and int(first_values[-1]) >= 0):
                touched.append((first, first.text, last, last.text))
                first.text = (first.text or '').rstrip() + ' -1'
                last.text = (last.text or '').rstrip() + ' -1'
        try:
            return original(self)
        finally:
            for first, first_text, last, last_text in touched:
                first.text = first_text
                last.text = last_text

    Sentinel1.populateBurstSpecificMetadata = compatible_populate
    Sentinel1._all_valid_line_compat_installed = True

def createParser():
    import argparse

    parser = argparse.ArgumentParser( description = 'Sentinel parser' )

    parser.add_argument('-d', '--dirname', dest='dirname', type=str,
            default=None, help='SAFE format directory. (Recommended)')

    parser.add_argument('-o', '--outdir', dest='outdir', type=str,
            required=True, help='Output SLC prefix.')

    parser.add_argument('-p', '--orbit', dest='orbit', type=str,
            default=None, help='Precise orbit file, Use of orbitdir preferred')

    parser.add_argument('-a', '--aux', dest='auxprod', type=str,
            default=None, help='Auxiliary product with antenna gains, Use of auxdir preferred')

    parser.add_argument('--orbitdir', dest='orbitdir', type=str,
            default=None, help = 'Directory with all the orbits')

    parser.add_argument('--auxdir', dest='auxdir', type=str,
            default=None, help = 'Directory with all the aux products')

    parser.add_argument('--pol', dest='polid', type=str,
            default='vv', help = 'Polarization of interest. Default: vv')

    parser.add_argument('-b', '--bbox', dest='bbox', type=str,
            default=None, help='Lat/Lon Bounding SNWE')

    parser.add_argument('-s', '--swaths', dest='swaths', type=str,
            default=None, help='list pf swaths')
    return parser


def main(iargs=None):
    from isceobj.Sensor.TOPS.Sentinel1 import Sentinel1

    _install_all_valid_line_compat(Sentinel1)

    
    inps = cmdLineParse(iargs)

    if inps.swaths is None:
       inps.swaths = [1,2,3]
    else:
       inps.swaths = [int(i) for i in inps.swaths.split()]

    for swath in inps.swaths:

       obj = Sentinel1()
       obj.configure()
       obj.safe = inps.dirname.split()
       obj.swathNumber = swath
       obj.output = os.path.join(inps.outdir, 'IW{0}'.format(swath))
       obj.orbitFile = inps.orbit
       obj.auxFile = inps.auxprod
       obj.orbitDir = inps.orbitdir
       obj.auxDir = inps.auxdir
       obj.polarization = inps.polid
       if inps.bbox is not None:
          obj.regionOfInterest = [float(x) for x in inps.bbox.split()]
       try:
          obj.parse()
          obj.extractImage(virtual=True)
       except Exception as e:
          message = str(e)
          # burst2safe/裁切 SAFE 可以只包含一个实际 IW。配置仍列出 1 2 3 时，
          # 对缺失 IW 正常跳过；实际存在的 IW 发生任何错误则必须向上传播，
          # 避免 run_00 看似成功、run_01 最后才以空 boxes 难以定位地失败。
          if ('No annotation xml file found' in message or
                  'No tiff file found' in message):
             print('Skipping unavailable IW{0}: {1}'.format(swath, message))
             continue
          raise

def cmdLineParse(iargs=None):
    '''
    Command Line Parser.
    '''
    parser = createParser()
    inps = parser.parse_args(args=iargs)
    return inps

if __name__ == '__main__':
     
    main()


