# LM-TS-InSAR
Long Mini Time-Series InSAR

长时序、小区域时序InSAR处理包

本软件包基于ISCE、Mintpy、MiaplPy开发。

修改了ISCE的stack处理流程，能够仅下载单个burst进行处理，且实现了先裁剪后处理的功能。

单burst数据下载→裁剪→配准→干涉→解缠→时序处理  （ISCE+MintPy）

单burst数据下载→裁剪→配准→相位优化→干涉→解缠→时序处理 （ISCE+MiaplPy+MintPy）

相较于以往处理上百景影像需要数TB的存储空间且耗费大量时间，本软件包对上百景影像的小区域子集进行处理，仅需数十GB存储空间，避免了非研究区域配准、干涉、解缠等处理所造成的无意义时间损耗。
