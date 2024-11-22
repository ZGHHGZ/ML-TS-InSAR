# ML-TS-InSAR
Mini Long Time-Series InSAR

小区域、长时序InSAR处理包

本软件包基于ISCE、MintPy、MiaplPy开发。

修改了ISCE的stack处理流程，能够仅下载单个burst进行处理，且实现了先裁剪后处理的功能。

SBAS-InSAR：单burst数据下载→裁剪→配准→干涉→解缠→时序处理  （ISCE+MintPy）

DS-InSAR：单burst数据下载→裁剪→配准→相位优化→干涉→解缠→时序处理 （ISCE+MiaplPy+MintPy）

相较于以往处理上百景影像需要数TB的存储空间且耗费大量时间而言，本软件包仅对影像的小区域子集进行处理，仅需数十GB存储空间，避免了非研究区域数据下载、配准、干涉、解缠等处理所造成的无意义空间、时间损耗。


#####################################################################################

操作流程待完善。。。。。。


#####################################################################################
<center class="half">
<img src="https://github.com/ZGHHGZ/Single-Burst-Processing-Flow/blob/main/as.svg" width=400/>
</center>

![image](https://github.com/ZGHHGZ/Single-Burst-Processing-Flow/blob/main/as.svg)
![image](https://github.com/user-attachments/assets/0e5edf1d-3a4d-4669-8b6e-d78ab2c695a3)
