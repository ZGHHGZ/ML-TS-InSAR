# ML-TS-InSAR
Mini Long Time-Series InSAR

小区域、长时序InSAR处理包

本软件包基于ISCE、MintPy、MiaplPy开发。

修改了ISCE的stack处理流程，能够仅下载单个burst进行处理，且实现了先裁剪后处理的功能。

SBAS-InSAR：单burst数据下载→裁剪→配准→干涉→解缠→时序处理  （ISCE+MintPy）

DS-InSAR：单burst数据下载→裁剪→配准→相位优化→干涉→解缠→时序处理 （ISCE+MiaplPy+MintPy）

相较于以往处理上百景影像需要数TB的存储空间且耗费大量时间而言，本软件包仅对影像的小区域子集进行处理，仅需数十GB存储空间，避免了非研究区域数据下载、配准、干涉、解缠等处理所造成的无意义空间、时间损耗。


#####################################################################################

使用说明待完善。。。。。。


#####################################################################################
<div style="display: flex; justify-content: center; gap: 20px;">
    <img src="https://raw.githubusercontent.com/ZGHHGZ/Single-Burst-Processing-Flow/main/as.svg" width="400" alt="升轨图像" />
    <img src="https://raw.githubusercontent.com/ZGHHGZ/Single-Burst-Processing-Flow/main/des.svg" width="400" alt="降轨图像" />
</div>

![image](https://github.com/user-attachments/assets/0e5edf1d-3a4d-4669-8b6e-d78ab2c695a3)








#####################################################################################
#########################               教程1               #########################
#####################################################################################
# Sentinel-1小区域长时序InSAR本地高效处理案例

## 1 SARGet数据搜索及下载（Sentinel-1）
SARGet软件参考链接：  
<https://github.com/ZGHHGZ/SARGet.exe>  

以某研究区域的中心经纬度（经度：118.0971；纬度：36.7292）为例进行搜索：  
- 搜索结果：共3组数据覆盖该位置；  
- 实验选取：分组编号为069_146274_IW2的时序影像；  
- 时间范围：2022年1月1日至2026年3月5日；  
- 筛选结果：共95景影像。  

<center>图1 筛选数据结果</center>  

将95景影像下载至自定义空文件夹，本案例下载路径为：  
`C:\InSAR_process\S1_process`  

## 2 时序InSAR预处理
### 2.1 处理环境配置
本流程适用于小区域长时序InSAR高效处理，对计算机性能要求较低，可在个人电脑运行，环境配置步骤如下：  

| 步骤 | 操作说明 | 参考链接 |
|------|----------|----------|
| 1 | 安装WSL（Windows Subsystem for Linux） | <https://learn.microsoft.com/zh-cn/windows/wsl/install> |
| 2 | 在WSL下安装conda | - |
| 3 | 通过hyp3-isce2安装ISCE2软件 | <https://github.com/ASFHyP3/hyp3-isce2> |

> 注意：安装ISCE2后，需手动进入conda环境下的ISCE2包路径，将`contrib`文件夹复制到对应位置（`contrib`文件夹来源：<https://github.com/isce-framework/isce2>）。  

<center>图2 contrib文件夹</center>  

### 2.2 数据处理
ISCE2处理Sentinel-1数据的最小单位为Burst，默认流程会产生大量冗余数据；ML-TS-InSAR支持单Burst自定义裁剪，仅处理研究区范围内数据，大幅提升效率、降低存储占用。  

#### 2.2.1 工具准备
ML-TS-InSAR链接：<https://github.com/ZGHHGZ/ML-TS-InSAR>  

文件放置规则：  
1. 将ML-TS-InSAR的`code`文件夹置于SARGet下载文件夹的同级目录：  
   - SARGet下载路径（案例）：`C:\InSAR_process\S1_process`  
   - `code`文件夹路径（案例）：`C:\InSAR_process\`  
2. 将ML-TS-InSAR内的`pre_parameter.txt`文件放入SARGet下载文件夹内，并按需求设定处理参数。  

<center>图3 处理参数设定</center>  

#### 2.2.2 命令行执行
1. 打开WSL命令行窗口，进入SARGet下载路径的WSL映射目录：  
   案例Windows路径：`C:\InSAR_process\S1_process`  
   对应WSL路径：`/mnt/c/InSAR_process/S1_process`  

2. 执行时序预处理命令：  
   ```bash
   python ../code/pre_data_SARget.py
