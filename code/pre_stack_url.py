import numpy as np
url_split = np.loadtxt('Sentinel-1下载链接.txt',delimiter='/',dtype=str)
##写入stack.txt
with open('stack.txt','w') as f:
    for i in range(len(url_split)):
        f.write(url_split[i,3]+","+url_split[i,4]+","+url_split[i,6][:-5]+'\n')