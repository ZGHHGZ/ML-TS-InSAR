import os
import isce 
import shutil
from topsStack import stackSentinel

path= "./run_files"
if  os.path.exists(path):
   shutil.rmtree(path)
path= "./configs"
if  os.path.exists(path):
   shutil.rmtree(path)
os.system('python /home/jovyan/iscepredata/code/topsStack/stackSentinel.py -s ../iscepredata/ -d ./dem/full_res.dem.wgs84  -a ./aux_cal/ -o ./orbits -C NESD  -W slc')

run=[]
run.append("start")
run.append('run_01_unpack_topo_reference')
run.append('run_02_unpack_secondary_slc')
run.append('run_03_average_baseline')
run.append('run_04_extract_burst_overlaps')
run.append('run_05_overlap_geo2rdr')
run.append('run_06_overlap_resample')
run.append('run_07_pairs_misreg')
run.append('run_08_timeseries_misreg')
run.append('run_09_fullBurst_geo2rdr')
run.append('run_10_fullBurst_resample')
run.append('run_11_extract_stack_valid_region')
run.append('run_12_merge_reference_secondary_slc')
run.append('run_13_grid_baseline')
              
#要运行的步骤
runstep=run[7]
              
with open('./run_files/'+str(runstep),"r") as f:
   a=f.readlines()
   print(str(a[0])[0:6])
if str(a[0])[0:6]!="python":
    for i in range(len(a)):
       a[i]="python /home/jovyan/iscepredata/code/topsStack/"+a[i]
    with open('./run_files/'+str(runstep),"w") as f:
        for i in range(len(a)):
            f.write(str(a[i]))
os.system('python /home/jovyan/iscepredata/code/topsStack/run.py -i ./run_files/'+str(runstep)+' -p 8')