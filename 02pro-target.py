###                                   by scc in 2023/10/9
########################################################################################################################
#######                                                                                                     ############
#######                                produce target                                                     ############
#######                                                                                                     ############
########################################################################################################################

import os
import time
from pdb2sql import StructureSimilarity

filepath='/data2/data_home/ccsun/link/scc_neuralnetwork/Recognition_of_antigens_by_antibodies/try1/train/1ADQ/structures'   #改！！！
os.mkdir('/data2/data_home/ccsun/link/scc_neuralnetwork/Recognition_of_antigens_by_antibodies/try1/train/1ADQ/txt')       #改！！！
num2=0  
num=0
begin_time = time.time()

for i in range(1000):
    path_list=os.listdir(filepath)
    path_list.sort()
    pdb=path_list[i]
    pdb2=path_list[i]
    pdb2=pdb2.strip('.pdb')
    #print(pdb2)
    aa=os.path.join(filepath,pdb)      
    pdb=aa

    sim=StructureSimilarity(aa,'/data2/data_home/ccsun/link/deeprank/data1/1ADQ/native/1ADQ.pdb') #改！！！
    irmsd_fast=sim.compute_irmsd_fast(method='svd',izone=None)
    print(irmsd_fast)
    cc=pdb2+'    '+str(irmsd_fast)
    print(cc)
    irmsd_file='/data2/data_home/ccsun/link/scc_neuralnetwork/Recognition_of_antigens_by_antibodies/try1/train/1ADQ/irmsd_file_1ADQ.txt'
    with open(irmsd_file,'a') as ab:
        ab.write(cc+'\n')
        
    
    print('%s的irmsd为:%f'%(pdb2,irmsd_fast))
    b=i+1
    print('processing....⌚..,Wait a moment😀,total 1000,and %d have been processed'%b)
    
    
    if irmsd_fast<4:
        os.chdir('/data2/data_home/ccsun/link/scc_neuralnetwork/Recognition_of_antigens_by_antibodies/try1/train/1ADQ/txt')   #改！！！
        txt_file_name='%s'%pdb2+'.txt'
        with open(txt_file_name,'a') as rw:
            rw.write('1')
        num+=1
     
    else:
        os.chdir('/data2/data_home/ccsun/link/scc_neuralnetwork/Recognition_of_antigens_by_antibodies/try1/train/1ADQ/txt')  #改！！！
        txt_file_name='%s'%pdb2+'.txt'
        with open(txt_file_name,'a') as rw:
            rw.write('0')
        num2+=1    

end_time = time.time()
run_time = round(end_time-begin_time)
hour = run_time//3600
minute = (run_time-3600*hour)//60
second = run_time-3600*hour-60*minute
print('程序运行时长为：%d分钟'%minute)
print('在1000个decoys中，irmsd小于4埃的有%d个'%num)
print('在1000个decoys中，irmsd小于4埃的有%d个'%num2)