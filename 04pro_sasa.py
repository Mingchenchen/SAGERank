###                                             by scc in 2023/10/10
########################################################################################################################
#######                                                                                                     ############
#######                                       produce sasa for inteface                                     ############
#######                                                                                                     ############
########################################################################################################################
import os
a1='/data4_large1/home_data/ccsun/scc_neuralnetwork/02/pro-pro/zxc/name.txt'
with open(a1,'r') as rw:
    file =rw.readlines()
    #file=file.strip()
    for i in file:
        i=i.rstrip('\n')
        print(i)
        a='/data4_large1/home_data/ccsun/scc_neuralnetwork/02/pro-pro/zxc/pdb2/%s/inter/rec'%i
        path_list=os.listdir(a)
        path_list.sort()
        for ii in path_list:
            os.system('freesasa /data4_large1/home_data/ccsun/scc_neuralnetwork/02/pro-pro/zxc/pdb2/%s/inter/rec/%s -n 100 --depth=residue --format=seq -o /data4_large1/home_data/ccsun/scc_neuralnetwork/02/pro-pro/zxc/pdb2/%s/sasa/rec/%s'%(i,ii,i,ii))