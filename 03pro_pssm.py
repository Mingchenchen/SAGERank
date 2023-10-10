###                                        by scc in 2032/10/10
########################################################################################################################
#######                                                                                                     ############
#######                                        produce pssm                                                 ############
#######                                                                                                     ############
########################################################################################################################
import os 
from pssmgen import PSSM

path='/data4_large1/home_data/ccsun/scc_neuralnetwork/02/pro-pro/zxc/name.txt'

with open(path,'r') as rw:
    file =rw.readlines()
    #file=file.strip()
    for i in file[:1]:
        i=i.rstrip('\n')
        print(i)
        #i=i.strip('.pdb')
        gen = PSSM(work_dir='/data4_large1/home_data/ccsun/scc_neuralnetwork/02/pro-pro/zxc/pdb2/%s'%i)

        gen.configure(blast_exe='/data/mab_lab/blast/ncbi-blast-2.11.0+/bin/psiblast',
                database='/data/mab_lab/blast/ncbi-blast-2.11.0+/db/nr',
                num_threads = 8, evalue=0.0001, comp_based_stats='T',
                max_target_seqs=2000, num_iterations=3, outfmt=7,
                save_each_pssm=True, save_pssm_after_last_round=True)

        gen.get_pssm(fasta_dir='fasta', out_dir='pssm_raw', run=True, save_all_psiblast_output=True)

        gen.map_pssm(pssm_dir='pssm_raw', pdb_dir='./preprocess/inter/rec', out_dir='./preprocess/rec_pssm', chain=('C'))
