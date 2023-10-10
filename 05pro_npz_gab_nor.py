###                                   by scc in 2023/10/9
########################################################################################################################
#######                                                                                                     ############
#######                               produce gab-nor npz                                   ############
#######                                                                                                     ############
########################################################################################################################


import os
import numpy as np
import os
from rdkit.Chem.rdmolfiles import MolFromPDBFile
from data_processing.Feature_Processing import get_atom_feature
import numpy as np
#np.set_printoptions(threshold = np.inf)
from rdkit.Chem.rdmolops import GetAdjacencyMatrix
from scipy.spatial import distance_matrix


def Prepare_Input(receptor_path,ligand_path,txtpath,savepath):
    # extract the interface region
    #root_path=os.path.split(structure_path)[0]  #
    interface_rec_list=os.listdir(receptor_path)
    interface_rec_list.sort()
    for i,ii in enumerate(interface_rec_list):
        pdbid=ii.replace('.pdb','')
        rec_pdb=os.path.join(receptor_path,ii)
       
        lig_pdb=os.path.join(ligand_path,ii)
        
        receptor_mol = MolFromPDBFile(rec_pdb, sanitize=False)
        ligand_mol = MolFromPDBFile(lig_pdb, sanitize=False)   
    
    
        receptor_count = receptor_mol.GetNumAtoms()
        ligand_count = ligand_mol.GetNumAtoms()  
    
        receptor_feature = get_atom_feature(receptor_mol, is_ligand=False) 
        ligand_feature = get_atom_feature(ligand_mol, is_ligand=True)
   
        # get receptor adj matrix
        c1 = receptor_mol.GetConformers()[0]
        d1 = np.array(c1.GetPositions()) 
    
    
        adj1 = GetAdjacencyMatrix(receptor_mol) + np.eye(receptor_count)
   
        #np.shape(GetAdjacencyMatrix(receptor_mol)) 为（387*387），np.shape(np.eye(receptor_count)) 为（387*387）
    
 
        # get ligand adj matrix
        c2 = ligand_mol.GetConformers()[0]
        d2 = np.array(c2.GetPositions())
        adj2 = GetAdjacencyMatrix(ligand_mol) + np.eye(ligand_count)
    
        # combine analys                                                                                                                                                                    is
        H = np.concatenate([receptor_feature, ligand_feature], 0)
   
    
        agg_adj1 = np.zeros((receptor_count + ligand_count, receptor_count + ligand_count))
        agg_adj1[:receptor_count, :receptor_count] = adj1
        agg_adj1[receptor_count:, receptor_count:] = adj2  # array without r-l interaction
        
        dm = distance_matrix(d1, d2)
        
        dm_shape=np.shape(dm)
    
        for ii in np.nditer(dm,op_flags=['readwrite']):

        #print(i)
            if ii >8:
                ii[...]=0
       
            else:
                ii[...]=1
    
    
        agg_adj2 = np.copy(agg_adj1)
        agg_adj2[:receptor_count, receptor_count:] = np.copy(dm)
        agg_adj2[receptor_count:, :receptor_count] = np.copy(np.transpose(dm))  # with interaction array
        # node indice for aggregation
        #print(agg_adj2)
    
        valid = np.zeros((receptor_count + ligand_count,))  
        valid[:receptor_count] = 1
    
 
        #aim_file_list = [x for x in os.listdir(txtpath) if 'txt' in x]
        #aim_file_list.sort()
        
        
        
        target=[]
        path_list1=os.listdir(txtpath)
        path_list1.sort()
        txt=path_list1[i]
        bb=os.path.join(txtpath,txt)
        with open(bb,'r') as tmp_file:
            line=tmp_file.readline()
            line=line.strip()  
            aim_tmp=int(line)
            target.append(aim_tmp)
        
        
        #aim_path = os.path.join(save_path, 'aimset.npy')
        #np.save(aim_path, atom40_output)
    
        input_file=os.path.join(savepath,"%s.npz"%pdbid)
    
        np.savez(input_file,  H=H, A1=agg_adj2, T=target)
    #return input_file
