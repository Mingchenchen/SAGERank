###                                   by scc in 2023/10/9
########################################################################################################################
#######                                                                                                     ############
#######                               produce grb-nor npz                                  ############
#######                                                                                                     ############
########################################################################################################################

import os
import numpy as np
np.set_printoptions(suppress=True) 
#np.set_printoptions(threshold = np.inf)
import math
import operator
from functools import reduce
from scipy.spatial import distance_matrix
ee=math.e

def Prepare_Input(pssm_rec_path,pssm_lig_path,sasa_rec_path,sasa_lig_path,rec_path,lig_path,txt_path,savepath):
    for num in range(1):
        
        interface_rec_list=os.listdir(rec_path)
        interface_rec_list.sort()
        interface_pdb=interface_rec_list[num]
        pdbid=interface_pdb.replace('.pdb','')
        rec_pdb=os.path.join(rec_path,interface_pdb)
        
        interface_lig_list=os.listdir(lig_path)
        interface_lig_list.sort()
        interface_pdb1=interface_lig_list[num]
        lig_pdb=os.path.join(lig_path,interface_pdb1)
        
        path_list=os.listdir(pssm_rec_path)
        path_list.sort()
        pdb=path_list[num]
        aa=os.path.join(pssm_rec_path,pdb)
        
        path_list1=os.listdir(pssm_lig_path)
        path_list1.sort()
        pdb1=path_list1[num]
        aaa=os.path.join(pssm_lig_path,pdb1)
        
        sasa_rec_list=os.listdir(sasa_rec_path)
        sasa_rec_list.sort()
        sasa_rec_list1=sasa_rec_list[num]
        sasa_rec=os.path.join(sasa_rec_path,sasa_rec_list1)
        
        sasa_lig_list=os.listdir(sasa_lig_path)
        sasa_lig_list.sort()
        sasa_lig_list1=sasa_lig_list[num]
        sasa_lig=os.path.join(sasa_lig_path,sasa_lig_list1)
        
        # 归一化1
        def dic_normalize(dic): #计算方法为：该值-最小值/最大值-最小值
            #print(dic)
            max_value = dic[max(dic, key=dic.get)]
            min_value = dic[min(dic, key=dic.get)]
            #print(max_value)
            #print(min_value)
            interval = float(max_value) - float(min_value)
            #print('interval',interval)
            for key in dic.keys():
                dic[key] = round((dic[key] - min_value) / interval,4)
            dic['X'] = (max_value + min_value) / 2.0
            return dic
        
        # 归一化2
        def normalize_pssm(value):
            pssm=1/(1+ee**-value)
            pssm=round(pssm,4)
            return pssm
        
        
        pro_res_table = ['ALA ', 'CYS ', 'ASP ', 'GLU ', 'PHE ', 'GLY ', 'HIS ', 'ILE ', 'LYS ', 'LEU ', 'MET ', 'ASN ', 'PRO ', 'GLN ', 
                             'ARG ', 'SER ', 
                             'THR ', 'VAL ', 'TRP ', 'TYR ']
        
        #Whether the residue is aliphatic(脂肪族氨基酸)
        pro_res_aliphatic_table = ['ALA ', 'ILE ', 'LEU ', 'MET ', 'VAL ']
        
        #Whether the residue is aromatic（芳香性氨基酸）
        pro_res_aromatic_table = ['PHE ', 'TRP ', 'TYR ']
        
        #Whether the residue is polar neutral（中性氨基酸）
        pro_res_polar_neutral_table = ['CYS ', 'ASN ', 'GLN ', 'SER ', 'THR ']
        
        #Whether the residue is acidic charged（酸性氨基酸）
        pro_res_acidic_charged_table = ['ASP ', 'GLU ']
        
        #Whether the residue is basic charged（碱性氨基酸）
        pro_res_basic_charged_table = ['HIS ', 'LYS ', 'ARG ']
        
        #Residue weight（分子量）
        res_weight_table = {'ALA ': 71.08, 'CYS ': 103.15, 'ASP ': 115.09, 'GLU ': 129.12, 'PHE ': 147.18, 'GLY ': 57.05, 'HIS ': 137.14,
                                        'ILE ': 113.16, 'LYS ': 128.18, 'LEU ': 113.16, 'MET ': 131.20, 'ASN ': 114.11, 'PRO ': 97.12, 'GLN ': 128.13,
                                        'ARG ': 156.19, 'SER ': 87.08, 'THR ': 101.11, 'VAL ': 99.13, 'TRP ': 186.22, 'TYR ': 163.18}
        
        #The negative of the logarithm of the dissociation constant for the –COOH group（羧基的解离常数的负对数pK1）
        res_pka_table = {'ALA ': 2.34, 'CYS ': 1.96, 'ASP ': 1.88, 'GLU ': 2.19, 'PHE ': 1.83, 'GLY ': 2.34, 'HIS ': 1.82, 'ILE ': 2.36,
                                     'LYS ': 2.18, 'LEU ': 2.36, 'MET ': 2.28, 'ASN ': 2.02, 'PRO ': 1.99, 'GLN ': 2.17, 'ARG ': 2.17, 'SER ': 2.21,
                                     'THR ': 2.09, 'VAL ': 2.32, 'TRP ': 2.83, 'TYR ': 2.32}
        
        #The negative of the logarithm of the dissociation constant for the –NH3 group（氨基的解离常数的负对数pK2）
        res_pkb_table = {'ALA ': 9.69, 'CYS ': 10.28, 'ASP ': 9.60, 'GLU ': 9.67, 'PHE ': 9.13, 'GLY ': 9.60, 'HIS ': 9.17,
                 'ILE ': 9.60, 'LYS ': 8.95, 'LEU ': 9.60, 'MET ': 9.21, 'ASN ': 8.80, 'PRO ': 10.60, 'GLN ': 9.13,
                 'ARG ': 9.04, 'SER ': 9.15, 'THR ': 9.10, 'VAL ': 9.62, 'TRP ': 9.39, 'TYR ': 9.62}

        #The pH at the isoelectric point（等电点）
        res_pl_table = {'ALA ': 6.00, 'CYS ': 5.07, 'ASP ': 2.77, 'GLU ': 3.22, 'PHE ': 5.48, 'GLY ': 5.97, 'HIS ': 7.59,
                        'ILE ': 6.02, 'LYS ': 9.74, 'LEU ': 5.98, 'MET ': 5.74, 'ASN ': 5.41, 'PRO ': 6.30, 'GLN ': 5.65,
                        'ARG ': 10.76, 'SER ': 5.68, 'THR ': 5.60, 'VAL ': 5.96, 'TRP ': 5.89, 'TYR ': 5.96}

        #Hydrophobicity of residue (pH=2)
        res_hydrophobic_ph2_table = {'ALA ': 47, 'CYS ': 52, 'ASP ': -18, 'GLU ': 8, 'PHE ': 92, 'GLY ': 0, 'HIS ': -42, 'ILE ': 100,
                                     'LYS ': -37, 'LEU ': 100, 'MET ': 74, 'ASN ': -41, 'PRO ': -46, 'GLN ': -18, 'ARG ': -26, 'SER ': -7,
                                     'THR ': 13, 'VAL ': 79, 'TRP ': 84, 'TYR ': 49}

        #Hydrophobicity of residue (pH=7)
        res_hydrophobic_ph7_table = {'ALA ': 41, 'CYS ': 49, 'ASP ': -55, 'GLU ': -31, 'PHE ': 100, 'GLY ': 0, 'HIS ': 8, 'ILE ': 99,
                                     'LYS ': -23, 'LEU ': 97, 'MET ': 74, 'ASN ': -28, 'PRO ': -46, 'GLN ': -10, 'ARG ': -14, 'SER ': -5,
                                     'THR ': 13, 'VAL ': 76, 'TRP ': 97, 'TYR ': 63}
        
        res_weight_table = dic_normalize(res_weight_table)
        res_pka_table = dic_normalize(res_pka_table)
        res_pkb_table =dic_normalize(res_pkb_table)
        res_pl_table =dic_normalize(res_pl_table)
        res_hydrophobic_ph2_table =dic_normalize(res_hydrophobic_ph2_table)
        res_hydrophobic_ph7_table = dic_normalize(res_hydrophobic_ph7_table)
        
        def one_of_k_encoding(x, allowable_set):
            if x not in allowable_set:
                # print(x)
                raise Exception('input {0} not in allowable set{1}:'.format(x, allowable_set))
            return list(map(lambda s: x == s, allowable_set))

        
        def residue_features(residue,pssm1,pssm2,pssm3,pssm4,pssm5,pssm6,pssm7,pssm8,pssm9,pssm10,pssm11,pssm12,pssm13,pssm14,pssm15,pssm16,pssm17,pssm18,pssm19,pssm20):
            res_property1 = [1 if residue in pro_res_aliphatic_table else 0, 
                                1 if residue in pro_res_aromatic_table else 0,
                                1 if residue in pro_res_polar_neutral_table else 0,
                                1 if residue in pro_res_acidic_charged_table else 0,
                                1 if residue in pro_res_basic_charged_table else 0]
            
            res_property2 = [res_weight_table[residue], res_pka_table[residue], res_pkb_table[residue], res_pl_table[residue]]
            
            res_property3 =[normalize_pssm(pssm1),normalize_pssm(pssm2),normalize_pssm(pssm3),normalize_pssm(pssm4),normalize_pssm(pssm5)
                           ,normalize_pssm(pssm6),normalize_pssm(pssm7),normalize_pssm(pssm8),normalize_pssm(pssm9),normalize_pssm(pssm10)
                           ,normalize_pssm(pssm11),normalize_pssm(pssm12),normalize_pssm(pssm13),normalize_pssm(pssm14),normalize_pssm(pssm15)
                           ,normalize_pssm(pssm16),normalize_pssm(pssm17),normalize_pssm(pssm18),normalize_pssm(pssm19),normalize_pssm(pssm20)]
            #print(res_property1+res_property2+res_property3)
            #print(np.array(res_property1 + res_property2+res_property3).shape)
            #return np.array(res_property1 + res_property2+res_property3)
            features=(res_property1 + res_property2+res_property3)
            return features
        
        
        res_list=[]
        res_list1=[]
        value1=[]
        with open(sasa_rec,'r') as rec_file:
            line = rec_file.readlines() 
            for x in line[1:-1]: 
                x=x.strip()
                resid=x[12:16]
                sasa_value=x[19:26]
                value1.append(sasa_value)
                res_list.append(resid)
        maxnum=float(max(value1))
        minnum=float(min(value1))
        internal_value=maxnum-minnum
        value2=[]
        for x1 in value1:
            x1=float(x1)
            final_value=round((x1-minnum)/internal_value,4)
            #print(final_value)
            value2.append(final_value)
        
        value3=[]
        with open(sasa_lig,'r') as lig_file:
            line = lig_file.readlines() 
            for y in line[1:-1]:
                y=y.strip()
                resid=y[12:16]
                sasa_value=y[19:26]
                #name.append(resid)
                value3.append(sasa_value)
                res_list1.append(resid)
        maxnum1=float(max(value3))
        minnum1=float(min(value3))
        internal_value1=maxnum1-minnum1
        value4=[]
        for y1 in value3:
            y1=float(y1)
            final_value=round((y1-minnum)/internal_value1,4)
            #print(final_value)
            value4.append(final_value)        
        

        rec_fea=[]
        lig_fea=[]
        with open(aa,'r') as rec_file:
            line = rec_file.readlines() 
            for m,i in enumerate(line):
                i=i.strip()
                resid=i[0:4]
                pssm1=int(i[5:8])
                pssm2=int(i[9:12])
                pssm3=int(i[13:16])
                pssm4=int(i[17:20])
                pssm5=int(i[21:24])
                pssm6=int(i[25:28])
                pssm7=int(i[29:32])
                pssm8=int(i[33:36])
                pssm9=int(i[37:40])
                pssm10=int(i[41:44])
                pssm11=int(i[45:48])
                pssm12=int(i[49:52])
                pssm13=int(i[53:56])
                pssm14=int(i[57:60])
                pssm15=int(i[61:64])
                pssm16=int(i[65:68])
                pssm17=int(i[69:72])
                pssm18=int(i[73:76])
                pssm19=int(i[77:80])
                pssm20=int(i[81:84])
                
                sasa= [value2[m]]
                
                pro_seq=[res_list[m]]
                pro_hot = np.zeros((1, len(pro_res_table)))
                for i in range(len(pro_seq)):
                    pro_hot[i,] = one_of_k_encoding(pro_seq[i], pro_res_table)
                pro_hot1=reduce(operator.add, pro_hot)  
                pro_hot1=pro_hot1.tolist() 
                
                rec_feature=residue_features(resid,pssm1,pssm2,pssm3,pssm4,pssm5,pssm6,pssm7,pssm8,pssm9,pssm10,pssm11,pssm12,pssm13,pssm14,pssm15,pssm16,pssm17,pssm18,pssm19,pssm20)
                
                rec_fea.append(pro_hot1+rec_feature+sasa)
                
            rec_fea=np.array(rec_fea) 
            
        with open(aaa,'r') as lig_file:
            line = lig_file.readlines() 
            for z,n in enumerate(line):
                n=n.strip()
                resid=n[0:4]
                pssm1=int(n[5:8])
                pssm2=int(n[9:12])
                pssm3=int(n[13:16])
                pssm4=int(n[17:20])
                pssm5=int(n[21:24])
                pssm6=int(n[25:28])
                pssm7=int(n[29:32])
                pssm8=int(n[33:36])
                pssm9=int(n[37:40])
                pssm10=int(n[41:44])
                pssm11=int(n[45:48])
                pssm12=int(n[49:52])
                pssm13=int(n[53:56])
                pssm14=int(n[57:60])
                pssm15=int(n[61:64])
                pssm16=int(n[65:68])
                pssm17=int(n[69:72])
                pssm18=int(n[73:76])
                pssm19=int(n[77:80])
                pssm20=int(n[81:84])
                
                sasa= [value4[z]]
                
                pro_seq=[res_list1[z]]
                pro_hot = np.zeros((1, len(pro_res_table)))
                for i in range(len(pro_seq)):
                    pro_hot[i,] = one_of_k_encoding(pro_seq[i], pro_res_table)
                pro_hot1=reduce(operator.add, pro_hot)  
                pro_hot1=pro_hot1.tolist()
                
                lig_feature=residue_features(resid,pssm1,pssm2,pssm3,pssm4,pssm5,pssm6,pssm7,pssm8,pssm9,pssm10,pssm11,pssm12,pssm13,pssm14,pssm15,pssm16,pssm17,pssm18,pssm19,pssm20)
                
                lig_fea.append(pro_hot1+lig_feature+sasa)
                
            lig_fea=np.array(lig_fea)
            
        H = np.concatenate([rec_fea, lig_fea], 0)
        
        ca_accord_lig=[]
        ca_accord_rec=[]
        with open(lig_pdb,'r') as file:
            file=file.readlines()
            for i in file:
                x=[float(i[31:38])]
                y=[float(i[39:46])]
                z=[float(i[47:54])]
                atom=i[13:16]
                #print(atom)
                if atom=='CA ':
                    ca_accord_lig.append(x+y+z)
            ca_accord_lig=np.array(ca_accord_lig)
        ligand_count=len(ca_accord_lig)
        with open(rec_pdb,'r') as file:
            file=file.readlines()
            for i in file:
                x=[float(i[31:38])]
                y=[float(i[39:46])]
                z=[float(i[47:54])]
                atom=i[13:16]
                #print(atom)
                if atom=='CA ':
                    ca_accord_rec.append(x+y+z)
            ca_accord_rec=np.array(ca_accord_rec)
        receptor_count=len(ca_accord_rec)
        #print(ca_accord_rec)
        #print(ca_accord_lig)
        dis_rec_rec = distance_matrix(ca_accord_rec, ca_accord_rec)
        dis_lig_lig = distance_matrix(ca_accord_lig, ca_accord_lig)
        dis_rec_lig =distance_matrix(ca_accord_rec,ca_accord_lig)
        dis_lig_rec =distance_matrix(ca_accord_lig,ca_accord_rec)
            
        for ii in np.nditer(dis_rec_rec,op_flags=['readwrite']):
            if ii >4.5:
                ii[...]=0
            else:
                ii[...]=1

        for ii in np.nditer(dis_lig_lig,op_flags=['readwrite']):
            if ii >4.5:
                ii[...]=0
            else:
                ii[...]=1

        for ii in np.nditer(dis_rec_lig,op_flags=['readwrite']):
            if ii >10:
                ii[...]=0
            else:
                ii[...]=1

        for ii in np.nditer(dis_lig_rec,op_flags=['readwrite']):
            if ii >10:
                ii[...]=0
            else:
                ii[...]=1
            
        adj= np.zeros((receptor_count + ligand_count, receptor_count + ligand_count))
        adj[:receptor_count, :receptor_count] = dis_rec_rec
        adj[receptor_count:, receptor_count:] = dis_lig_lig
        adj[:receptor_count, receptor_count:] = dis_rec_lig
        adj[receptor_count:, :receptor_count] = dis_lig_rec
        

        target=[]
        txt_list=os.listdir(txt_path)
        txt_list.sort()
        txt=txt_list[num]
        txt_name=os.path.join(txt_path,txt)
        with open(txt_name,'r') as tmp_file:
            line=tmp_file.readline()
            line=line.strip()  
            aim_tmp=int(line)
            target.append(aim_tmp)
        
        
        input_file=os.path.join(savepath,"%s.npz"%pdbid)
        np.savez(input_file,  H=H, A1=adj, T=target)