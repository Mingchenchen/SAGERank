###                                   by scc in 2023/10/9
########################################################################################################################
#######                                                                                                     ############
#######                               produce grb-cdr npz                                   ############
#######                                                                                                     ############
########################################################################################################################

import os
import numpy as np
np.set_printoptions(suppress=True)  # 取消科学计数法输出
#np.set_printoptions(threshold = np.inf)
import math
import operator
from functools import reduce
from scipy.spatial import distance_matrix
ee=math.e

import requests

class annotate():
    """
    class `annotate`. 
    
    Initiator `__init__` has 2 parameters:
    
    :param aaseq: STRING: A single-letter, amino acid sequence corresponding to the complete VH or VL chain. Both uppercase and lowercase are accepted. 
    
    :param scheme: STRING: "kabat", "chothia", "contact", or "imgt". Must be in lowercase
    
    Class has 3 methods. `retrieve()`: retrieves numbered seqs from Abnum website, then sends it to method `analyze` to determine the FR and CDR regions, and to `output() ` to print the result and return a list of 2 dictionaries, the first of which contains to region:seq pairs, the second of which contains number:residue pairs. 
    
    """
    
    def __init__(self, aaseq, scheme):
        
        self.aaseq=aaseq
        self.scheme=scheme
    
    def __repr__(self):
        return "Annotation of VH or VL sequence using Kabat, Chothia, Contact, or IMGT scheme"
    
    def output(self, chain, lst, regionlst):
        """
        Prints the FR and CDR regions and their corresponding seq. It returns a `list` of 2 `dict`. 
        
        :param chain: STRING, either "H" or "L" in uppercase
        :param lst:  LIST, a list of residue and their corresponding numbers in kabat or chothia scheme
        :param regionlst: LIST, a list of peptides, each corresponds to a FR or CDR region
        :return: LIST, a list of 2 `dict`, The first dict consists of region: seq pairs. The second dict consists of number:residue pairs.
        
        """
        
        self.chain=chain
        self.lst=lst
        self.regionlst=regionlst

        self.regiondict, self.numberdict={}, {}
        print('------------',self.chain)
        for i in range (0, len(self.lst), 2):
            self.numberdict[self.lst[i]]=self.lst[i+1]
        
        
        if self.scheme=="kabat":
            print("Annotation in Kabat scheme:")
        elif self.scheme=="chothia":
            print("Annotation in Chothia scheme:")
        elif self.scheme=="contact":
            print("Annotation in Contact scheme:")
        else:
            print("Annotation in IMGT scheme:")
        
        if self.chain=="L":
            print("L-FR1:  ", self.regionlst[0])
            print("L-CDR1: ", self.regionlst[1])
            print("L-FR2:  ", self.regionlst[2])
            print("L-CDR2: ", self.regionlst[3])
            print("L-FR3:  ", self.regionlst[4])
            print("L-CDR3: ", self.regionlst[5])
            print("L-FR4:  ", self.regionlst[6])
            
            for region, seq in zip(["L-FR1", "L-CDR1", "L-FR2","L-CDR2", "L-FR3", "L-CDR3", "L-FR4"], self.regionlst):
                self.regiondict[region]=seq
            
            return [self.regiondict, self.numberdict]
                
        else:
            print("H-FR1:  ", self.regionlst[0])
            print("H-CDR1: ", self.regionlst[1])
            print("H-FR2:  ", self.regionlst[2])
            print("H-CDR2: ", self.regionlst[3])
            print("H-FR3:  ", self.regionlst[4])
            print("H-CDR3: ", self.regionlst[5])
            print("H-FR4:  ", self.regionlst[6])
            
            for region, seq in zip(["H-FR1", "H-CDR1", "H-FR2","H-CDR2", "H-FR3", "H-CDR3", "H-FR4"], self.regionlst):
                self.regiondict[region]=seq
            
            return [self.regiondict, self.numberdict]
            
        
        
        

    
    def analyze(self,chain, lst):
        """
        Define CDR and FR regions based on the numbered sequence returned from website
        
        :param chain: STRING, "H" or "L" in uppercase
        :param lst: LIST, a list of residue and their corresponding numbers in kabat or chothia scheme
        :return: LIST, a list of strings, where each string is a peptide corresponding to the a region, in the order of: FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4
        
        :raises: `ValueError` if any of the FR or CDR region is missing
        
        """
        
        self.chain=chain
        self.lst=lst
        if self.chain=="L":
            self.L_FR1, self.L_CDR1, self.L_FR2, self.L_CDR2, self.L_FR3, self.L_CDR3, self.L_FR4=["" for i in range (0, 7)]
            
            try:
                if self.scheme in ["kabat", "chothia"]:
                    self.L_FR1="".join([self.lst[i+1] for i in range (0, self.lst.index("L24"), 2)])
                    self.L_CDR1="".join([self.lst[i+1] for i in range (self.lst.index("L24"), self.lst.index("L35"), 2)])
                    self.L_FR2="".join([self.lst[i+1] for i in range (self.lst.index("L35"), self.lst.index("L50"), 2)])
                    self.L_CDR2="".join([self.lst[i+1] for i in range (self.lst.index("L50"), self.lst.index("L57"), 2)])
                    self.L_FR3="".join([self.lst[i+1] for i in range (self.lst.index("L57"), self.lst.index("L89"), 2)])
                    self.L_CDR3="".join([self.lst[i+1] for i in range (self.lst.index("L89"), self.lst.index("L98"), 2)])
                    self.L_FR4="".join([self.lst[i+1] for i in range (self.lst.index("L98"), len(self.lst), 2)])
                                    
                elif self.scheme =="contact": 
                    self.L_FR1="".join([self.lst[i+1] for i in range (0, self.lst.index("L30"), 2)])
                    self.L_CDR1="".join([self.lst[i+1] for i in range (self.lst.index("L30"), self.lst.index("L37"), 2)])
                    self.L_FR2="".join([self.lst[i+1] for i in range (self.lst.index("L37"), self.lst.index("L46"), 2)])
                    self.L_CDR2="".join([self.lst[i+1] for i in range (self.lst.index("L46"), self.lst.index("L56"), 2)])
                    self.L_FR3="".join([self.lst[i+1] for i in range (self.lst.index("L56"), self.lst.index("L89"), 2)])
                    self.L_CDR3="".join([self.lst[i+1] for i in range (self.lst.index("L89"), self.lst.index("L97"), 2)])
                    self.L_FR4="".join([self.lst[i+1] for i in range (self.lst.index("L97"), len(self.lst), 2)])
                                    
                else: #IMGT scheme
                    self.L_FR1="".join([self.lst[i+1] for i in range (0, self.lst.index("L27"), 2)])
                    self.L_CDR1="".join([self.lst[i+1] for i in range (self.lst.index("L27"), self.lst.index("L33"), 2)])
                    self.L_FR2="".join([self.lst[i+1] for i in range (self.lst.index("L33"), self.lst.index("L50"), 2)])
                    self.L_CDR2="".join([self.lst[i+1] for i in range (self.lst.index("L50"), self.lst.index("L52"), 2)])
                    self.L_FR3="".join([self.lst[i+1] for i in range (self.lst.index("L52"), self.lst.index("L89"), 2)])
                    self.L_CDR3="".join([self.lst[i+1] for i in range (self.lst.index("L89"), self.lst.index("L98"), 2)])
                    self.L_FR4="".join([self.lst[i+1] for i in range (self.lst.index("L98"), len(self.lst), 2)])
                
                return [self.L_FR1, self.L_CDR1, self.L_FR2, self.L_CDR2, self.L_FR3, self.L_CDR3, self.L_FR4] 

            except ValueError:
                print("Unable to retrieve complete V region. Make sure the sequence has complete V region")
            except:
                print("An error occured")
        else:
            self.H_FR1, self.H_CDR1, self.H_FR2, self.H_CDR2, self.H_FR3, self.H_CDR3, self.H_FR4=["" for i in range (0, 7)]
            try:
                if self.scheme=="kabat":
                    self.H_FR1="".join([self.lst[i+1] for i in range (0, self.lst.index("H31"), 2)])
                    self.H_CDR1="".join([self.lst[i+1] for i in range (self.lst.index("H31"), self.lst.index("H36"), 2)])
                    self.H_FR2="".join([self.lst[i+1] for i in range (self.lst.index("H36"), self.lst.index("H50"), 2)])
                    self.H_CDR2="".join([self.lst[i+1] for i in range (self.lst.index("H50"), self.lst.index("H66"), 2)])
                    self.H_FR3="".join([self.lst[i+1] for i in range (self.lst.index("H66"), self.lst.index("H95"), 2)])
                    self.H_CDR3="".join([self.lst[i+1] for i in range (self.lst.index("H95"), self.lst.index("H103"), 2)])
                    self.H_FR4="".join([self.lst[i+1] for i in range (self.lst.index("H103"), len(self.lst), 2)])            
            
                elif self.scheme=="chothia":
                    self.H_FR1="".join([self.lst[i+1] for i in range (0, self.lst.index("H26"), 2)])
                    self.H_CDR1="".join([self.lst[i+1] for i in range (self.lst.index("H26"), self.lst.index("H33"), 2)])
                    self.H_FR2="".join([self.lst[i+1] for i in range (self.lst.index("H33"), self.lst.index("H52"), 2)])
                    self.H_CDR2="".join([self.lst[i+1] for i in range (self.lst.index("H52"), self.lst.index("H57"), 2)])
                    self.H_FR3="".join([self.lst[i+1] for i in range (self.lst.index("H57"), self.lst.index("H95"), 2)])
                    self.H_CDR3="".join([self.lst[i+1] for i in range (self.lst.index("H95"), self.lst.index("H103"), 2)])
                    self.H_FR4="".join([self.lst[i+1] for i in range (self.lst.index("H103"), len(self.lst), 2)])

                elif self.scheme=="contact":
                    self.H_FR1="".join([self.lst[i+1] for i in range (0, self.lst.index("H30"), 2)])
                    self.H_CDR1="".join([self.lst[i+1] for i in range (self.lst.index("H30"), self.lst.index("H36"), 2)])
                    self.H_FR2="".join([self.lst[i+1] for i in range (self.lst.index("H36"), self.lst.index("H47"), 2)])
                    self.H_CDR2="".join([self.lst[i+1] for i in range (self.lst.index("H47"), self.lst.index("H59"), 2)])
                    self.H_FR3="".join([self.lst[i+1] for i in range (self.lst.index("H59"), self.lst.index("H93"), 2)])
                    self.H_CDR3="".join([self.lst[i+1] for i in range (self.lst.index("H93"), self.lst.index("H102"), 2)])
                    self.H_FR4="".join([self.lst[i+1] for i in range (self.lst.index("H102"), len(self.lst), 2)])
                                        
                else: #IMGT scheme
                    self.H_FR1="".join([self.lst[i+1] for i in range (0, self.lst.index("H26"), 2)])
                    self.H_CDR1="".join([self.lst[i+1] for i in range (self.lst.index("H26"), self.lst.index("H34"), 2)])
                    self.H_FR2="".join([self.lst[i+1] for i in range (self.lst.index("H34"), self.lst.index("H51"), 2)])
                    self.H_CDR2="".join([self.lst[i+1] for i in range (self.lst.index("H51"), self.lst.index("H58"), 2)]) #51>57 (instead of 56)
                    self.H_FR3="".join([self.lst[i+1] for i in range (self.lst.index("H58"), self.lst.index("H93"), 2)])
                    self.H_CDR3="".join([self.lst[i+1] for i in range (self.lst.index("H93"), self.lst.index("H103"), 2)])
                    self.H_FR4="".join([self.lst[i+1] for i in range (self.lst.index("H103"), len(self.lst), 2)])                    
                
                return [self.H_FR1, self.H_CDR1, self.H_FR2, self.H_CDR2, self.H_FR3, self.H_CDR3, self.H_FR4]                    

            except ValueError:
                print("Unable to retrieve complete V region. Make sure the sequence has complete V region")
            except:
                print("An error occured in the `analyze()` method")
        
    def retrieve (self):
        """
        Retrieve numbered residues from Abnum website
        
        :return: returns same object from the `output()` method. 
        
        :raises: `ValueError` if input scheme is not among "kabat", "chothia", "contact", and "imgt"
        
        """
        
        self.url="http://www.bioinf.org.uk/abs/abnum/abnum.cgi"
        
        try: 
            if self.scheme not in ["kabat", "chothia", "contact", "imgt"]:
                raise Exception
            
        except ValueError:
            print("Incorrect scheme mode. Must be one of the following (lowercase): kabat, chothia, contact, imgt")
        
        else:
            if self.scheme=="kabat":
                self.sche="-k"
            else:
                self.sche="-c"
        
        try:
            self.d={"plain":1, "scheme":self.sche, "aaseq":self.aaseq}
            self.myPage=requests.get(self.url, params=self.d)
            self.text=self.myPage.text
            self.lst=self.text.split()
                
            if len(self.lst)>1:
                self.chain=self.lst[0][0]
                self.result=self.output(self.chain, self.lst, self.analyze(self.chain, self.lst))
                return self.result
            else:
                print("No annotation retrieved. Did you enter the complete VH or VL sequence?")
        except:
            print("An error occured in the `retrieve()` method")



def Prepare_Input(pssm_rec_path,pssm_lig_path,sasa_rec_path,sasa_lig_path,rec_path,lig_path,txt_path,savepath):
   
    interface_rec_list=os.listdir(rec_path)
    interface_rec_list.sort()
    
    for num,numm in enumerate(interface_rec_list):
        #num=num+360000
        nu=num
        interface_pdb=numm
        pdbid=interface_pdb.replace('.pdb','')
        #print(numm)
        rec_pdb=os.path.join(rec_path,interface_pdb)
        interface_pdb1=numm
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

        def dic_normalize(dic): 
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
        
        def normalize_pssm(value):
            pssm=1/(1+ee**-value)
            pssm=round(pssm,4)
            return pssm
        
        
        pro_res_table = ['ALA ', 'CYS ', 'CSD ','CSX ', 'ASP ', 'GLU ', 'CGU ', 'PHE ', 'GLY ', 'HIS ', 'ILE ', 'LYS ', 'MLY ','LEU ', 'MET ', 'MSE ', 'ASN ', 'PRO ', 'GLN ', 
                             'ARG ', 'SER ', 
                             'THR ', 'VAL ', 'TRP ', 'TYR ', 'ACE ', 'PCA ']
        
        #Whether the residue is aliphatic(脂肪族氨基酸)
        pro_res_aliphatic_table = ['ALA ', 'ILE ', 'LEU ', 'MET ', 'MSE', 'VAL ']
        
        #Whether the residue is aromatic（芳香性氨基酸）
        pro_res_aromatic_table = ['PHE ', 'TRP ', 'TYR ']
        
        #Whether the residue is polar neutral（中性氨基酸）
        pro_res_polar_neutral_table = ['CYS ', 'CSD ','CSX ', 'ASN ', 'GLN ',  'PCA ', 'SER ', 'THR ']
        
        #Whether the residue is acidic charged（酸性氨基酸）
        pro_res_acidic_charged_table = ['ASP ', 'GLU ', 'ACE ', 'CGU ']
        
        #Whether the residue is basic charged（碱性氨基酸）
        pro_res_basic_charged_table = ['HIS ', 'LYS ', 'MLY ', 'ARG ']
        
        #Residue weight（分子量）
        res_weight_table = {'ALA ': 71.08, 'CYS ': 103.15, 'CSX ': 103.15, 'CSD ': 103.15, 'ASP ': 115.09, 'ACE ': 115.09, 'GLU ': 129.12, 'CGU ': 129.12, 'PHE ': 147.18, 'GLY ': 57.05, 'HIS ': 137.14,
                                        'ILE ': 113.16, 'LYS ': 128.18, 'MLY ': 128.18, 'LEU ': 113.16, 'MET ': 131.20, 'MSE ': 131.20, 'ASN ': 114.11, 'PRO ': 97.12, 'GLN ': 128.13, 'PCA ': 128.13,
                                        'ARG ': 156.19, 'SER ': 87.08, 'THR ': 101.11, 'VAL ': 99.13, 'TRP ': 186.22, 'TYR ': 163.18}
        
        #The negative of the logarithm of the dissociation constant for the –COOH group（羧基的解离常数的负对数pK1）
        res_pka_table = {'ALA ': 2.34, 'CYS ': 1.96, 'CSX ': 1.96, 'CSD ': 1.96, 'ASP ': 1.88, 'ACE ': 1.88, 'GLU ': 2.19, 'CGU ': 2.19, 'PHE ': 1.83, 'GLY ': 2.34, 'HIS ': 1.82, 'ILE ': 2.36,
                                     'LYS ': 2.18, 'MLY ': 2.18,'LEU ': 2.36, 'MET ': 2.28, 'MSE ': 2.28, 'ASN ': 2.02, 'PRO ': 1.99, 'GLN ': 2.17, 'PCA ': 2.17, 'ARG ': 2.17, 'SER ': 2.21,
                                     'THR ': 2.09, 'VAL ': 2.32, 'TRP ': 2.83, 'TYR ': 2.32}
        
        #The negative of the logarithm of the dissociation constant for the –NH3 group（氨基的解离常数的负对数pK2）
        res_pkb_table = {'ALA ': 9.69, 'CYS ': 10.28, 'CSX ': 10.28, 'CSD ': 10.28, 'ASP ': 9.60, 'ACE ': 9.60, 'GLU ': 9.67, 'CGU ': 9.67, 'PHE ': 9.13, 'GLY ': 9.60, 'HIS ': 9.17,
                 'ILE ': 9.60, 'LYS ': 8.95, 'MLY ': 8.95, 'LEU ': 9.60, 'MET ': 9.21, 'MSE ': 9.21, 'ASN ': 8.80, 'PRO ': 10.60, 'GLN ': 9.13, 'PCA ': 9.13,
                 'ARG ': 9.04, 'SER ': 9.15, 'THR ': 9.10, 'VAL ': 9.62, 'TRP ': 9.39, 'TYR ': 9.62}

        #The pH at the isoelectric point（等电点）
        res_pl_table = {'ALA ': 6.00, 'CYS ': 5.07, 'CSX ': 5.07, 'CSD ': 5.07, 'ASP ': 2.77, 'ACE ': 2.77, 'GLU ': 3.22, 'CGU ': 3.22, 'PHE ': 5.48, 'GLY ': 5.97, 'HIS ': 7.59,
                        'ILE ': 6.02, 'LYS ': 9.74, 'MLY ': 9.74, 'LEU ': 5.98, 'MET ': 5.74, 'MSE ': 5.74, 'ASN ': 5.41, 'PRO ': 6.30, 'GLN ': 5.65, 'PCA ': 5.65,
                        'ARG ': 10.76, 'SER ': 5.68, 'THR ': 5.60, 'VAL ': 5.96, 'TRP ': 5.89, 'TYR ': 5.96}

        #Hydrophobicity of residue (pH=2)
        res_hydrophobic_ph2_table = {'ALA ': 47, 'CYS ': 52, 'CSX ': 52, 'CSD ': 52, 'ASP ': -18, 'ACE ': -18, 'GLU ': 8, 'CGU ': 8, 'PHE ': 92, 'GLY ': 0, 'HIS ': -42, 'ILE ': 100,
                                     'LYS ': -37, 'MLY ': -37, 'LEU ': 100, 'MET ': 74, 'MSE ': 74, 'ASN ': -41, 'PRO ': -46, 'GLN ': -18, 'PCA ': -18, 'ARG ': -26, 'SER ': -7, 'THR ': 13, 'VAL ': 79, 'TRP ': 84, 'TYR ': 49}

        #Hydrophobicity of residue (pH=7)
        res_hydrophobic_ph7_table = {'ALA ': 41, 'CYS ': 49, 'CSX ': 49, 'CSD ': 49, 'ASP ': -55, 'ACE ': -55, 'GLU ': -31, 'CGU ': -31, 'PHE ': 100, 'GLY ': 0, 'HIS ': 8, 'ILE ': 99,
                                     'LYS ': -23, 'MLY ': -23, 'LEU ': 97, 'MET ': 74, 'MSE ': 74, 'ASN ': -28, 'PRO ': -46, 'GLN ': -10, 'PCA ': -10, 'ARG ': -14, 'SER ': -5,
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
        scc=[]
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
                residue_num=i[23:26]
                atom=i[13:16]
                #print(atom)
                if atom=='CA ':
                    ca_accord_rec.append(x+y+z)
                    scc.append(residue_num)
            ca_accord_rec=np.array(ca_accord_rec)
        receptor_count=len(ca_accord_rec)
        #print(scc)
        #print(ca_accord_lig)
        dis_rec_rec = distance_matrix(ca_accord_rec, ca_accord_rec)
        dis_lig_lig = distance_matrix(ca_accord_lig, ca_accord_lig)
        dis_rec_lig =distance_matrix(ca_accord_rec,ca_accord_lig)
        dis_lig_rec =distance_matrix(ca_accord_lig,ca_accord_rec)
        
        
        pdbidd=numm[:4]
        fasta_path='/data4_large1/home_data/ccsun/scc_neuralnetwork/Recognition_of_antigens_by_antibodies/fasta/all-fasta-c/%s.C.fasta'%pdbidd
        print('fasta is '：',fasta_path)
        fasta_data = open(fasta_path,'r')
        fasta_seq = fasta_data.readlines()
        #print(fasta_seq)
        seq_all = ''
        for i in fasta_seq[1:]:
            seq_all += i[:-1]
        
        num = len(seq_all)//2
        #print(num)
        seq_1 = seq_all[:num]
        seq_2 = seq_all[num:]
        #print(1,seq_1,2,seq_2)
        
        scheme="chothia"
        seq1=annotate(seq_1,scheme)
        result1=seq1.retrieve()
        seq2=annotate(seq_2,scheme)
        result2=seq2.retrieve()
        
        def Find_cdr(result):
            list_1 = []
            list_2=[]
            list_3=[]
            all_list = list(result[0].values())
            #print('all_list是：',all_list)
            
            cdr1 = all_list[1]
            #print(Get_Cdr_num(cdr1,seq_all))
            a1=Get_Cdr_num(cdr1,seq_all)
            for i in range(a1[0],a1[1]+1):
                list_1.append(i)
                  
            cdr2 = all_list[3]
            a2=Get_Cdr_num(cdr2,seq_all)
            for i in range(a2[0],a2[1]+1):
                list_2.append(i)
                #print(i)
        
            cdr3 = all_list[5]
            a3=Get_Cdr_num(cdr3,seq_all)
            for i in range(a3[0],a3[1]+1):
                list_3.append(i)
                #print(i)
                
            return list_1,list_2,list_3
        
        def Get_Cdr_num(cdr,seq_all):
            for i in range(len(seq_all)-len(cdr)+1):
            #     print(i)
                cdr_num_last = i+len(cdr)
            #     print(cdr_num_last)
            #     print(seq_all[i:cdr_num_last])
                if seq_all[i:cdr_num_last] == cdr:
                    return i+1,cdr_num_last
        print('-------------')
        
        list_all = Find_cdr(result1) + Find_cdr(result2)
        cdrh1=list_all[0]
        cdrh2=list_all[1]
        cdrh3=list_all[2]
        cdrl1=list_all[3]
        cdrl2=list_all[4]
        cdrl3=list_all[5]
        #print(cdrh1,cdrh2,cdrh3)
        #print(cdrl1,cdrl2,cdrl3)
        
        adj11=np.zeros((receptor_count,receptor_count))  
        cdrh1_num=[]
        for zz in cdrh1:
            for p1,p2 in enumerate(scc):
                p3=int(p2.strip())
                #print(type(p2))
                if p3==zz:
                    cdrh1_num.append(p1)
       # print(cdrh1_num)
        
        cdrh2_num=[]
        for zz in cdrh2:
            for p1,p2 in enumerate(scc):
                p3=int(p2.strip())
                #print(type(p2))
                if p3==zz:
                    cdrh2_num.append(p1)
        #print(cdrh2_num)
        
        cdrh3_num=[]
        for zz in cdrh3:
            for p1,p2 in enumerate(scc):
                p3=int(p2.strip())
                #print(type(p2))
                if p3==zz:
                    cdrh3_num.append(p1)
        #print(cdrh3_num)
        
        cdrl1_num=[]
        for zz in cdrl1:
            for p1,p2 in enumerate(scc):
                p3=int(p2.strip())
                #print(type(p2))
                if p3==zz:
                    cdrl1_num.append(p1)
        #print(cdrl1_num)
        
        cdrl2_num=[]
        for zz in cdrl2:
            for p1,p2 in enumerate(scc):
                p3=int(p2.strip())
                #print(type(p2))
                if p3==zz:
                    cdrl2_num.append(p1)
        #print(cdrl2_num)
        
        cdrl3_num=[]
        for zz in cdrl3:
            for p1,p2 in enumerate(scc):
                p3=int(p2.strip())
                #print(type(p2))
                if p3==zz:
                    cdrl3_num.append(p1)
        #print(cdrl3_num)
          
        #[15, 16, 17, 18, 19, 20, 21, 22, 30, 31, 32, 33, 37, 38, 39, 40] cdrh1
        #[48, 49, 50, 51, 54, 55, 56, 57, 68, 69, 70, 71, 76, 77, 78, 79, 83, 84, 85, 86]  cdrh2
        #[99, 100, 101, 102]  cdrh3
        #[]    cdr11
        #[]    cdr12
        #[]    cdr13
        
        #cdrh1与cdrh2
        if len(cdrh1_num)>0 and len(cdrh2_num)>0:
            for pp in cdrh1_num:
                #print(pp)
                for pp1 in cdrh2_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                
        #cdrh1与cdrh3
        if len(cdrh1_num)>0 and len(cdrh3_num)>0:
            for pp in cdrh1_num:
                #print(pp)
                for pp1 in cdrh3_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1       
        
        #cdrh1与cdrl1
        if len(cdrh1_num)>0 and len(cdrl1_num)>0:
            for pp in cdrh1_num:
                #print(pp)
                for pp1 in cdrl1_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
        
       #cdrh1与cdrl2
        if len(cdrh1_num)>0 and len(cdrl2_num)>0:
            for pp in cdrh1_num:
                #print(pp)
                for pp1 in cdrl2_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #cdrh1与cdrl3
        if len(cdrh1_num)>0 and len(cdrl3_num)>0:
            for pp in cdrh1_num:
                #print(pp)
                for pp1 in cdrl3_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #####################################################
    
        #cdrh2与cdrh1
        if len(cdrh2_num)>0 and len(cdrh1_num)>0:
            for pp in cdrh2_num:
                #print(pp)
                for pp1 in cdrh1_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #cdrh2与cdrh3
        if len(cdrh2_num)>0 and len(cdrh3_num)>0:
            for pp in cdrh2_num:
                #print(pp)
                for pp1 in cdrh3_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #cdrh2与cdrl1
        if len(cdrh2_num)>0 and len(cdrl1_num)>0:
            for pp in cdrh2_num:
                #print(pp)
                for pp1 in cdrl1_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #cdrh2与cdrl2
        if len(cdrh2_num)>0 and len(cdrl2_num)>0:
            for pp in cdrh2_num:
                #print(pp)
                for pp1 in cdrl2_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #cdrh2与cdrl3
        if len(cdrh2_num)>0 and len(cdrl3_num)>0:
            for pp in cdrh2_num:
                #print(pp)
                for pp1 in cdrl3_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        ###########################################################
        
        #cdrh3与cdrh1
        if len(cdrh3_num)>0 and len(cdrh1_num)>0:
            for pp in cdrh3_num:
                #print(pp)
                for pp1 in cdrh1_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #cdrh3与cdrh2
        if len(cdrh3_num)>0 and len(cdrh2_num)>0:
            for pp in cdrh3_num:
                #print(pp)
                for pp1 in cdrh2_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #cdrh3与cdrl1
        if len(cdrh3_num)>0 and len(cdrl1_num)>0:
            for pp in cdrh3_num:
                #print(pp)
                for pp1 in cdrl1_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #cdrh3与cdrl2
        if len(cdrh3_num)>0 and len(cdrl2_num)>0:
            for pp in cdrh3_num:
                #print(pp)
                for pp1 in cdrl2_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #cdrh3与cdrl3
        if len(cdrh3_num)>0 and len(cdrl3_num)>0:
            for pp in cdrh3_num:
                #print(pp)
                for pp1 in cdrl3_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
        
        ##########################################################
                    
        #cdrl1与cdrh1
        if len(cdrl1_num)>0 and len(cdrh1_num)>0:
            for pp in cdrl1_num:
                #print(pp)
                for pp1 in cdrh1_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #cdrl1与cdrh2
        if len(cdrl1_num)>0 and len(cdrh2_num)>0:
            for pp in cdrl1_num:
                #print(pp)
                for pp1 in cdrh2_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #cdrl1与cdrh3
        if len(cdrl1_num)>0 and len(cdrh3_num)>0:
            for pp in cdrl1_num:
                #print(pp)
                for pp1 in cdrh3_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #cdrl1与cdrl2
        if len(cdrl1_num)>0 and len(cdrl2_num)>0:
            for pp in cdrl1_num:
                #print(pp)
                for pp1 in cdrl2_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #cdrl1与cdrl3
        if len(cdrl1_num)>0 and len(cdrl3_num)>0:
            for pp in cdrl1_num:
                #print(pp)
                for pp1 in cdrl3_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
        
        ##########################################################
        
        #cdrl2与cdrh1
        if len(cdrl2_num)>0 and len(cdrh1_num)>0:
            for pp in cdrl2_num:
                #print(pp)
                for pp1 in cdrh1_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #cdrl2与cdrh2
        if len(cdrl2_num)>0 and len(cdrh2_num)>0:
            for pp in cdrl2_num:
                #print(pp)
                for pp1 in cdrh2_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #cdrl2与cdrh3
        if len(cdrl2_num)>0 and len(cdrh3_num)>0:
            for pp in cdrl2_num:
                #print(pp)
                for pp1 in cdrh3_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #cdrl2与cdrl1
        if len(cdrl2_num)>0 and len(cdrl1_num)>0:
            for pp in cdrl2_num:
                #print(pp)
                for pp1 in cdrl1_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #cdrl2与cdrl3
        if len(cdrl2_num)>0 and len(cdrl3_num)>0:
            for pp in cdrl2_num:
                #print(pp)
                for pp1 in cdrl3_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #######################################################
        
        #cdrl3与cdrh1
        if len(cdrl3_num)>0 and len(cdrh1_num)>0:
            for pp in cdrl3_num:
                #print(pp)
                for pp1 in cdrh1_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #cdrl3与cdrh2
        if len(cdrl3_num)>0 and len(cdrh2_num)>0:
            for pp in cdrl3_num:
                #print(pp)
                for pp1 in cdrh2_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #cdrl3与cdrh3
        if len(cdrl3_num)>0 and len(cdrh3_num)>0:
            for pp in cdrl3_num:
                #print(pp)
                for pp1 in cdrh3_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #cdrl3与cdrl1
        if len(cdrl3_num)>0 and len(cdrl1_num)>0:
            for pp in cdrl3_num:
                #print(pp)
                for pp1 in cdrl1_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        #cdrl3与cdrl2
        if len(cdrl3_num)>0 and len(cdrl2_num)>0:
            for pp in cdrl3_num:
                #print(pp)
                for pp1 in cdrl2_num:
                    #print(pp,pp1)
                    adj11[pp,pp1]=1
                    
        ##################################################### 
        #print(adj11)
        #for ii in np.nditer(dis_rec_rec,op_flags=['readwrite']):
        #    if ii >4.5:
        #        ii[...]=0
        #    else:
        #        ii[...]=1           

        for ii in np.nditer(dis_lig_lig,op_flags=['readwrite']):
            if ii >4.5:
                ii[...]=0
            if 4.5>=ii>0:
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
            
        shu1=receptor_count
        shu2=ligand_count
        shu3=shu1+shu2
        #print(shu1,shu2)

        adj= np.zeros((shu3, shu3))
        
        
        adj[:receptor_count, :receptor_count] = adj11
        adj[receptor_count:, receptor_count:] = dis_lig_lig
        adj[:receptor_count, receptor_count:] = dis_rec_lig
        adj[receptor_count:, :receptor_count] = dis_lig_rec
        #print(adj)
        adj=adj+np.eye(shu3)
            
        target=[]
        txt_list=os.listdir(txt_path)
        txt_list.sort()
        txt=txt_list[nu]
        print(txt)
        txt_name=os.path.join(txt_path,txt)
        with open(txt_name,'r') as tmp_file:
            line=tmp_file.readline()
            line=line.strip()  
            aim_tmp=int(line)
            target.append(aim_tmp)
        
        input_file=os.path.join(savepath,"%s.npz"%pdbid)
        np.savez(input_file,  H=H, A1=adj, T=target)
        print('--------------------------end----------------------')
        
