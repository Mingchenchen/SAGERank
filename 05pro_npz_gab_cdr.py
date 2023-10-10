###                                   by scc in 2023/10/9
########################################################################################################################
#######                                                                                                     ############
#######                               produce gab-cdr npz                                   ############
#######                                                                                                     ############
########################################################################################################################


import os
import numpy as np
import os
from rdkit.Chem.rdmolfiles import MolFromPDBFile
from data_processing.Feature_Processing import get_atom_feature
import numpy as np
np.set_printoptions(threshold = np.inf)
from rdkit.Chem.rdmolops import GetAdjacencyMatrix
from scipy.spatial import distance_matrix
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



def Prepare_Input(receptor_path,ligand_path,txtpath,savepath,fasta_path):
    # extract the interface region
    #root_path=os.path.split(structure_path)[0]  
    interface_rec_list=os.listdir(receptor_path)
    interface_rec_list.sort()
    
    #fasta_path='/data2/data_home/ccsun/link/scc_neuralnetwork/Recognition_of_antigens_by_antibodies/fasta/all-fasta-c/1ADQ.C.fasta'
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
    
    
    
    for i,ii in enumerate(interface_rec_list):
        pdbid=ii.replace('.pdb','')
        rec_pdb=os.path.join(receptor_path,ii)
        
        lig_pdb=os.path.join(ligand_path,ii)
        
        #receptor_path, ligand_path = Extract_Interface(aa)  
    
        receptor_mol = MolFromPDBFile(rec_pdb, sanitize=False)
        ligand_mol = MolFromPDBFile(lig_pdb, sanitize=False)   
    
    
        receptor_count = receptor_mol.GetNumAtoms()
        ligand_count = ligand_mol.GetNumAtoms()  
    
    
        receptor_feature = get_atom_feature(receptor_mol, is_ligand=False) 
        ligand_feature = get_atom_feature(ligand_mol, is_ligand=True)
   
    
        # get receptor adj matrix
        c1 = receptor_mol.GetConformers()[0]
        d1 = np.array(c1.GetPositions()) 
    
    
        #adj1 = GetAdjacencyMatrix(receptor_mol) + np.eye(receptor_count)   
        scc1=GetAdjacencyMatrix(receptor_mol)
        print(np.shape(scc1))
        #np.shape(GetAdjacencyMatrix(receptor_mol)) 为（387*387），np.shape(np.eye(receptor_count)) 为（387*387）
    
    
        adj11=np.zeros((receptor_count,receptor_count))  
        #print(adj1)
        cdrh1=list_all[0]
        cdrh2=list_all[1]
        cdrh3=list_all[2]
        cdrl1=list_all[3]
        cdrl2=list_all[4]
        cdrl3=list_all[5]
        with open(rec_pdb,'r') as s1:
            line = s1.readlines() 
            #CDRH1:
            cdrh1_num=[]
            for zz in cdrh1:
                #print(zz)
                #print('---------------')
                for p1,p2 in enumerate(line):
                    p3=int(p2[23:26])
                    p4=str(p2[13:16]).strip() 
                    if  p3==zz and p4=='N':
                        h1_n=p1
                        h1_ca=p1+1
                        h1_c=p1+2
                        h1_o=p1+3
                        #print(h1_n,h1_ca,h1_c,h1_o)
                        cdrh1_num.append(h1_n)
                        cdrh1_num.append(h1_ca)
                        cdrh1_num.append(h1_c)
                        cdrh1_num.append(h1_o)

            cdrh2_num=[]
            for zz in cdrh2:
                #print(zz)
                #print('---------------')
                for p1,p2 in enumerate(line):
                    p3=int(p2[23:26])
                    p4=str(p2[13:16]).strip() 
                    if  p3==zz and p4=='N':
                        h1_n=p1
                        h1_ca=p1+1
                        h1_c=p1+2
                        h1_o=p1+3            
                        #print(h1_n,h1_ca,h1_c,h1_o)
                        cdrh2_num.append(h1_n)
                        cdrh2_num.append(h1_ca)
                        cdrh2_num.append(h1_c)
                        cdrh2_num.append(h1_o)
            #print(cdrh2_num)
            
            #CDRH3:
            cdrh3_num=[]
            for zz in cdrh3:
                #print(zz)
                #print('---------------')
                for p1,p2 in enumerate(line):
                    p3=int(p2[23:26])
                    p4=str(p2[13:16]).strip() 
                    if  p3==zz and p4=='N':
                        h1_n=p1
                        h1_ca=p1+1
                        h1_c=p1+2
                        h1_o=p1+3
                        #print(h1_n,h1_ca,h1_c,h1_o)
                        cdrh3_num.append(h1_n)
                        cdrh3_num.append(h1_ca)
                        cdrh3_num.append(h1_c)
                        cdrh3_num.append(h1_o)
            #print(cdrh3_num)
            
            #CDRL1:
            cdrl1_num=[]
            for zz in cdrl1:
                #print(zz)
                #print('---------------')
                for p1,p2 in enumerate(line):
                    p3=int(p2[23:26])
                    p4=str(p2[13:16]).strip() 
                    if  p3==zz and p4=='N':
                        h1_n=p1
                        h1_ca=p1+1
                        h1_c=p1+2
                        h1_o=p1+3
                        #print(h1_n,h1_ca,h1_c,h1_o)
                        cdrl1_num.append(h1_n)
                        cdrl1_num.append(h1_ca)
                        cdrl1_num.append(h1_c)
                        cdrl1_num.append(h1_o)
            #print(cdrl1_num)
            
            #CDRL2:
            cdrl2_num=[]
            for zz in cdrl1:
                #print(zz)
                #print('---------------')
                for p1,p2 in enumerate(line):
                    p3=int(p2[23:26])
                    p4=str(p2[13:16]).strip() 
                    if  p3==zz and p4=='N':
                        h1_n=p1
                        h1_ca=p1+1
                        h1_c=p1+2
                        h1_o=p1+3
                        #print(h1_n,h1_ca,h1_c,h1_o)
                        cdrl2_num.append(h1_n)
                        cdrl2_num.append(h1_ca)
                        cdrl2_num.append(h1_c)
                        cdrl2_num.append(h1_o)
            #print(cdrl2_num)
            
            #CDRL3:
            cdrl3_num=[]
            for zz in cdrl3:
                #print(zz)
                #print('---------------')
                for p1,p2 in enumerate(line):
                    p3=int(p2[23:26])
                    p4=str(p2[13:16]).strip() 
                    if  p3==zz and p4=='N':
                        h1_n=p1
                        h1_ca=p1+1
                        h1_c=p1+2
                        h1_o=p1+3
                        #print(h1_n,h1_ca,h1_c,h1_o)
                        cdrl3_num.append(h1_n)
                        cdrl3_num.append(h1_ca)
                        cdrl3_num.append(h1_c)
                        cdrl3_num.append(h1_o)
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
        
        #print('-----------------')
        #print(adj11)
        adj111=adj11 + np.eye(receptor_count)             
        #print('-----------------')
        #print(adj111)
        
        
        # get ligand adj matrix
        c2 = ligand_mol.GetConformers()[0]
        d2 = np.array(c2.GetPositions())
        adj2 = GetAdjacencyMatrix(ligand_mol) + np.eye(ligand_count)
    
        # combine                                                                                                                                                                     is
        H = np.concatenate([receptor_feature, ligand_feature], 0)
   
    
        agg_adj1 = np.zeros((receptor_count + ligand_count, receptor_count + ligand_count))
        agg_adj1[:receptor_count, :receptor_count] = adj111
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

        valid = np.zeros((receptor_count + ligand_count,))  #
   
        valid[:receptor_count] = 1
    

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
        
    
        input_file=os.path.join(savepath,"%s.npz"%pdbid)
        #np.savez(input_file,  H=H, A1=agg_adj2)
        np.savez(input_file,  H=H, A1=agg_adj2, T=target)



