###                                   by scc in 2023/10/9
########################################################################################################################
#######                                                                                                     ############
#######                               produce interface                                                     ############
#######                                                                                                     ############
########################################################################################################################



import os
RESIDUE_Forbidden_SET={"FAD"}
def Extract_Interface(pdb_path,interfacepath,interfacepath1):
    """
    specially for 2 docking models
    :param pdb_path:docking model path
    :rcount: receptor atom numbers
    :return:
    extract a receptor and ligand, meanwhile, write two files of the receptor interface part, ligand interface part
    """
    (filepath, filename) = os.path.split(pdb_path)
    #interfacepath='/data2/data_home/ccsun/link/scc_neuralnetwork/test_set/negative_above_positive/3lev/preprocessed/inter/rec'
    #interfacepath1='/data2/data_home/ccsun/link/scc_neuralnetwork/test_set/negative_above_positive/3lev/preprocessed/inter/lig'

    receptor_list=[]
    ligand_list=[]
    rlist=[]
    llist=[]
    count_r=0
    count_l=0
    with open(pdb_path,'r') as file:
        line = file.readline()               # call readline()
        while line[0:4]!='ATOM':
            line=file.readline()
        atomid = 0
        count = 1
        goon = False
        chain_id = line[21]
        residue_type = line[17:20]
        pre_residue_type = residue_type
        tmp_list = []
        pre_residue_id = 0
        pre_chain_id = line[21]
        first_change=True
        while line:

            dat_in = line[0:80].split()
            if len(dat_in) == 0:
                line = file.readline()
                continue

            if (dat_in[0] == 'ATOM'):
                chain_id = line[21]
                residue_id = int(line[23:26])

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                residue_type = line[17:20]
                # First try CA distance of contact map
                atom_type = line[13:16].strip()
                if chain_id=="C":
                    goon=True
                if (goon):
                    if first_change:
                        rlist.append(tmp_list)
                        tmp_list = []
                        tmp_list.append([x, y, z, atom_type, count_l])
                        count_l += 1
                        ligand_list.append(line)
                        first_change=False
                    else:
                        ligand_list.append(line)  # used to prepare write interface region
                        if pre_residue_type == residue_type:
                            tmp_list.append([x, y, z, atom_type, count_l])
                        else:
                            llist.append(tmp_list)
                            tmp_list = []
                            tmp_list.append([x, y, z, atom_type, count_l])
                        count_l += 1
                else:
                    receptor_list.append(line)
                    if pre_residue_type == residue_type:
                        tmp_list.append([x, y, z, atom_type, count_r])
                    else:
                        rlist.append(tmp_list)
                        tmp_list = []
                        tmp_list.append([x, y, z, atom_type, count_r])
                    count_r += 1

                atomid = int(dat_in[1])
                chain_id = line[21]
                count = count + 1
                pre_residue_type = residue_type
                pre_residue_id = residue_id
                pre_chain_id = chain_id
            line = file.readline()
    print("Extracting %d/%d atoms for receptor, %d/%d atoms for ligand"%(len(receptor_list),count_r,len(ligand_list),count_l))
    final_receptor, final_ligand=Form_interface(rlist,llist,receptor_list,ligand_list)
    #write that into our path
    rpath=Write_Interface(final_receptor,interfacepath,"%s"%filename)
    lpath=Write_Interface(final_ligand, interfacepath1,"%s"%filename)
    return rpath,lpath
#@set_timeout(100000, after_timeout)
def Form_interface(rlist,llist,receptor_list,ligand_list,cut_off=10):

    cut_off=cut_off**2
    r_index=set()
    l_index=set()
    for rindex,item1 in enumerate(rlist):
        for lindex,item2 in enumerate(llist):
            min_distance=1000000
            residue1_len=len(item1)
            residue2_len=len(item2)
            for m in range(residue1_len):
                atom1=item1[m]
                for n in range(residue2_len):
                    atom2=item2[n]
                    distance=0
                    for k in range(3):
                        distance+=(atom1[k]-atom2[k])**2
                    #distance=np.linalg.norm(atom1[:3]-atom2[:3])
                    if distance<=min_distance:
                        min_distance=distance
            if min_distance<=cut_off:
                if rindex not in r_index:
                    r_index.add(rindex)
                if lindex not in l_index:
                    l_index.add(lindex)
    r_index=list(r_index)
    l_index=list(l_index)
    newrlist=[]
    for k in range(len(r_index)):
        newrlist.append(rlist[r_index[k]])
    newllist=[]
    for k in range(len(l_index)):
        newllist.append(llist[l_index[k]])
    #print("After filtering the interface region, %d/%d residue in receptor, %d/%d residue in ligand" % (len(newrlist),len(rlist), len(newllist),len(llist)))
    #get the line to write new interface file
    final_receptor=[]
    final_ligand=[]
    #print(newllist)
    list1=[]
    list2=[]
    #print(newrlist)
    for residue in newrlist:
        for tmp_atom in residue:
            our_index=tmp_atom[4]
            list1.append(our_index)
    list1.sort()
    for index in list1:
        #print(index)
        final_receptor.append(receptor_list[index])

    for residue in newllist:
        for tmp_atom in residue:
            our_index=tmp_atom[4]
            list2.append(our_index)
    list2.sort()
    for index in list2:
        final_ligand.append(ligand_list[index])
    
    #print(list1)
    print("In the interface region, %d atoms in receptor, %d atoms in ligand"%(len(final_receptor),len(final_ligand)))

    return final_receptor,final_ligand

def Write_Interface(line_list,pdb_path,ext_file):
    new_path=os.path.join(pdb_path,ext_file)
    #print(line_list)
    with open(new_path,'w') as file:
        for line in line_list:
            #print(line)
            #check residue in the common residue or not. If not, no write for this residue
            residue_type = line[17:20]
            if residue_type in RESIDUE_Forbidden_SET:
                continue
            file.write(line)
    return new_path

def produce_interface(pdb_path,r_path,l_path):
    path_list=os.listdir(pdb_path)
    path_list.sort()
    for i in path_list:
        pdb=i
        pdbid=pdb.replace('.pdb','')
        aa=os.path.join(pdb_path,pdb)
        rec=r_path
        lig=l_path
        Extract_Interface(aa,rec,lig)
        