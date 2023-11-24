###                                   by scc in 2023/10/9
########################################################################################################################
#######                                                                                                     ############
#######                               using trained model to test new data                                  ############
#######                                                                                                     ############
########################################################################################################################
import matplotlib.pyplot as plt
import networkx as nx
from torch_geometric.utils.convert import to_networkx
import torch
from torch_geometric.data import Dataset
from torch_geometric.data import download_url
import os
from torch_geometric.io import read_planetoid_data
from torch_geometric.datasets import Planetoid
import numpy as np
np.set_printoptions(threshold=np.inf)
from torch_geometric.data import Data
from torch.nn import Linear
import numpy as np
import torch.nn.functional as F
from torch_geometric.nn import GCNConv,SAGEConv
import os.path as osp
from torch_geometric.nn import global_mean_pool
import scipy.sparse as sp
import torch_geometric.nn as pyg_nn
from torch_geometric.data import DataLoader
import warnings
warnings.filterwarnings("ignore", category=Warning)
from sklearn.metrics import roc_auc_score
np.set_printoptions(suppress=True)

import matplotlib.pyplot as plt
import networkx as nx
from torch_geometric.utils.convert import to_networkx
import torch
from torch_geometric.data import Dataset
from torch_geometric.data import download_url
import os
from torch_geometric.io import read_planetoid_data
from torch_geometric.datasets import Planetoid
import numpy as np
from torch_geometric.data import Data
from torch.nn import Linear
import numpy as np
import torch.nn.functional as F
from torch_geometric.nn import GCNConv
import os.path as osp
from torch_geometric.nn import global_mean_pool
import scipy.sparse as sp
import torch_geometric.nn as pyg_nn
from torch_geometric.data import DataLoader
import warnings
warnings.filterwarnings("ignore", category=Warning)



file='ditectory/npz'

class MyOwnDataset(Dataset):
    def __init__(self, root, transform=None, pre_transform=None):
        super().__init__(root, transform, pre_transform)
    @property
    def raw_file_names(self):
        return []
    @property
    def processed_file_names(self):
        return ['datas0.pt','datas1.pt']
    def download(self):
        pass
    def process(self):
        #
        data_list=[]
        for i in range(number):
            path_list=os.listdir(file)
            path_list.sort()
            npz=path_list[i]
            print(npz)
            #pdb2=path_list[i]
            aa=os.path.join(file,npz)
            fileload=np.load(aa)  
            key=list(fileload.keys())
            node=fileload['H'] 
            #node = torch.tensor(x) 
            adj=fileload['A1'] 
            
            edge_index_temp = sp.coo_matrix(adj)  
    
            indices = np.vstack((edge_index_temp.row, edge_index_temp.col))
            edge_index = torch.LongTensor(indices)
            
            x = node
            x = torch.FloatTensor(x)
            
            data=Data(x=x, edge_index=edge_index)
            if self.pre_filter is not None and not self.pre_filter(data):
                continue
            if self.pre_transform is not None:
                data = self.pre_transform(data)
            torch.save(data,osp.join(self.processed_dir,'datas{}.pt'.format(i)))
    
    def len(self):
        return number
    

    def get(self, idx):
        data = torch.load(osp.join(self.processed_dir, 'datas{}.pt'.format(idx)))
        return data


scc_data=MyOwnDataset('directory/11') 

     
class Net(torch.nn.Module):
    """构造GCN模型网络"""
    def __init__(self):
        super(Net, self).__init__()
        self.conv1 = SAGEConv(56, 16) 
        self.conv2 = SAGEConv(16, 32) 
        self.conv3 = SAGEConv(32, 64)
        self.conv4 = SAGEConv(64,2)
    def forward(self, x, edge_index, batch):
        # 1. Obtain node embeddings 
        x = self.conv1(x, edge_index)
        x = x.relu()
        x = F.dropout(x, training=self.training)
        x = self.conv2(x, edge_index)
        x = x.relu()
        x = F.dropout(x, training=self.training)
        x = self.conv3(x, edge_index)
        x = x.relu()
        x = F.dropout(x, training=self.training)
        x = self.conv4(x, edge_index)
        x = x.relu()
        x = F.dropout(x, training=self.training)
        x = pyg_nn.global_mean_pool(x, batch) 
        return F.softmax(x, dim=1) 
test_loader = DataLoader(scc_data, batch_size=4000, shuffle=False)
preds=[]
prob_all=[]
label_all=[]
def evaluate(loader): 
    sccmodel=Net()
    sccmodel=torch.load('trained model/application1/2/3.pt',map_location='cpu') 
    sccmodel.eval()
    correct = 0
    with torch.no_grad():
        for data in loader:
            out = sccmodel(data.x, data.edge_index, data.batch)
            prob_all.extend(out[:,1].detach().cpu().numpy()) 
            bb=out.detach().numpy()
            a2=bb[:,1]
            for i in a2:
                print(i) 
testacc=evaluate(test_loader)
