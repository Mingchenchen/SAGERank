###                                   by scc in 2023/10/9
########################################################################################################################
#######                                                                                                     ############
#######                                 start training                                                      ############
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
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import precision_recall_curve, average_precision_score
warnings.filterwarnings("ignore", category=Warning)


#file='/data4_large1/home_data/ccsun/scc_neuralnetwork/Recognition_of_antigens_by_antibodies/database/train200/preprocess/npz/gab'
#file1='/data4_large1/home_data/ccsun/scc_neuralnetwork/Recognition_of_antigens_by_antibodies/database/test30/preprocess/npz/gab'

class MyOwnDataset(Dataset):
    def __init__(self, root, transform=None, pre_transform=None):
        super().__init__(root, transform, pre_transform)
   
    @property
    def raw_file_names(self):
        # pass 
        return []
    @property
    def processed_file_names(self):
        return ['datas0.pt','datas1.pt']

    def download(self):
        pass
    
    def len(self):
        return 369932
        #return 10
    def get(self, idx):
        data = torch.load(osp.join(self.processed_dir, 'datas{}.pt'.format(idx)))
        return data 

class MyOwnDataset1(Dataset):
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
    
    
    def len(self):
        return 46835
        #return 5
    def get(self, idx):
        data = torch.load(osp.join(self.processed_dir, 'datas{}.pt'.format(idx)))
        return data

dataset=MyOwnDataset('/data4_large1/home_data/ccsun/scc_neuralnetwork/02/database/train-pt/gab-normal')
dataset1=MyOwnDataset1('/data4_large1/home_data/ccsun/scc_neuralnetwork/02/database/test-pt/gab-normal')
class Net(torch.nn.Module):
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
        return F.log_softmax(x, dim=1) 

model = Net()
device = torch.device('cuda:0') 
model = Net().to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=0.0003) 
train_loader = DataLoader(dataset, batch_size=512, shuffle=True) 
#for data in train_loader:
#    print(data)
#valid_loader = DataLoader(valid_dataset, batch_size=512, shuffle=False)
test_loader = DataLoader(dataset1, batch_size=1000, shuffle=False)
criterion = torch.nn.CrossEntropyLoss()
criterion = criterion.to(device)

def train():
    model.train() 
    loss_all = 0
    for data in train_loader: 
        #print(data.y)
        # print(data.batch.shape)
        # print(data.x.shape)
        data=data.to(device)
        output = model(data.x, data.edge_index, data.batch) 
        label = data.y 
        # print(label)
        loss = criterion(output,label) 
        loss.backward() 
        loss_all += loss.item() 
        optimizer.step() 
        optimizer.zero_grad() 
    train_loss = (loss_all / len(dataset)) 
 
    return train_loss
    
preds=[]
def evaluate(loader): 
    model.eval()
    
    correct = 0
    for data in loader:
        # Iterate in batches over the training/test dataset.
        data=data.to(device)
        out = model(data.x, data.edge_index, data.batch)
        #print(out)
        pred = out.argmax(dim=1)  
        #print(pred)
        aa=data.y
        #print(aa)
        preds.append(pred)
        correct += int((pred == data.y).sum())  # Check against ground-truth labels.
        #print(correct)
    return correct / len(loader.dataset)  # Derive ratio of correct predictions.

prob_all=[]
label_all=[]
def evaluate1(loader): 
    model.eval()
    
    correct = 0
    for data in loader:
        # Iterate in batches over the training/test dataset.
        data=data.to(device)
        out = model(data.x, data.edge_index, data.batch)
        
        aa=data.y.detach().cpu()
        prob_all.extend(out[:,1].detach().cpu().numpy()) 
        label_all.extend(aa)
        
    return roc_auc_score(label_all,prob_all)# Derive ratio of correct pre

a1=[]
a2=[]
def evaluate2(loader): 
    model.eval()
    
    correct = 0
    for data in loader:
        # Iterate in batches over the training/test dataset.
        data=data.to(device)
        out = model(data.x, data.edge_index, data.batch)
        #print(out)
        pred = out.argmax(dim=1).detach().cpu()  
        #print(pred)
        a1.extend(pred)
        
        #print(a1)
        aa=data.y.detach().cpu()
        a2.extend(aa)
    #print(a1)
	    
    return precision_score(a2,a1), recall_score(a2,a1),  f1_score(aa, pred)  
train_loss_all = []
#valid_loss_all = []

train_acc_all=[]
#valid_acc_all=[]
test_acc_all=[]

for epoch in range(100):
    train_loss=train()	
    train_acc = evaluate(train_loader)
    train_auc = evaluate1(train_loader)
    a3=evaluate2(train_loader)
    train_recision = a3[0]
    train_recall = a3[1]
    train_f1_score = a3[2]
    test_acc = evaluate(test_loader)
    test_auc = evaluate1(test_loader)
    a4=evaluate2(test_loader)
    test_recision = a4[0]
    test_recall = a4[1]
    test_f1_score = a4[2]
    t_acc=str(test_acc)
    if test_acc >= 0.7:
        torch.save(model, "/data4_large1/home_data/ccsun/scc_neuralnetwork/02/train/gab/savemodel/%d_model.pt"%epoch)
    print(f'Epoch: {epoch:03d}, train loss: {train_loss:.4f},train Acc: {train_acc:.4f},train AUC: {train_auc:.4f},test Acc: {test_acc:.4f},test AUC: {test_auc:.4f}')



