"""Shared window GAT, mean/std/trend pooling; fixed before evaluation."""
import torch
from torch import nn
from torch.nn import functional as F
from torch_geometric.nn import GATConv,global_mean_pool,global_max_pool
class TemporalEncoder(nn.Module):
 def __init__(self):
  super().__init__();self.ep=nn.Linear(3,64);self.c1=GATConv(34,16,heads=4,edge_dim=64,dropout=.3);self.b1=nn.BatchNorm1d(64);self.c2=GATConv(64,64,heads=1,edge_dim=64,dropout=.3,concat=False);self.b2=nn.BatchNorm1d(64);self.drop=nn.Dropout(.3);self.head=nn.Sequential(nn.Linear(384,64),nn.ReLU(),nn.Linear(64,1))
 def forward(self,d):
  ea=F.relu(self.ep(d.edge_attr));x=self.drop(F.relu(self.b1(self.c1(d.x,d.edge_index,ea))));x=F.relu(self.b2(self.c2(x,d.edge_index,ea)));emb=torch.cat([global_mean_pool(x,d.batch),global_max_pool(x,d.batch)],1)
  assert emb.shape==(4,128),'Exactly four ordered windows per target required'
  return self.head(torch.cat([emb.mean(0),emb.std(0,unbiased=False),emb[-1]-emb[0]])).reshape(1)
