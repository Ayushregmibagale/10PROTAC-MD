"""Exact MDEncoder AST from original 04h MD-only model; no original entrypoint."""
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GATConv,global_mean_pool,global_max_pool
HIDDEN=64
HEADS=4
DROP=0.3

class MDEncoder(nn.Module):

    def __init__(self):
        super().__init__()
        self.ep = nn.Linear(3, HIDDEN)
        self.c1 = GATConv(3, HIDDEN // HEADS, heads=HEADS, edge_dim=HIDDEN, dropout=DROP, concat=True)
        self.b1 = nn.BatchNorm1d(HIDDEN)
        self.c2 = GATConv(HIDDEN, HIDDEN, heads=1, edge_dim=HIDDEN, dropout=DROP, concat=False)
        self.b2 = nn.BatchNorm1d(HIDDEN)
        self.drop = nn.Dropout(DROP)
        self.head = nn.Sequential(nn.Linear(HIDDEN * 2, HIDDEN), nn.ReLU(), nn.Linear(HIDDEN, 1))

    def forward(self, d):
        x, ei, ea = (d.x, d.edge_index, d.edge_attr)
        b = getattr(d, 'batch', None)
        if b is None:
            b = torch.zeros(x.size(0), dtype=torch.long, device=x.device)
        ea = F.relu(self.ep(ea))
        x = F.relu(self.b1(self.c1(x, ei, ea)))
        x = self.drop(x)
        x = F.relu(self.b2(self.c2(x, ei, ea)))
        emb = torch.cat([global_mean_pool(x, b), global_max_pool(x, b)], dim=1)
        return self.head(emb).squeeze(1)
