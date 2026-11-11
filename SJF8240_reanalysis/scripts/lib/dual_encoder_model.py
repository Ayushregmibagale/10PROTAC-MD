"""Original 04d class definitions; source retained in provenance. No training entrypoint."""
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.nn import GATConv,global_mean_pool,global_max_pool

class EnsembleEncoder(nn.Module):

    def __init__(self, in_ch: int=9, edge_dim: int=3, hidden: int=64, heads: int=4, dropout: float=0.3):
        super().__init__()
        self.edge_proj = nn.Linear(edge_dim, hidden)
        self.conv1 = GATConv(in_ch, hidden // heads, heads=heads, edge_dim=hidden, dropout=dropout, concat=True)
        self.bn1 = nn.BatchNorm1d(hidden)
        self.conv2 = GATConv(hidden, hidden // heads, heads=heads, edge_dim=hidden, dropout=dropout, concat=True)
        self.bn2 = nn.BatchNorm1d(hidden)
        self.conv3 = GATConv(hidden, hidden, heads=1, edge_dim=hidden, dropout=dropout, concat=False)
        self.bn3 = nn.BatchNorm1d(hidden)
        self.dropout = nn.Dropout(dropout)
        self.out_dim = hidden * 2

    def forward(self, data: Data) -> torch.Tensor:
        x, ei, ea = (data.x, data.edge_index, data.edge_attr)
        batch = getattr(data, 'batch', None)
        if batch is None:
            batch = torch.zeros(x.size(0), dtype=torch.long, device=x.device)
        ea_emb = F.relu(self.edge_proj(ea))
        x = F.relu(self.bn1(self.conv1(x, ei, ea_emb)))
        x = self.dropout(x)
        x = F.relu(self.bn2(self.conv2(x, ei, ea_emb)))
        x = self.dropout(x)
        x = F.relu(self.bn3(self.conv3(x, ei, ea_emb)))
        return torch.cat([global_mean_pool(x, batch), global_max_pool(x, batch)], dim=1)

class MDEncoder(nn.Module):

    def __init__(self, in_ch: int=3, edge_dim: int=3, hidden: int=64, heads: int=4, dropout: float=0.3):
        super().__init__()
        self.edge_proj = nn.Linear(edge_dim, hidden)
        self.conv1 = GATConv(in_ch, hidden // heads, heads=heads, edge_dim=hidden, dropout=dropout, concat=True)
        self.bn1 = nn.BatchNorm1d(hidden)
        self.conv2 = GATConv(hidden, hidden, heads=1, edge_dim=hidden, dropout=dropout, concat=False)
        self.bn2 = nn.BatchNorm1d(hidden)
        self.dropout = nn.Dropout(dropout)
        self.out_dim = hidden * 2

    def forward(self, data: Data) -> torch.Tensor:
        x, ei, ea = (data.x, data.edge_index, data.edge_attr)
        batch = getattr(data, 'batch', None)
        if batch is None:
            batch = torch.zeros(x.size(0), dtype=torch.long, device=x.device)
        ea_emb = F.relu(self.edge_proj(ea))
        x = F.relu(self.bn1(self.conv1(x, ei, ea_emb)))
        x = self.dropout(x)
        x = F.relu(self.bn2(self.conv2(x, ei, ea_emb)))
        return torch.cat([global_mean_pool(x, batch), global_max_pool(x, batch)], dim=1)

class DualEncoderGNN(nn.Module):
    """
    Combined pRosettaC + MD dual encoder.

    For POIs without MD (has_md=False), Encoder B output is replaced with zeros.
    The fusion MLP sees [emb_A || emb_B] where emb_B = 0 for non-MD POIs.
    The model learns to rely on emb_A alone in that case.

    Gating: an additional learned gate (sigmoid) weighs the MD contribution.
    When has_md=False, gate is forced to 0 during training/inference.
    """

    def __init__(self, static_in: int=9, static_edge: int=3, dynamic_in: int=3, dynamic_edge: int=3, hidden: int=64, heads: int=4, dropout: float=0.3):
        super().__init__()
        self.enc_a = EnsembleEncoder(static_in, static_edge, hidden, heads, dropout)
        self.enc_b = MDEncoder(dynamic_in, dynamic_edge, hidden, heads, dropout)
        self.md_gate = nn.Parameter(torch.tensor(0.5))
        fused_dim = self.enc_a.out_dim + self.enc_b.out_dim
        self.fusion = nn.Sequential(nn.Linear(fused_dim, hidden * 2), nn.LayerNorm(hidden * 2), nn.ReLU(), nn.Dropout(dropout), nn.Linear(hidden * 2, hidden), nn.ReLU(), nn.Linear(hidden, 1))

    def forward(self, static_graph: Data, dynamic_graph: Data | None) -> torch.Tensor:
        emb_a = self.enc_a(static_graph)
        if dynamic_graph is not None:
            emb_b = self.enc_b(dynamic_graph)
            gate = torch.sigmoid(self.md_gate)
            emb_b = gate * emb_b
        else:
            emb_b = torch.zeros(1, self.enc_b.out_dim, device=emb_a.device)
        fused = torch.cat([emb_a, emb_b], dim=1)
        return self.fusion(fused).squeeze(1)
