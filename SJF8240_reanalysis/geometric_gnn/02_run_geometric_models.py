#!/usr/bin/env python3
"""
Geometric GNN benchmark for PROTAC-induced kinase degradation.

Compares 4 architectures on 8-target LOPO-CV (corrected labels: EPHA2=1):
  1. BASELINE   — GATConv (contact topology only, no 3D)
  2. GeomGNN    — SchNet-style invariant GNN (Gaussian RBF on Cα distances)
  3. EGNN       — E(3)-equivariant GNN (Satorras et al. 2021)
  4. ViSNet     — SE(3)-equivariant (PyG ViSNet representation block)

Input graphs: corrected_md_only_v1/inputs/elob/*.npz (EPHA2=1 corrected labels)
Positions:    geometric_gnn/positions/*.npz  (Cα XYZ from pRosettaC rank1 PDBs)

Outputs:
  geometric_gnn/results/geometric_comparison.tsv  — per-target scores, all models
  geometric_gnn/results/geometric_summary.tsv     — per-model AUROC/accuracy
  geometric_gnn/logs/run.log
"""
import sys, json, time, random, math, warnings
import numpy as np
import pandas as pd
from pathlib import Path

warnings.filterwarnings('ignore')

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    from torch_geometric.data import Data
    from torch_geometric.nn import GATConv, global_mean_pool, global_max_pool
    from torch_geometric.nn.models import ViSNet
    from torch_geometric.utils import scatter
    from sklearn.metrics import roc_auc_score
    HAS_DEPS = True
except ImportError as e:
    print(f"[ERROR] Missing dependency: {e}")
    print("Run with: ~/.conda/envs/embeddings/bin/python3 02_run_geometric_models.py")
    sys.exit(1)

SCRIPT_DIR  = Path(__file__).resolve().parent
GRAPH_DIR   = SCRIPT_DIR.parent / 'corrected_md_only_v1/inputs/elob'
POS_DIR     = SCRIPT_DIR / 'positions'
RESULTS_DIR = SCRIPT_DIR / 'results'
LOG_PATH    = SCRIPT_DIR / 'logs/run.log'
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
LOG_PATH.parent.mkdir(parents=True, exist_ok=True)

LABELS = {
    'DDR2': 1, 'ABL1': 0, 'AXL': 0, 'EPHA2': 1,
    'EPHB2': 1, 'MAPK14_full': 1, 'MET': 1, 'RIPK2_full': 1
}
POI_ORDER = ['DDR2', 'ABL1', 'AXL', 'EPHA2', 'EPHB2', 'MAPK14_full', 'MET', 'RIPK2_full']

EPOCHS      = 400
LR          = 1e-3
HIDDEN      = 64
N_RUNS      = 5          # multiple runs per fold to reduce variance
RBF_CUTOFF  = 15.0       # Å — build edges for all Cα pairs within this distance
N_RBF       = 20         # Gaussian RBF basis functions for distance encoding
SEED        = 42

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


# ── Gaussian RBF edge builder ─────────────────────────────────────────────────

def gaussian_rbf(d, centers, widths):
    """Expand distances d (E,) into RBF features (E, n_rbf)."""
    return torch.exp(-widths * (d.unsqueeze(1) - centers.unsqueeze(0))**2)


def build_rbf_params(cutoff=RBF_CUTOFF, n_rbf=N_RBF):
    centers = torch.linspace(0, cutoff, n_rbf)
    widths  = torch.full((n_rbf,), (2 * cutoff / n_rbf) ** -2)
    return centers, widths


def build_distance_edges(pos, cutoff=RBF_CUTOFF):
    """
    Build edge_index + pairwise distances for all node pairs within cutoff.
    Returns: edge_index (2, E), distances (E,)
    """
    n = pos.size(0)
    diff = pos.unsqueeze(1) - pos.unsqueeze(0)   # (N, N, 3)
    dist = diff.norm(dim=-1)                      # (N, N)
    mask = (dist < cutoff) & (dist > 0)
    src, dst = mask.nonzero(as_tuple=True)
    return torch.stack([src, dst], dim=0), dist[src, dst]


# ── Graph loading + augmentation ──────────────────────────────────────────────

def load_graphs():
    """
    Load 8-target corrected contact graphs + 3D Cα positions.
    Returns list of dicts with keys: poi, label, x, edge_index, edge_attr, pos
    """
    graphs = []
    for poi in POI_ORDER:
        g_path = GRAPH_DIR / f'{poi}.npz'
        p_path = POS_DIR  / f'{poi}.npz'
        if not g_path.exists() or not p_path.exists():
            print(f"[WARN] Missing graph or positions for {poi}")
            continue

        g = np.load(g_path)
        p = np.load(p_path)

        x          = torch.tensor(g['x'],          dtype=torch.float32)
        edge_index = torch.tensor(g['edge_index'],  dtype=torch.long)
        edge_attr  = torch.tensor(g['edge_attr'],   dtype=torch.float32)
        label      = int(g['label'][0])

        # Assemble positions: VHL nodes first, then kinase nodes
        pos = np.concatenate([p['vhl_pos'], p['kin_pos']], axis=0)
        pos = torch.tensor(pos, dtype=torch.float32)

        # Center positions at their centroid (improves numerical stability)
        pos = pos - pos.mean(dim=0, keepdim=True)

        # z vector: 1=VHL, 2=kinase (for ViSNet atomic-number embedding)
        n_vhl = len(p['vhl_resids'])
        n_kin = len(p['kin_resids'])
        z = torch.cat([torch.ones(n_vhl, dtype=torch.long),
                       torch.full((n_kin,), 2, dtype=torch.long)])

        # Node feature augmentation: add per-residue PROTAC contact frequency
        # from edge_attr (mean over all edges touching each node)
        n = x.size(0)
        contact_freq_node = torch.zeros(n, 1)
        if edge_index.size(1) > 0:
            src, dst = edge_index
            freq = edge_attr[:, 0:1]  # first column = contact_freq
            contact_freq_node.scatter_add_(0, dst.unsqueeze(1).expand_as(freq), freq)
            cnt = torch.zeros(n, 1)
            cnt.scatter_add_(0, dst.unsqueeze(1), torch.ones(freq.size(0), 1))
            cnt = cnt.clamp(min=1)
            contact_freq_node = contact_freq_node / cnt

        x_aug = torch.cat([x, contact_freq_node], dim=1)   # (N, 4)

        graphs.append({
            'poi':        poi,
            'label':      label,
            'x':          x_aug,
            'x_base':     x,
            'z':          z,
            'pos':        pos,
            'edge_index': edge_index,
            'edge_attr':  edge_attr,
            'n_vhl':      n_vhl,
            'n_kin':      n_kin,
        })

    print(f"Loaded {len(graphs)} graphs")
    return graphs


# ── Model 1: Baseline GATConv (topology only) ─────────────────────────────────

class BaselineGAT(nn.Module):
    def __init__(self, in_ch=3, edge_dim=3, hidden=HIDDEN, dropout=0.3):
        super().__init__()
        heads = 4
        self.edge_proj = nn.Linear(edge_dim, hidden)
        self.conv1 = GATConv(in_ch, hidden // heads, heads=heads,
                             edge_dim=hidden, dropout=dropout)
        self.bn1 = nn.BatchNorm1d(hidden)
        self.conv2 = GATConv(hidden, hidden, heads=1,
                             edge_dim=hidden, dropout=dropout, concat=False)
        self.bn2 = nn.BatchNorm1d(hidden)
        self.dropout = nn.Dropout(dropout)
        self.clf = nn.Sequential(
            nn.Linear(hidden * 2, hidden), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(hidden, 1)
        )

    def forward(self, data):
        x, ei, ea = data.x_base, data.edge_index, data.edge_attr
        batch = data.batch if hasattr(data, 'batch') and data.batch is not None \
            else torch.zeros(x.size(0), dtype=torch.long, device=x.device)
        ea = F.relu(self.edge_proj(ea))
        x  = F.relu(self.bn1(self.conv1(x, ei, ea)))
        x  = self.dropout(x)
        x  = F.relu(self.bn2(self.conv2(x, ei, ea)))
        g  = torch.cat([global_mean_pool(x, batch), global_max_pool(x, batch)], 1)
        return self.clf(g).squeeze(-1)


# ── Model 2: GeomGNN (SchNet-style, invariant) ────────────────────────────────

class RBFLayer(nn.Module):
    def __init__(self, cutoff=RBF_CUTOFF, n_rbf=N_RBF):
        super().__init__()
        self.register_buffer('centers', torch.linspace(0, cutoff, n_rbf))
        gamma = (cutoff / n_rbf) ** -2
        self.register_buffer('widths', torch.full((n_rbf,), gamma))

    def forward(self, d):
        return torch.exp(-self.widths * (d.unsqueeze(1) - self.centers.unsqueeze(0))**2)


class GeomConv(nn.Module):
    """Single SchNet-style interaction block using distance + contact features."""
    def __init__(self, hidden, n_rbf=N_RBF, edge_feat_dim=1):
        super().__init__()
        edge_in = n_rbf + edge_feat_dim
        self.edge_mlp = nn.Sequential(
            nn.Linear(edge_in, hidden), nn.SiLU(),
            nn.Linear(hidden, hidden)
        )
        self.node_mlp = nn.Sequential(
            nn.Linear(hidden * 2, hidden), nn.SiLU(),
            nn.Linear(hidden, hidden)
        )
        self.res_proj = nn.Linear(hidden, hidden)

    def forward(self, h, edge_index, rbf_feat, edge_feat):
        row, col = edge_index
        ef = torch.cat([rbf_feat, edge_feat], dim=-1)
        msg = self.edge_mlp(ef) * h[row]           # filter-based gating
        agg = scatter(msg, col, dim=0, dim_size=h.size(0), reduce='sum')
        h_new = self.node_mlp(torch.cat([h, agg], dim=-1))
        return h_new + self.res_proj(h)             # residual


class GeomGNN(nn.Module):
    def __init__(self, in_ch=4, hidden=HIDDEN, n_rbf=N_RBF, n_layers=3, dropout=0.3):
        super().__init__()
        self.rbf = RBFLayer(n_rbf=n_rbf)
        self.node_embed = nn.Linear(in_ch, hidden)
        self.convs = nn.ModuleList([
            GeomConv(hidden, n_rbf=n_rbf, edge_feat_dim=1) for _ in range(n_layers)
        ])
        self.dropout = nn.Dropout(dropout)
        self.clf = nn.Sequential(
            nn.Linear(hidden * 2, hidden), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(hidden, 1)
        )
        self.cutoff = RBF_CUTOFF

    def forward(self, data):
        pos, x = data.pos, data.x
        batch = data.batch if hasattr(data, 'batch') and data.batch is not None \
            else torch.zeros(x.size(0), dtype=torch.long, device=x.device)

        # Build distance-based edges
        edge_index, dists = build_distance_edges(pos, self.cutoff)
        edge_index = edge_index.to(device)
        dists = dists.to(device)

        rbf_feat = self.rbf(dists)                  # (E, n_rbf)

        # Contact frequency per edge from original edges (if available)
        n = x.size(0)
        contact_node = x[:, 3:4]                    # contact_freq augmented in load
        row, col = edge_index
        edge_contact = (contact_node[row] + contact_node[col]) / 2  # (E, 1)

        h = F.relu(self.node_embed(x))
        for conv in self.convs:
            h = F.relu(conv(h, edge_index, rbf_feat, edge_contact))
            h = self.dropout(h)

        g = torch.cat([global_mean_pool(h, batch), global_max_pool(h, batch)], 1)
        return self.clf(g).squeeze(-1)


# ── Model 3: EGNN (E(3)-equivariant) ─────────────────────────────────────────

class EGNNLayer(nn.Module):
    """EGNN layer: Satorras et al., NeurIPS 2021."""
    def __init__(self, hidden, edge_feat_dim=1):
        super().__init__()
        self.edge_mlp = nn.Sequential(
            nn.Linear(2 * hidden + 1 + edge_feat_dim, hidden), nn.SiLU(),
            nn.Linear(hidden, hidden), nn.SiLU()
        )
        self.node_mlp = nn.Sequential(
            nn.Linear(2 * hidden, hidden), nn.SiLU(),
            nn.Linear(hidden, hidden)
        )
        self.coord_mlp = nn.Sequential(
            nn.Linear(hidden, hidden), nn.SiLU(),
            nn.Linear(hidden, 1, bias=False)
        )
        self.res = nn.Linear(hidden, hidden)
        # Gate for coordinate updates (suppress large moves)
        self.coord_gate = nn.Linear(hidden, 1)

    def forward(self, h, x, edge_index, edge_attr):
        row, col = edge_index
        d_sq = ((x[row] - x[col]) ** 2).sum(dim=-1, keepdim=True)   # (E, 1)
        edge_in = torch.cat([h[row], h[col], d_sq, edge_attr], dim=-1)
        m = self.edge_mlp(edge_in)                                     # (E, hidden)

        # Coordinate update (equivariant: weighted sum of displacement vectors)
        coord_diff  = x[row] - x[col]                                  # (E, 3)
        norm        = coord_diff.norm(dim=-1, keepdim=True).clamp(min=1e-6)
        unit_diff   = coord_diff / norm
        coord_w     = torch.tanh(self.coord_mlp(m))                   # (E, 1)
        x_update    = scatter(coord_w * unit_diff, col, dim=0,
                              dim_size=h.size(0), reduce='mean')
        x_new       = x + x_update

        # Node feature update
        m_agg = scatter(m, col, dim=0, dim_size=h.size(0), reduce='sum')
        h_new = self.node_mlp(torch.cat([h, m_agg], dim=-1))
        h_new = h_new + self.res(h)
        return h_new, x_new


class EGNN(nn.Module):
    def __init__(self, in_ch=4, hidden=HIDDEN, n_layers=3, dropout=0.3):
        super().__init__()
        self.node_embed = nn.Linear(in_ch, hidden)
        self.layers = nn.ModuleList([EGNNLayer(hidden) for _ in range(n_layers)])
        self.dropout = nn.Dropout(dropout)
        self.clf = nn.Sequential(
            nn.Linear(hidden * 2, hidden), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(hidden, 1)
        )
        self.cutoff = RBF_CUTOFF

    def forward(self, data):
        pos, x = data.pos, data.x
        batch = data.batch if hasattr(data, 'batch') and data.batch is not None \
            else torch.zeros(x.size(0), dtype=torch.long, device=x.device)

        edge_index, dists = build_distance_edges(pos, self.cutoff)
        edge_index = edge_index.to(device)
        dists = dists.to(device)

        n = x.size(0)
        contact_node = x[:, 3:4]
        row, col = edge_index
        edge_contact = (contact_node[row] + contact_node[col]) / 2

        h = F.relu(self.node_embed(x))
        cur_pos = pos.clone()
        for layer in self.layers:
            h, cur_pos = layer(h, cur_pos, edge_index, edge_contact)
            h = self.dropout(h)

        g = torch.cat([global_mean_pool(h, batch), global_max_pool(h, batch)], 1)
        return self.clf(g).squeeze(-1)


# ── Model 4: ViSNet-based (SE(3)-equivariant via PyG) ─────────────────────────

class ViSNetClassifier(nn.Module):
    """
    Wraps PyG's ViSNet representation block + classification head.
    Uses ViSNet internals for SE(3)-equivariant node features, then
    applies global mean pooling + linear classifier.
    """
    def __init__(self, hidden=HIDDEN, num_layers=3, max_z=3, cutoff=12.0, dropout=0.3):
        super().__init__()
        self.visnet = ViSNet(
            lmax=1,
            num_layers=num_layers,
            hidden_channels=hidden,
            cutoff=cutoff,
            max_z=max_z,
            reduce_op='sum',
            derivative=False,
        )
        self.clf = nn.Sequential(
            nn.Linear(hidden, hidden), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(hidden, 1)
        )
        self.dropout = nn.Dropout(dropout)

    def forward(self, data):
        z   = data.z.to(device)
        pos = data.pos.to(device)
        batch = data.batch if hasattr(data, 'batch') and data.batch is not None \
            else torch.zeros(z.size(0), dtype=torch.long, device=z.device)

        # Extract node features from ViSNet representation block
        x, v = self.visnet.representation_model(z, pos, batch)
        # x: (N, hidden) scalar features; v: (N, 3, hidden) equivariant vectors
        # For classification, use mean + norm of vectors as graph repr
        v_norm = v.norm(dim=1)                           # (N, hidden)
        node_feat = x + v_norm                           # merge scalar + vector info
        node_feat = self.dropout(node_feat)
        g = global_mean_pool(node_feat, batch)           # (1, hidden) for single graph
        return self.clf(g).squeeze(-1)


# ── Training utils ────────────────────────────────────────────────────────────

def to_device(graph, dev):
    """Move graph tensors to device, return as PyG Data-like namespace."""
    class G: pass
    g = G()
    for k, v in graph.items():
        if isinstance(v, torch.Tensor):
            setattr(g, k, v.to(dev))
        else:
            setattr(g, k, v)
    g.batch = torch.zeros(g.x.size(0), dtype=torch.long, device=dev)
    return g


def train_one_fold(model_cls, train_graphs, test_graph, seed=42, model_kwargs=None):
    """
    Train model_cls on train_graphs, evaluate on test_graph.
    Returns (pred_score, true_label).
    """
    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)

    model_kwargs = model_kwargs or {}
    model = model_cls(**model_kwargs).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, EPOCHS, eta_min=1e-5)

    # Class weights: 6:2 imbalance
    pos_weight = torch.tensor([3.0], device=device)  # weight positives less

    model.train()
    for epoch in range(EPOCHS):
        for g in train_graphs:
            gd = to_device(g, device)
            optimizer.zero_grad()
            logit = model(gd)
            label = torch.tensor([float(g['label'])], device=device)
            loss = F.binary_cross_entropy_with_logits(logit, label, pos_weight=pos_weight)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
        scheduler.step()

    model.eval()
    with torch.no_grad():
        gd = to_device(test_graph, device)
        logit = model(gd)
        score = torch.sigmoid(logit).item()

    return score, test_graph['label']


def run_lopo_cv(model_cls, graphs, model_kwargs=None, name=''):
    """LOPO-CV with N_RUNS per fold for variance reduction."""
    print(f"\n  Running {name} ({N_RUNS} runs/fold) ...")
    results = []
    t0 = time.time()
    for i, test_g in enumerate(graphs):
        train_gs = [g for j, g in enumerate(graphs) if j != i]
        scores = []
        for run in range(N_RUNS):
            sc, lab = train_one_fold(model_cls, train_gs, test_g,
                                     seed=SEED + run * 100, model_kwargs=model_kwargs)
            scores.append(sc)
        mean_score = float(np.mean(scores))
        results.append({'poi': test_g['poi'], 'label': test_g['label'],
                        'score': mean_score, 'model': name})
        correct = (round(mean_score) == test_g['label'])
        print(f"    {test_g['poi']:18s} y={test_g['label']}  "
              f"score={mean_score:.4f}  {'✓' if correct else '✗'}")

    dt = time.time() - t0
    scores_arr = [r['score'] for r in results]
    labels_arr = [r['label'] for r in results]
    n_correct = sum(round(s) == l for s, l in zip(scores_arr, labels_arr))

    try:
        auroc = roc_auc_score(labels_arr, scores_arr)
    except Exception:
        auroc = float('nan')

    print(f"    AUROC={auroc:.4f}  Accuracy={n_correct}/{len(results)}  "
          f"({dt:.0f}s)")
    return results, auroc, n_correct


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("GEOMETRIC GNN BENCHMARK — PROTAC LOPO-CV")
    print(f"Device: {device}  |  Epochs: {EPOCHS}  |  Runs/fold: {N_RUNS}")
    print("=" * 70)

    graphs = load_graphs()
    if len(graphs) < 8:
        print(f"[ERROR] Expected 8 graphs, got {len(graphs)}")
        sys.exit(1)

    # Verify labels
    for g in graphs:
        assert g['label'] == LABELS[g['poi']], f"Label mismatch for {g['poi']}"
    print(f"Labels verified: "
          f"{sum(g['label']==1 for g in graphs)} degraders, "
          f"{sum(g['label']==0 for g in graphs)} non-degraders\n")

    all_results = []
    summary_rows = []

    # ── 1. Baseline GATConv ──
    res, auroc, n_ok = run_lopo_cv(
        BaselineGAT, graphs,
        model_kwargs={'in_ch': 3, 'edge_dim': 3},
        name='Baseline_GAT'
    )
    all_results.extend(res)
    summary_rows.append({'model': 'Baseline_GAT', 'auroc': auroc,
                         'n_correct': n_ok, 'n_total': 8})

    # ── 2. GeomGNN (SchNet-style) ──
    res, auroc, n_ok = run_lopo_cv(
        GeomGNN, graphs,
        model_kwargs={'in_ch': 4},
        name='GeomGNN_SchNet'
    )
    all_results.extend(res)
    summary_rows.append({'model': 'GeomGNN_SchNet', 'auroc': auroc,
                         'n_correct': n_ok, 'n_total': 8})

    # ── 3. EGNN ──
    res, auroc, n_ok = run_lopo_cv(
        EGNN, graphs,
        model_kwargs={'in_ch': 4},
        name='EGNN'
    )
    all_results.extend(res)
    summary_rows.append({'model': 'EGNN', 'auroc': auroc,
                         'n_correct': n_ok, 'n_total': 8})

    # ── 4. ViSNet ──
    res, auroc, n_ok = run_lopo_cv(
        ViSNetClassifier, graphs,
        model_kwargs={'hidden': HIDDEN},
        name='ViSNet_SE3'
    )
    all_results.extend(res)
    summary_rows.append({'model': 'ViSNet_SE3', 'auroc': auroc,
                         'n_correct': n_ok, 'n_total': 8})

    # ── Save results ──
    df_all = pd.DataFrame(all_results)
    df_sum = pd.DataFrame(summary_rows)

    # Pivot to wide format: rows=POI, cols=model
    df_wide = df_all.pivot(index='poi', columns='model', values='score')
    df_wide['label'] = [LABELS[p] for p in df_wide.index]
    df_wide = df_wide[['label'] + [c for c in df_wide.columns if c != 'label']]

    out_comp = RESULTS_DIR / 'geometric_comparison.tsv'
    out_summ = RESULTS_DIR / 'geometric_summary.tsv'
    df_wide.to_csv(out_comp, sep='\t', float_format='%.4f')
    df_sum.to_csv(out_summ, sep='\t', index=False, float_format='%.4f')

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\n{'Model':<20} {'AUROC':>8} {'Correct':>10}")
    print('-' * 42)
    for r in summary_rows:
        print(f"{r['model']:<20} {r['auroc']:>8.4f} {r['n_correct']:>5}/{r['n_total']}")

    print(f"\nReference (E2-GNN, historical best): AUROC ~0.85, 7/8 correct")
    print(f"\nResults:    {out_comp}")
    print(f"Summary:    {out_summ}")

    # Save log
    with open(LOG_PATH, 'w') as f:
        json.dump({
            'run_date': time.strftime('%Y-%m-%dT%H:%M:%S'),
            'epochs': EPOCHS, 'n_runs': N_RUNS, 'hidden': HIDDEN,
            'device': str(device),
            'summary': summary_rows
        }, f, indent=2)


if __name__ == '__main__':
    main()
