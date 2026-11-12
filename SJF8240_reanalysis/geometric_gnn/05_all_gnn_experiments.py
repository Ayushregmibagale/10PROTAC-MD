#!/usr/bin/env python3
"""
Five GNN experiments for PROTAC degrader vs non-degrader classification.
Pure GNN approach: only MD contact graphs and pRosettaC static graphs. No external features.

Experiments:
  E1  Dual encoder (pRosettaC static 17-POI + MD 8-POI, corrected EPHA2=1)
  E2  Per-pose MD augmentation (16-19 graphs/POI, same GATConv)
  E3  Bridging-conditional dual-view (G_bridge + G_nonbr → fusion)
  E4  Richer MD node features (aa-type std-contact from per-pose variation)
  E5  E2+E3+E4 combined (all graph types, all features)

Run: ~/.conda/envs/embeddings/bin/python3 05_all_gnn_experiments.py
"""
# RETIRED_BY_FOLDER_AUDIT_20261112: original source is preserved in the audit archive.
if __name__ == "__main__":
    raise SystemExit("Retired exploratory entry point: see folder_audit_20261112/AUDIT.md and use RUN_ANALYSIS.py. Prior outputs are preserved but not certified by this script.")

import sys, json, time, random, warnings
import numpy as np
import pandas as pd
from pathlib import Path
from itertools import combinations

warnings.filterwarnings('ignore')

try:
    import torch, torch.nn as nn, torch.nn.functional as F
    from torch_geometric.data import Data
    from torch_geometric.nn import GATConv, global_mean_pool, global_max_pool
    from sklearn.metrics import roc_auc_score
except ImportError as e:
    print(f"[ERROR] {e}\nRun with ~/.conda/envs/embeddings/bin/python3")
    sys.exit(1)

SCRIPT_DIR  = Path(__file__).resolve().parent
GRAPH_DIR_ORIG = SCRIPT_DIR.parent / 'corrected_md_only_v1/inputs/elob'
STATIC_PT   = Path.home() / 'protac-md/analysis/gnn/ensemble_graphs.pt'
BRIDGE_DIR  = SCRIPT_DIR / 'graphs/bridge'
PP_DIR      = SCRIPT_DIR / 'graphs/per_pose'
RESULTS_DIR = SCRIPT_DIR / 'results'
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

LABELS = {
    'DDR2':1,'ABL1':0,'AXL':0,'EPHA2':1,
    'EPHB2':1,'MAPK14_full':1,'MET':1,'RIPK2_full':1,
}
POI_ORDER = ['DDR2','ABL1','AXL','EPHA2','EPHB2','MAPK14_full','MET','RIPK2_full']

# Non-MD POIs with pRosettaC static graphs (for 17-POI dual encoder)
# Excluding TNIK (artifact), MAPK14 (artifact), CDK17 (relabeled), RPS6KB2/SIK3 (unknown)
STATIC_ONLY = {
    'MAP4K5':0,'SLK':0,'CDK2':0,'CDK4':0,'CDK6':0,'DAPK1':0,'IRAK1':0,'TLK2':0
}

EPOCHS  = 500
LR      = 1e-3
HIDDEN  = 64
N_RUNS  = 5
SEED    = 42
device  = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


# ── Data loaders ──────────────────────────────────────────────────────────────

def load_md_graphs(label_override=None):
    """Load 8-target corrected MD contact graphs."""
    graphs = []
    for poi in POI_ORDER:
        p = GRAPH_DIR_ORIG / f'{poi}.npz'
        if not p.exists():
            continue
        d = np.load(p)
        lab = label_override.get(poi, int(d['label'][0])) if label_override else int(d['label'][0])
        graphs.append({'poi': poi, 'label': lab,
                       'x':   torch.tensor(d['x'],          dtype=torch.float32),
                       'ei':  torch.tensor(d['edge_index'],  dtype=torch.long),
                       'ea':  torch.tensor(d['edge_attr'],   dtype=torch.float32)})
    return graphs


def load_static_graphs():
    """Load pRosettaC ensemble graphs for all labelled POIs, correcting EPHA2=1."""
    static = torch.load(STATIC_PT, weights_only=False, map_location='cpu')
    # Correct EPHA2 label and filter to clean POIs
    keep_pois = set(list(LABELS.keys()) + list(STATIC_ONLY.keys()))
    out = {}
    for g in static:
        if g.poi not in keep_pois:
            continue
        lab = LABELS.get(g.poi, STATIC_ONLY.get(g.poi, int(g.y)))
        out[g.poi] = {'poi': g.poi, 'label': lab,
                      'x': g.x, 'ei': g.edge_index, 'ea': g.edge_attr}
    return out


def load_per_pose_graphs():
    """Load per-pose MD graphs. Returns dict: poi → list of graphs."""
    out = {poi: [] for poi in POI_ORDER}
    if not PP_DIR.exists():
        return out
    for f in sorted(PP_DIR.glob('*.npz')):
        # filename: <POI>_pose<NNN>.npz
        parts = f.stem.split('_pose')
        if len(parts) != 2:
            continue
        poi = parts[0]
        if poi not in LABELS:
            continue
        d = np.load(f)
        out[poi].append({'poi': poi, 'label': LABELS[poi],
                         'x':  torch.tensor(d['x'],          dtype=torch.float32),
                         'ei': torch.tensor(d['edge_index'],  dtype=torch.long),
                         'ea': torch.tensor(d['edge_attr'],   dtype=torch.float32)})
    return out


def load_bridge_graphs():
    """Load bridging-conditional graphs. Returns (bridge_dict, nonbr_dict)."""
    br, nb = {}, {}
    if not BRIDGE_DIR.exists():
        return br, nb
    for poi in POI_ORDER:
        fb = BRIDGE_DIR / f'{poi}_bridge.npz'
        fn = BRIDGE_DIR / f'{poi}_nonbr.npz'
        if fb.exists():
            d = np.load(fb)
            br[poi] = {'poi': poi, 'label': LABELS[poi],
                       'x':  torch.tensor(d['x'],          dtype=torch.float32),
                       'ei': torch.tensor(d['edge_index'],  dtype=torch.long),
                       'ea': torch.tensor(d['edge_attr'],   dtype=torch.float32)}
        if fn.exists():
            d = np.load(fn)
            nb[poi] = {'poi': poi, 'label': LABELS[poi],
                       'x':  torch.tensor(d['x'],          dtype=torch.float32),
                       'ei': torch.tensor(d['edge_index'],  dtype=torch.long),
                       'ea': torch.tensor(d['edge_attr'],   dtype=torch.float32)}
    return br, nb


def compute_per_pose_node_std(per_pose_graphs):
    """
    For each POI, compute per-node contact-frequency std across poses.
    Returns dict: poi → std_tensor (N,)
    """
    out = {}
    for poi, poses in per_pose_graphs.items():
        if not poses:
            continue
        n = poses[0]['x'].shape[0]
        # Reconstruct per-node contact freq from edge_attr (col 0 = contact_freq)
        # Average edge freq touching each dst node → contact score
        stacked = []
        for g in poses:
            ei, ea = g['ei'], g['ea']
            node_cf = torch.zeros(n)
            if ei.shape[1] > 0:
                cnt = torch.zeros(n)
                node_cf.scatter_add_(0, ei[1], ea[:, 0])
                cnt.scatter_add_(0, ei[1], torch.ones(ea.shape[0]))
                node_cf = node_cf / cnt.clamp(min=1)
            stacked.append(node_cf)
        out[poi] = torch.stack(stacked, 0).std(0)  # (N,)
    return out


# ── GNN Models ────────────────────────────────────────────────────────────────

class GATEncoder(nn.Module):
    """Standard 2-layer GATConv encoder; in_ch and edge_dim configurable."""
    def __init__(self, in_ch, edge_dim=3, hidden=HIDDEN, heads=4, dropout=0.3):
        super().__init__()
        self.edge_proj = nn.Linear(edge_dim, hidden)
        self.conv1 = GATConv(in_ch, hidden//heads, heads=heads,
                             edge_dim=hidden, dropout=dropout)
        self.bn1   = nn.BatchNorm1d(hidden)
        self.conv2 = GATConv(hidden, hidden, heads=1,
                             edge_dim=hidden, dropout=dropout, concat=False)
        self.bn2   = nn.BatchNorm1d(hidden)
        self.drop  = nn.Dropout(dropout)
        self.out_dim = hidden * 2

    def forward(self, x, ei, ea, batch=None):
        if batch is None:
            batch = torch.zeros(x.size(0), dtype=torch.long, device=x.device)
        ea = F.relu(self.edge_proj(ea))
        x  = F.relu(self.bn1(self.conv1(x, ei, ea)));  x = self.drop(x)
        x  = F.relu(self.bn2(self.conv2(x, ei, ea)))
        return torch.cat([global_mean_pool(x, batch),
                          global_max_pool(x, batch)], dim=1)


def to_dev(g, dev=device):
    return {k: v.to(dev) if isinstance(v, torch.Tensor) else v for k, v in g.items()}


def run_fold(model_fn, train_gs, test_g, seed, epochs=EPOCHS):
    torch.manual_seed(seed); random.seed(seed); np.random.seed(seed)
    model = model_fn().to(device)
    opt   = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, epochs, eta_min=1e-5)
    pw    = torch.tensor([3.0], device=device)

    model.train()
    for _ in range(epochs):
        for g in train_gs:
            gd = to_dev(g)
            opt.zero_grad()
            logit = model(gd)
            lab   = torch.tensor([float(g['label'])], device=device)
            loss  = F.binary_cross_entropy_with_logits(logit, lab, pos_weight=pw)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
        sched.step()

    model.eval()
    with torch.no_grad():
        return torch.sigmoid(model(to_dev(test_g))).item()


def lopo_cv(model_fn, graphs_by_poi, name, epochs=EPOCHS):
    """LOPO-CV. graphs_by_poi: {poi: [graph, ...]}. Single graph or list per poi."""
    pois = list(graphs_by_poi.keys())
    all_labels   = [LABELS[p] for p in pois]
    all_scores   = []

    print(f"\n  [{name}] {N_RUNS} runs/fold ...")
    t0 = time.time()
    for i, test_poi in enumerate(pois):
        # Test: mean over all poses of this POI
        test_gs = graphs_by_poi[test_poi] if isinstance(graphs_by_poi[test_poi], list) \
                  else [graphs_by_poi[test_poi]]

        # Train: all poses of all other POIs
        train_gs = []
        for p in pois:
            if p == test_poi:
                continue
            gs = graphs_by_poi[p]
            if isinstance(gs, list):
                train_gs.extend(gs)
            else:
                train_gs.append(gs)

        scores = []
        for run in range(N_RUNS):
            fold_scores = [run_fold(model_fn, train_gs, tg, seed=SEED+run*100, epochs=epochs)
                           for tg in test_gs]
            scores.append(float(np.mean(fold_scores)))

        sc = float(np.mean(scores))
        all_scores.append(sc)
        correct = round(sc) == LABELS[test_poi]
        print(f"    {test_poi:18s} y={LABELS[test_poi]}  score={sc:.4f}  "
              f"{'✓' if correct else '✗'}")

    n_correct = sum(round(s)==l for s,l in zip(all_scores, all_labels))
    try:
        auroc = roc_auc_score(all_labels, all_scores)
    except Exception:
        auroc = float('nan')
    print(f"    AUROC={auroc:.4f}  Accuracy={n_correct}/{len(pois)}  ({time.time()-t0:.0f}s)")
    return dict(zip(pois, all_scores)), auroc, n_correct


# ═══════════════════════════════════════════════════════════════════════════════
# Experiment 1: Dual Encoder (pRosettaC static 17-POI + MD 8-POI)
# ═══════════════════════════════════════════════════════════════════════════════

class DualEncoder(nn.Module):
    def __init__(self):
        super().__init__()
        self.static_enc = GATEncoder(in_ch=9, edge_dim=3)
        self.md_enc     = GATEncoder(in_ch=3, edge_dim=3)
        fused_dim = self.static_enc.out_dim + self.md_enc.out_dim
        self.clf = nn.Sequential(
            nn.Linear(fused_dim, HIDDEN), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(HIDDEN, 1)
        )

    def forward(self, g):
        # static graph
        s_emb = self.static_enc(g['sx'], g['sei'], g['sea'])
        # MD graph (may be zero-padded)
        m_emb = self.md_enc(g['mx'], g['mei'], g['mea'])
        return self.clf(torch.cat([s_emb, m_emb], dim=1)).squeeze(-1)


def run_e1(static_graphs, md_graphs):
    """Experiment 1: Corrected dual encoder on 17-POI training set."""
    print("\n" + "="*60)
    print("E1: Dual Encoder (static 17-POI + MD 8-POI, EPHA2=1 corrected)")
    print("="*60)

    # Build paired graphs dict: poi → paired_graph
    # For non-MD POIs: zero MD graph matching the MD encoder input size
    def zero_md(n=200):
        x  = torch.zeros(n, 3)
        ei = torch.zeros(2, 0, dtype=torch.long)
        ea = torch.zeros(0, 3)
        return x, ei, ea

    all_pois_gs = {}
    for poi, sg in static_graphs.items():
        if poi in md_graphs:
            md = md_graphs[poi]
            mx, mei, mea = md['x'], md['ei'], md['ea']
        else:
            mx, mei, mea = zero_md()
        all_pois_gs[poi] = {'poi': poi, 'label': sg['label'],
                            'sx': sg['x'], 'sei': sg['ei'], 'sea': sg['ea'],
                            'mx': mx, 'mei': mei, 'mea': mea}

    # LOPO-CV restricted to 8 MD POIs as test, all 17 as training pool
    pois_test  = POI_ORDER   # 8 MD POIs
    pois_train = list(all_pois_gs.keys())   # 17 total

    print(f"  Training pool: {len(pois_train)} POIs  "
          f"Test targets: {len(pois_test)} MD POIs")

    all_scores = {}
    t0 = time.time()
    for i, test_poi in enumerate(pois_test):
        test_g  = all_pois_gs[test_poi]
        train_gs = [g for p, g in all_pois_gs.items() if p != test_poi]

        scores = []
        for run in range(N_RUNS):
            sc = run_fold(DualEncoder, train_gs, test_g, seed=SEED+run*100)
            scores.append(sc)
        sc = float(np.mean(scores))
        all_scores[test_poi] = sc
        correct = round(sc) == LABELS[test_poi]
        print(f"    {test_poi:18s} y={LABELS[test_poi]}  score={sc:.4f}  "
              f"{'✓' if correct else '✗'}")

    labels = [LABELS[p] for p in pois_test]
    scores = [all_scores[p] for p in pois_test]
    n_ok   = sum(round(s)==l for s,l in zip(scores, labels))
    try:
        auroc = roc_auc_score(labels, scores)
    except Exception:
        auroc = float('nan')
    print(f"    AUROC={auroc:.4f}  Accuracy={n_ok}/8  ({time.time()-t0:.0f}s)")
    return all_scores, auroc, n_ok


# ═══════════════════════════════════════════════════════════════════════════════
# Experiment 2: Per-pose augmentation
# ═══════════════════════════════════════════════════════════════════════════════

def run_e2(per_pose_graphs):
    print("\n" + "="*60)
    print("E2: Per-pose MD augmentation (16-19 graphs/POI)")
    print("="*60)

    n_per_poi = {p: len(gs) for p, gs in per_pose_graphs.items()}
    total = sum(n_per_poi.values())
    print(f"  Total graphs: {total}  |  Per POI: {n_per_poi}")

    if total < 8:
        print("  [SKIP] Per-pose graphs not built yet — run 04_build_per_pose_and_bridge_graphs.py first")
        return {p: float('nan') for p in POI_ORDER}, float('nan'), 0

    def model_fn(): return nn.Sequential(
        _GATWrapper(in_ch=3, edge_dim=3)
    )

    # Use the GATEncoder wrapped in a classifier
    class SingleGAT(nn.Module):
        def __init__(self):
            super().__init__()
            self.enc = GATEncoder(in_ch=3, edge_dim=3)
            self.clf = nn.Sequential(
                nn.Linear(HIDDEN*2, HIDDEN), nn.ReLU(), nn.Dropout(0.3),
                nn.Linear(HIDDEN, 1))
        def forward(self, g):
            emb = self.enc(g['x'], g['ei'], g['ea'])
            return self.clf(emb).squeeze(-1)

    return lopo_cv(SingleGAT, per_pose_graphs, 'E2_PerPose')


# ═══════════════════════════════════════════════════════════════════════════════
# Experiment 3: Bridging-conditional dual-view
# ═══════════════════════════════════════════════════════════════════════════════

class BridgeEncoder(nn.Module):
    """Two GATEncoders: one for bridging contacts, one for non-bridging."""
    def __init__(self):
        super().__init__()
        self.br_enc = GATEncoder(in_ch=3, edge_dim=3)
        self.nb_enc = GATEncoder(in_ch=3, edge_dim=3)
        fused = self.br_enc.out_dim + self.nb_enc.out_dim
        self.clf = nn.Sequential(
            nn.Linear(fused, HIDDEN), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(HIDDEN, 1))

    def forward(self, g):
        br = self.br_enc(g['br_x'], g['br_ei'], g['br_ea'])
        nb = self.nb_enc(g['nb_x'], g['nb_ei'], g['nb_ea'])
        return self.clf(torch.cat([br, nb], dim=1)).squeeze(-1)


def run_e3(bridge_graphs, nonbr_graphs):
    print("\n" + "="*60)
    print("E3: Bridging-conditional dual-view (G_bridge vs G_nonbr)")
    print("="*60)

    if not bridge_graphs or not nonbr_graphs:
        print("  [SKIP] Bridge graphs not built yet — run 04_build_per_pose_and_bridge_graphs.py first")
        return {p: float('nan') for p in POI_ORDER}, float('nan'), 0

    # Build paired graphs
    paired = {}
    for poi in POI_ORDER:
        if poi not in bridge_graphs or poi not in nonbr_graphs:
            print(f"  [WARN] {poi} missing bridge or nonbr graph")
            continue
        bg = bridge_graphs[poi]
        ng = nonbr_graphs[poi]
        paired[poi] = {'poi': poi, 'label': LABELS[poi],
                       'br_x': bg['x'], 'br_ei': bg['ei'], 'br_ea': bg['ea'],
                       'nb_x': ng['x'], 'nb_ei': ng['ei'], 'nb_ea': ng['ea']}

    print(f"  POIs with bridge+nonbr graphs: {list(paired.keys())}")
    return lopo_cv(BridgeEncoder, paired, 'E3_BridgeView')


# ═══════════════════════════════════════════════════════════════════════════════
# Experiment 4: Richer MD node features (aa-type + contact std)
# ═══════════════════════════════════════════════════════════════════════════════

class RichGAT(nn.Module):
    """GATConv with extended node features (base 3 + contact_std 1)."""
    def __init__(self):
        super().__init__()
        self.enc = GATEncoder(in_ch=4, edge_dim=3)
        self.clf = nn.Sequential(
            nn.Linear(HIDDEN*2, HIDDEN), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(HIDDEN, 1))
    def forward(self, g):
        emb = self.enc(g['x_rich'], g['ei'], g['ea'])
        return self.clf(emb).squeeze(-1)


def run_e4(md_graphs, per_pose_graphs):
    print("\n" + "="*60)
    print("E4: Richer MD node features (+ per-pose contact std)")
    print("="*60)

    if not any(len(v) > 0 for v in per_pose_graphs.values()):
        print("  [SKIP] Per-pose graphs needed for std computation")
        return {p: float('nan') for p in POI_ORDER}, float('nan'), 0

    node_std = compute_per_pose_node_std(per_pose_graphs)

    rich_graphs = {}
    for poi in POI_ORDER:
        if poi not in md_graphs:
            continue
        g = md_graphs[poi]
        std_feat = node_std.get(poi)
        if std_feat is None:
            std_feat = torch.zeros(g['x'].shape[0])
        x_rich = torch.cat([g['x'], std_feat.unsqueeze(1)], dim=1)
        rich_graphs[poi] = {'poi': poi, 'label': LABELS[poi],
                            'x_rich': x_rich,
                            'ei': g['ei'], 'ea': g['ea']}

    return lopo_cv(RichGAT, rich_graphs, 'E4_RichFeatures')


# ═══════════════════════════════════════════════════════════════════════════════
# Experiment 5: Combined (dual encoder + bridge-view + rich features)
# ═══════════════════════════════════════════════════════════════════════════════

class CombinedEncoder(nn.Module):
    def __init__(self):
        super().__init__()
        self.static_enc = GATEncoder(in_ch=9, edge_dim=3)
        self.md_enc     = GATEncoder(in_ch=4, edge_dim=3)
        self.br_enc     = GATEncoder(in_ch=3, edge_dim=3)
        self.nb_enc     = GATEncoder(in_ch=3, edge_dim=3)
        fused = (self.static_enc.out_dim + self.md_enc.out_dim +
                 self.br_enc.out_dim + self.nb_enc.out_dim)
        self.clf = nn.Sequential(
            nn.Linear(fused, HIDDEN*2), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(HIDDEN*2, HIDDEN), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(HIDDEN, 1))

    def forward(self, g):
        s_emb  = self.static_enc(g['sx'],    g['sei'],   g['sea'])
        md_emb = self.md_enc(g['x_rich'],    g['ei'],    g['ea'])
        br_emb = self.br_enc(g['br_x'],      g['br_ei'], g['br_ea'])
        nb_emb = self.nb_enc(g['nb_x'],      g['nb_ei'], g['nb_ea'])
        return self.clf(torch.cat([s_emb, md_emb, br_emb, nb_emb], dim=1)).squeeze(-1)


def run_e5(static_graphs, md_graphs, bridge_graphs, nonbr_graphs,
           per_pose_graphs, node_std):
    print("\n" + "="*60)
    print("E5: Combined (dual-static + bridge-view + rich-features + 17-POI)")
    print("="*60)

    if not bridge_graphs:
        print("  [SKIP] Bridge graphs needed")
        return {p: float('nan') for p in POI_ORDER}, float('nan'), 0

    def zero_g3(n=200): # zero MD graph (3-dim)
        return torch.zeros(n,3), torch.zeros(2,0,dtype=torch.long), torch.zeros(0,3)
    def zero_g4(n=200):
        return torch.zeros(n,4), torch.zeros(2,0,dtype=torch.long), torch.zeros(0,3)

    all_pois_gs = {}
    for poi, sg in static_graphs.items():
        if poi not in LABELS and poi not in STATIC_ONLY:
            continue
        lab = sg['label']
        md = md_graphs.get(poi)
        br = bridge_graphs.get(poi)
        nb = nonbr_graphs.get(poi)

        if md:
            std_feat = node_std.get(poi, torch.zeros(md['x'].shape[0]))
            x_rich = torch.cat([md['x'], std_feat.unsqueeze(1)], dim=1)
            mx, mei, mea = x_rich, md['ei'], md['ea']
        else:
            mx, mei, mea = zero_g4()

        if br:
            brx, brei, brea = br['x'], br['ei'], br['ea']
        else:
            brx, brei, brea = zero_g3()

        if nb:
            nbx, nbei, nbea = nb['x'], nb['ei'], nb['ea']
        else:
            nbx, nbei, nbea = zero_g3()

        all_pois_gs[poi] = {'poi': poi, 'label': lab,
                            'sx': sg['x'], 'sei': sg['ei'], 'sea': sg['ea'],
                            'x_rich': mx, 'ei': mei, 'ea': mea,
                            'br_x': brx, 'br_ei': brei, 'br_ea': brea,
                            'nb_x': nbx, 'nb_ei': nbei, 'nb_ea': nbea}

    # LOPO-CV: test only on MD POIs, train on all
    pois_test = POI_ORDER
    t0 = time.time()
    all_scores = {}
    for test_poi in pois_test:
        if test_poi not in all_pois_gs:
            continue
        test_g  = all_pois_gs[test_poi]
        train_gs = [g for p,g in all_pois_gs.items() if p != test_poi]
        scores = [run_fold(CombinedEncoder, train_gs, test_g, seed=SEED+r*100)
                  for r in range(N_RUNS)]
        sc = float(np.mean(scores))
        all_scores[test_poi] = sc
        correct = round(sc) == LABELS[test_poi]
        print(f"    {test_poi:18s} y={LABELS[test_poi]}  score={sc:.4f}  "
              f"{'✓' if correct else '✗'}")

    labels = [LABELS[p] for p in pois_test if p in all_scores]
    scores = [all_scores[p] for p in pois_test if p in all_scores]
    n_ok   = sum(round(s)==l for s,l in zip(scores, labels))
    try:
        auroc = roc_auc_score(labels, scores)
    except Exception:
        auroc = float('nan')
    print(f"    AUROC={auroc:.4f}  Accuracy={n_ok}/{len(scores)}  ({time.time()-t0:.0f}s)")
    return all_scores, auroc, n_ok


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("="*70)
    print("GNN EXPERIMENTS — PROTAC DEGRADER vs NON-DEGRADER")
    print(f"Device: {device}  Epochs: {EPOCHS}  Runs/fold: {N_RUNS}")
    print("="*70)

    md_graphs       = {g['poi']: g for g in load_md_graphs()}
    static_graphs   = load_static_graphs()
    per_pose_graphs = load_per_pose_graphs()
    bridge_graphs, nonbr_graphs = load_bridge_graphs()

    n_pp = sum(len(v) for v in per_pose_graphs.values())
    print(f"\nLoaded: {len(md_graphs)} MD graphs | {len(static_graphs)} static graphs | "
          f"{n_pp} per-pose graphs | {len(bridge_graphs)} bridge | {len(nonbr_graphs)} nonbr")

    node_std = compute_per_pose_node_std(per_pose_graphs) if n_pp > 0 else {}

    summary = []

    # E1
    scores1, auroc1, nok1 = run_e1(static_graphs, md_graphs)
    summary.append({'exp':'E1_DualEncoder_17POI','auroc':auroc1,'n_correct':nok1,'n_total':8})

    # E2
    if n_pp > 0:
        scores2, auroc2, nok2 = run_e2(per_pose_graphs)
    else:
        print("\n  [E2 SKIP] Run graph builder first")
        scores2 = {p: float('nan') for p in POI_ORDER}; auroc2 = float('nan'); nok2 = 0
    summary.append({'exp':'E2_PerPoseAugment','auroc':auroc2,'n_correct':nok2,'n_total':8})

    # E3
    scores3, auroc3, nok3 = run_e3(bridge_graphs, nonbr_graphs)
    summary.append({'exp':'E3_BridgeView','auroc':auroc3,'n_correct':nok3,'n_total':8})

    # E4
    scores4, auroc4, nok4 = run_e4(md_graphs, per_pose_graphs)
    summary.append({'exp':'E4_RichFeatures','auroc':auroc4,'n_correct':nok4,'n_total':8})

    # E5
    scores5, auroc5, nok5 = run_e5(static_graphs, md_graphs, bridge_graphs,
                                     nonbr_graphs, per_pose_graphs, node_std)
    summary.append({'exp':'E5_Combined','auroc':auroc5,'n_correct':nok5,'n_total':8})

    # ── Save ──
    rows = []
    for poi in POI_ORDER:
        rows.append({'poi': poi, 'label': LABELS[poi],
                     'E1_DualEncoder': scores1.get(poi, float('nan')),
                     'E2_PerPose':     scores2.get(poi, float('nan')),
                     'E3_BridgeView':  scores3.get(poi, float('nan')),
                     'E4_RichFeat':    scores4.get(poi, float('nan')),
                     'E5_Combined':    scores5.get(poi, float('nan'))})
    df_comp = pd.DataFrame(rows).set_index('poi')
    df_summ = pd.DataFrame(summary)

    df_comp.to_csv(RESULTS_DIR / 'all_experiments.tsv', sep='\t', float_format='%.4f')
    df_summ.to_csv(RESULTS_DIR / 'all_summary.tsv', sep='\t', index=False, float_format='%.4f')

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"\n{'Experiment':<30} {'AUROC':>8} {'Correct':>10}")
    print('-'*52)
    for r in summary:
        flag = "  ← BEST" if r['auroc'] == max(x['auroc'] for x in summary
                                                if not np.isnan(x['auroc'])) else ""
        print(f"{r['exp']:<30} {r['auroc']:>8.4f} {r['n_correct']:>5}/{r['n_total']}{flag}")
    print(f"\nReference (geometric models): AUROC 0.00–0.25, 6/8 correct")
    print(f"Reference (E2-GNN, 17-POI):   AUROC ~0.85, 7/8 correct")
    print(f"\nResults: {RESULTS_DIR}")


if __name__ == '__main__':
    main()
