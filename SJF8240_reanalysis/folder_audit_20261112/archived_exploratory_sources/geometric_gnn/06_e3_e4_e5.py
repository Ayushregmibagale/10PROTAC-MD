#!/usr/bin/env python3
"""
Run E3, E4, E5 only (E1 already done, E2 skipped — too slow).

E1 result baked in from logs/experiments.log:
  AUROC=0.1667  Accuracy=6/8

E3  Bridging-conditional dual-view (G_bridge + G_nonbr)
E4  Richer MD node features (+ per-pose contact std)
E5  Combined (dual-static + bridge-view + rich-features + 17-POI)

Run: ~/.conda/envs/embeddings/bin/python3 -u 06_e3_e4_e5.py
"""
import sys, time, random, warnings
import numpy as np
import pandas as pd
from pathlib import Path

warnings.filterwarnings('ignore')

try:
    import torch, torch.nn as nn, torch.nn.functional as F
    from torch_geometric.nn import GATConv, global_mean_pool, global_max_pool
    from sklearn.metrics import roc_auc_score
except ImportError as e:
    print(f"[ERROR] {e}\nRun with ~/.conda/envs/embeddings/bin/python3")
    sys.exit(1)

SCRIPT_DIR     = Path(__file__).resolve().parent
GRAPH_DIR_ORIG = SCRIPT_DIR.parent / 'corrected_md_only_v1/inputs/elob'
STATIC_PT      = Path.home() / 'protac-md/analysis/gnn/ensemble_graphs.pt'
BRIDGE_DIR     = SCRIPT_DIR / 'graphs/bridge'
PP_DIR         = SCRIPT_DIR / 'graphs/per_pose'
RESULTS_DIR    = SCRIPT_DIR / 'results'
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

LABELS = {
    'DDR2':1,'ABL1':0,'AXL':0,'EPHA2':1,
    'EPHB2':1,'MAPK14_full':1,'MET':1,'RIPK2_full':1,
}
POI_ORDER = ['DDR2','ABL1','AXL','EPHA2','EPHB2','MAPK14_full','MET','RIPK2_full']

STATIC_ONLY = {
    'MAP4K5':0,'SLK':0,'CDK2':0,'CDK4':0,'CDK6':0,'DAPK1':0,'IRAK1':0,'TLK2':0
}

EPOCHS  = 500
LR      = 1e-3
HIDDEN  = 64
N_RUNS  = 5
SEED    = 42
device  = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# E1 results (from logs/experiments.log — 2503s run)
E1_SCORES = {
    'DDR2':0.9967,'ABL1':0.9995,'AXL':1.0000,'EPHA2':0.9988,
    'EPHB2':0.9994,'MAPK14_full':0.9987,'MET':1.0000,'RIPK2_full':0.9991,
}
E1_AUROC = 0.1667
E1_NOK   = 6


# ── Data loaders ─────────────────────────────────────────────────────────────

def load_md_graphs():
    graphs = {}
    for poi in POI_ORDER:
        p = GRAPH_DIR_ORIG / f'{poi}.npz'
        if not p.exists():
            print(f"  [WARN] MD graph missing: {p}")
            continue
        d = np.load(p)
        graphs[poi] = {
            'poi': poi, 'label': LABELS[poi],
            'x':  torch.tensor(d['x'],          dtype=torch.float32),
            'ei': torch.tensor(d['edge_index'],  dtype=torch.long),
            'ea': torch.tensor(d['edge_attr'],   dtype=torch.float32),
        }
    return graphs


def load_static_graphs():
    static = torch.load(STATIC_PT, weights_only=False, map_location='cpu')
    keep = set(list(LABELS.keys()) + list(STATIC_ONLY.keys()))
    out = {}
    for g in static:
        if g.poi not in keep:
            continue
        lab = LABELS.get(g.poi, STATIC_ONLY.get(g.poi, int(g.y)))
        out[g.poi] = {'poi': g.poi, 'label': lab,
                      'x': g.x, 'ei': g.edge_index, 'ea': g.edge_attr}
    return out


def load_per_pose_graphs():
    out = {poi: [] for poi in POI_ORDER}
    if not PP_DIR.exists():
        return out
    for f in sorted(PP_DIR.glob('*.npz')):
        parts = f.stem.split('_pose')
        if len(parts) != 2:
            continue
        poi = parts[0]
        if poi not in LABELS:
            continue
        d = np.load(f)
        out[poi].append({
            'poi': poi, 'label': LABELS[poi],
            'x':  torch.tensor(d['x'],          dtype=torch.float32),
            'ei': torch.tensor(d['edge_index'],  dtype=torch.long),
            'ea': torch.tensor(d['edge_attr'],   dtype=torch.float32),
        })
    return out


def load_bridge_graphs():
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
    out = {}
    for poi, poses in per_pose_graphs.items():
        if not poses:
            continue
        n = poses[0]['x'].shape[0]
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
        out[poi] = torch.stack(stacked, 0).std(0)
    return out


# ── GNN Models ───────────────────────────────────────────────────────────────

class GATEncoder(nn.Module):
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


def to_dev(g):
    return {k: v.to(device) if isinstance(v, torch.Tensor) else v
            for k, v in g.items()}


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
    pois       = list(graphs_by_poi.keys())
    all_labels = [LABELS[p] for p in pois]
    all_scores = []

    print(f"\n  [{name}] {N_RUNS} runs/fold ...")
    t0 = time.time()
    for test_poi in pois:
        test_g   = graphs_by_poi[test_poi]
        train_gs = [graphs_by_poi[p] for p in pois if p != test_poi]

        scores = [run_fold(model_fn, train_gs, test_g,
                           seed=SEED+r*100, epochs=epochs)
                  for r in range(N_RUNS)]
        sc = float(np.mean(scores))
        all_scores.append(sc)
        correct = round(sc) == LABELS[test_poi]
        print(f"    {test_poi:18s} y={LABELS[test_poi]}  score={sc:.4f}  "
              f"{'✓' if correct else '✗'}", flush=True)

    n_ok = sum(round(s)==l for s,l in zip(all_scores, all_labels))
    try:
        auroc = roc_auc_score(all_labels, all_scores)
    except Exception:
        auroc = float('nan')
    print(f"    AUROC={auroc:.4f}  Accuracy={n_ok}/{len(pois)}  ({time.time()-t0:.0f}s)",
          flush=True)
    return dict(zip(pois, all_scores)), auroc, n_ok


# ── E3: Bridging-conditional dual-view ───────────────────────────────────────

class BridgeEncoder(nn.Module):
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
        print("  [SKIP] Bridge graphs missing")
        return {p: float('nan') for p in POI_ORDER}, float('nan'), 0

    paired = {}
    for poi in POI_ORDER:
        if poi not in bridge_graphs or poi not in nonbr_graphs:
            print(f"  [WARN] {poi} missing bridge or nonbr")
            continue
        bg, ng = bridge_graphs[poi], nonbr_graphs[poi]
        paired[poi] = {
            'poi': poi, 'label': LABELS[poi],
            'br_x': bg['x'], 'br_ei': bg['ei'], 'br_ea': bg['ea'],
            'nb_x': ng['x'], 'nb_ei': ng['ei'], 'nb_ea': ng['ea'],
        }

    print(f"  POIs: {list(paired.keys())}")
    return lopo_cv(BridgeEncoder, paired, 'E3_BridgeView')


# ── E4: Richer MD node features ───────────────────────────────────────────────

class RichGAT(nn.Module):
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

    n_pp = sum(len(v) for v in per_pose_graphs.values())
    if n_pp == 0:
        print("  [SKIP] Per-pose graphs needed for std computation")
        return {p: float('nan') for p in POI_ORDER}, float('nan'), 0

    node_std = compute_per_pose_node_std(per_pose_graphs)
    rich_graphs = {}
    for poi in POI_ORDER:
        if poi not in md_graphs:
            continue
        g = md_graphs[poi]
        std_feat = node_std.get(poi, torch.zeros(g['x'].shape[0]))
        x_rich = torch.cat([g['x'], std_feat.unsqueeze(1)], dim=1)
        rich_graphs[poi] = {
            'poi': poi, 'label': LABELS[poi],
            'x_rich': x_rich, 'ei': g['ei'], 'ea': g['ea'],
        }

    return lopo_cv(RichGAT, rich_graphs, 'E4_RichFeatures')


# ── E5: Combined ─────────────────────────────────────────────────────────────

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
        s_emb  = self.static_enc(g['sx'],     g['sei'],   g['sea'])
        md_emb = self.md_enc(g['x_rich'],     g['ei'],    g['ea'])
        br_emb = self.br_enc(g['br_x'],       g['br_ei'], g['br_ea'])
        nb_emb = self.nb_enc(g['nb_x'],       g['nb_ei'], g['nb_ea'])
        return self.clf(torch.cat([s_emb, md_emb, br_emb, nb_emb], dim=1)).squeeze(-1)


def run_e5(static_graphs, md_graphs, bridge_graphs, nonbr_graphs, node_std):
    print("\n" + "="*60)
    print("E5: Combined (dual-static + bridge-view + rich-features + 17-POI)")
    print("="*60)

    if not bridge_graphs:
        print("  [SKIP] Bridge graphs needed")
        return {p: float('nan') for p in POI_ORDER}, float('nan'), 0

    def zero3(n=200):
        return torch.zeros(n,3), torch.zeros(2,0,dtype=torch.long), torch.zeros(0,3)
    def zero4(n=200):
        return torch.zeros(n,4), torch.zeros(2,0,dtype=torch.long), torch.zeros(0,3)

    all_pois_gs = {}
    for poi, sg in static_graphs.items():
        lab = sg['label']
        md = md_graphs.get(poi)
        br = bridge_graphs.get(poi)
        nb = nonbr_graphs.get(poi)

        if md:
            sf = node_std.get(poi, torch.zeros(md['x'].shape[0]))
            x_rich = torch.cat([md['x'], sf.unsqueeze(1)], dim=1)
            mx, mei, mea = x_rich, md['ei'], md['ea']
        else:
            mx, mei, mea = zero4()

        brx, brei, brea = (br['x'], br['ei'], br['ea']) if br else zero3()
        nbx, nbei, nbea = (nb['x'], nb['ei'], nb['ea']) if nb else zero3()

        all_pois_gs[poi] = {
            'poi': poi, 'label': lab,
            'sx': sg['x'], 'sei': sg['ei'], 'sea': sg['ea'],
            'x_rich': mx, 'ei': mei, 'ea': mea,
            'br_x': brx, 'br_ei': brei, 'br_ea': brea,
            'nb_x': nbx, 'nb_ei': nbei, 'nb_ea': nbea,
        }

    t0 = time.time()
    all_scores = {}
    print(f"  [{N_RUNS} runs/fold ...]")
    for test_poi in POI_ORDER:
        if test_poi not in all_pois_gs:
            continue
        test_g   = all_pois_gs[test_poi]
        train_gs = [g for p,g in all_pois_gs.items() if p != test_poi]
        scores = [run_fold(CombinedEncoder, train_gs, test_g, seed=SEED+r*100)
                  for r in range(N_RUNS)]
        sc = float(np.mean(scores))
        all_scores[test_poi] = sc
        correct = round(sc) == LABELS[test_poi]
        print(f"    {test_poi:18s} y={LABELS[test_poi]}  score={sc:.4f}  "
              f"{'✓' if correct else '✗'}", flush=True)

    labels = [LABELS[p] for p in POI_ORDER if p in all_scores]
    scores = [all_scores[p] for p in POI_ORDER if p in all_scores]
    n_ok = sum(round(s)==l for s,l in zip(scores, labels))
    try:
        auroc = roc_auc_score(labels, scores)
    except Exception:
        auroc = float('nan')
    print(f"    AUROC={auroc:.4f}  Accuracy={n_ok}/{len(scores)}  ({time.time()-t0:.0f}s)",
          flush=True)
    return all_scores, auroc, n_ok


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    print("="*70)
    print("GNN EXPERIMENTS E3/E4/E5 — PROTAC DEGRADER vs NON-DEGRADER")
    print(f"Device: {device}  Epochs: {EPOCHS}  Runs/fold: {N_RUNS}")
    print("="*70)

    md_graphs          = load_md_graphs()
    static_graphs      = load_static_graphs()
    per_pose_graphs    = load_per_pose_graphs()
    bridge_graphs, nonbr_graphs = load_bridge_graphs()
    n_pp = sum(len(v) for v in per_pose_graphs.values())
    node_std = compute_per_pose_node_std(per_pose_graphs) if n_pp > 0 else {}

    print(f"\nLoaded: {len(md_graphs)} MD | {len(static_graphs)} static | "
          f"{n_pp} per-pose | {len(bridge_graphs)} bridge | {len(nonbr_graphs)} nonbr")

    summary = [
        {'exp':'E1_DualEncoder_17POI', 'auroc':E1_AUROC, 'n_correct':E1_NOK, 'n_total':8,
         'note':'from prior run'},
        {'exp':'E2_PerPoseAugment',    'auroc':float('nan'), 'n_correct':0, 'n_total':8,
         'note':'skipped (>100h runtime)'},
    ]

    scores3, auroc3, nok3 = run_e3(bridge_graphs, nonbr_graphs)
    summary.append({'exp':'E3_BridgeView',   'auroc':auroc3, 'n_correct':nok3, 'n_total':8, 'note':''})

    scores4, auroc4, nok4 = run_e4(md_graphs, per_pose_graphs)
    summary.append({'exp':'E4_RichFeatures', 'auroc':auroc4, 'n_correct':nok4, 'n_total':8, 'note':''})

    scores5, auroc5, nok5 = run_e5(static_graphs, md_graphs, bridge_graphs, nonbr_graphs, node_std)
    summary.append({'exp':'E5_Combined',     'auroc':auroc5, 'n_correct':nok5, 'n_total':8, 'note':''})

    # ── Save ──
    rows = []
    for poi in POI_ORDER:
        rows.append({
            'poi':   poi,
            'label': LABELS[poi],
            'E1_DualEncoder': E1_SCORES.get(poi, float('nan')),
            'E2_PerPose':     float('nan'),
            'E3_BridgeView':  scores3.get(poi, float('nan')),
            'E4_RichFeat':    scores4.get(poi, float('nan')),
            'E5_Combined':    scores5.get(poi, float('nan')),
        })
    df_comp = pd.DataFrame(rows).set_index('poi')
    df_summ = pd.DataFrame(summary)

    df_comp.to_csv(RESULTS_DIR / 'all_experiments.tsv', sep='\t', float_format='%.4f')
    df_summ.to_csv(RESULTS_DIR / 'all_summary.tsv',   sep='\t', index=False, float_format='%.4f')

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"\n{'Experiment':<30} {'AUROC':>8} {'Correct':>10}")
    print('-'*52)
    valid_aurocs = [r['auroc'] for r in summary if not np.isnan(r['auroc'])]
    best_auroc   = max(valid_aurocs) if valid_aurocs else float('nan')
    for r in summary:
        flag = "  ← BEST" if r['auroc'] == best_auroc else ""
        astr = f"{r['auroc']:8.4f}" if not np.isnan(r['auroc']) else "     n/a"
        print(f"{r['exp']:<30} {astr} {r['n_correct']:>5}/{r['n_total']}{flag}")
    print(f"\nReference (geometric models): AUROC 0.00–0.25, 6/8 correct")
    print(f"Reference (E2-GNN, 17-POI):   AUROC ~0.85, 7/8 correct")
    print(f"\nResults saved: {RESULTS_DIR}")


if __name__ == '__main__':
    main()
