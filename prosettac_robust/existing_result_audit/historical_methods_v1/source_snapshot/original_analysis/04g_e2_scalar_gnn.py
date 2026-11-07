#!/usr/bin/env python3
"""
04g — Dual-encoder GNN + E2 reachability scalars injected at fusion stage.

Identical to 04d (dual-encoder) except the fusion MLP receives additional
graph-level E2 scalar features:
  - e2_mean_min_dist_A      : mean E2–lysine distance during bridging (Å)
  - e2_delta_reach_50A      : Δ E2 reach at 50 Å (bridging vs non-bridging)
  - e2_delta_reach_80A      : Δ E2 reach at 80 Å
  - lk_lys_nearest_e2_dist_A: nearest E2–lysine distance across all frames

These are NaN for non-MD POIs (no bridging → no E2 estimate).
Imputation: per-fold z-normalise on training set; NaN → 0.0 after normalisation
(so non-MD POIs receive a zero-mean signal, not a spurious discriminating value).

Hypothesis: ABL1/AXL/EPHA2 may be separable from DEGs via E2 accessibility
patterns not captured by graph topology alone.

Usage:
  python3 scripts/02_analysis/04g_e2_scalar_gnn.py
  python3 scripts/02_analysis/04g_e2_scalar_gnn.py --epochs 400 --hidden 64
"""
import argparse
import sys
import numpy as np
import pandas as pd
from pathlib import Path

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    from torch_geometric.data import Data
    from torch_geometric.nn import GATConv, global_mean_pool, global_max_pool
    from sklearn.metrics import roc_auc_score, accuracy_score
    HAS_DEPS = True
except ImportError as e:
    HAS_DEPS = False
    MISSING = str(e)

GNN_DIR    = Path.home() / 'protac-md' / 'analysis' / 'gnn'
STATIC_PT  = GNN_DIR / 'ensemble_graphs.pt'
DYNAMIC_PT = GNN_DIR / 'md_contact_graphs.pt'
ML_TABLE   = Path.home() / 'protac-md' / 'analysis' / 'features' / 'ml_table.tsv'
OUT_PRED   = GNN_DIR / 'e2gnn_lopo_results.tsv'
OUT_SUMM   = GNN_DIR / 'e2gnn_lopo_summary.tsv'

E2_COLS = [
    'e2_mean_min_dist_A',
    'e2_delta_reach_50A',
    'e2_delta_reach_80A',
    'lk_lys_nearest_e2_dist_A',
]

POI_17 = {
    'DDR2', 'RIPK2_full', 'MET', 'MAPK14_full', 'EPHB2',
    'EPHA2', 'ABL1', 'AXL', 'TNIK', 'MAP4K5', 'SLK',
    'CDK2', 'CDK4', 'CDK6', 'DAPK1', 'IRAK1', 'TLK2',
}
DEGRADERS = {'DDR2', 'RIPK2_full', 'MET', 'MAPK14_full', 'EPHB2'}


# ── Encoders (identical to 04d) ───────────────────────────────────────────────

class EnsembleEncoder(nn.Module):
    def __init__(self, in_ch=9, edge_dim=3, hidden=64, heads=4, dropout=0.3):
        super().__init__()
        self.edge_proj = nn.Linear(edge_dim, hidden)
        self.conv1 = GATConv(in_ch, hidden // heads, heads=heads, edge_dim=hidden,
                             dropout=dropout, concat=True)
        self.bn1 = nn.BatchNorm1d(hidden)
        self.conv2 = GATConv(hidden, hidden // heads, heads=heads, edge_dim=hidden,
                             dropout=dropout, concat=True)
        self.bn2 = nn.BatchNorm1d(hidden)
        self.conv3 = GATConv(hidden, hidden, heads=1, edge_dim=hidden,
                             dropout=dropout, concat=False)
        self.bn3 = nn.BatchNorm1d(hidden)
        self.dropout = nn.Dropout(dropout)
        self.out_dim = hidden * 2

    def forward(self, data):
        x, ei, ea = data.x, data.edge_index, data.edge_attr
        batch = getattr(data, 'batch', None)
        if batch is None:
            batch = torch.zeros(x.size(0), dtype=torch.long, device=x.device)
        ea_emb = F.relu(self.edge_proj(ea))
        x = F.relu(self.bn1(self.conv1(x, ei, ea_emb)))
        x = self.dropout(x)
        x = F.relu(self.bn2(self.conv2(x, ei, ea_emb)))
        x = self.dropout(x)
        x = F.relu(self.bn3(self.conv3(x, ei, ea_emb)))
        return torch.cat([global_mean_pool(x, batch),
                          global_max_pool(x, batch)], dim=1)


class MDEncoder(nn.Module):
    def __init__(self, in_ch=3, edge_dim=3, hidden=64, heads=4, dropout=0.3):
        super().__init__()
        self.edge_proj = nn.Linear(edge_dim, hidden)
        self.conv1 = GATConv(in_ch, hidden // heads, heads=heads, edge_dim=hidden,
                             dropout=dropout, concat=True)
        self.bn1 = nn.BatchNorm1d(hidden)
        self.conv2 = GATConv(hidden, hidden, heads=1, edge_dim=hidden,
                             dropout=dropout, concat=False)
        self.bn2 = nn.BatchNorm1d(hidden)
        self.dropout = nn.Dropout(dropout)
        self.out_dim = hidden * 2

    def forward(self, data):
        x, ei, ea = data.x, data.edge_index, data.edge_attr
        batch = getattr(data, 'batch', None)
        if batch is None:
            batch = torch.zeros(x.size(0), dtype=torch.long, device=x.device)
        ea_emb = F.relu(self.edge_proj(ea))
        x = F.relu(self.bn1(self.conv1(x, ei, ea_emb)))
        x = self.dropout(x)
        x = F.relu(self.bn2(self.conv2(x, ei, ea_emb)))
        return torch.cat([global_mean_pool(x, batch),
                          global_max_pool(x, batch)], dim=1)


# ── E2-augmented dual encoder ─────────────────────────────────────────────────

class E2DualEncoderGNN(nn.Module):
    """
    Dual-encoder GNN with E2 reachability scalars at fusion.

    fused = concat(emb_A, gate * emb_B, e2_vec)
    where e2_vec is z-normalised per LOPO fold (NaN → 0.0 after normalisation).
    """
    def __init__(self, static_in=9, static_edge=3, dynamic_in=3, dynamic_edge=3,
                 n_e2_scalars=4, hidden=64, heads=4, dropout=0.3):
        super().__init__()
        self.enc_a = EnsembleEncoder(static_in, static_edge, hidden, heads, dropout)
        self.enc_b = MDEncoder(dynamic_in, dynamic_edge, hidden, heads, dropout)
        self.md_gate = nn.Parameter(torch.tensor(0.5))

        fused_dim = self.enc_a.out_dim + self.enc_b.out_dim + n_e2_scalars
        self.fusion = nn.Sequential(
            nn.Linear(fused_dim, hidden * 2),
            nn.LayerNorm(hidden * 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden * 2, hidden),
            nn.ReLU(),
            nn.Linear(hidden, 1),
        )

    def forward(self, static_graph, dynamic_graph, e2_vec):
        emb_a = self.enc_a(static_graph)

        if dynamic_graph is not None:
            emb_b = self.enc_b(dynamic_graph)
            gate  = torch.sigmoid(self.md_gate)
            emb_b = gate * emb_b
        else:
            emb_b = torch.zeros(1, self.enc_b.out_dim, device=emb_a.device)

        fused = torch.cat([emb_a, emb_b, e2_vec], dim=1)
        return self.fusion(fused).squeeze(1)


# ── E2 scalar loading and normalisation ───────────────────────────────────────

def load_e2_table(ml_table_path):
    df = pd.read_csv(ml_table_path, sep='\t')
    df = df[df['poi'].isin(POI_17) & (df['prostac_artifact'] == 0)].copy()
    df = df.set_index('poi')
    return df[E2_COLS]


def get_e2_vec(poi, e2_df, train_pois, device):
    """
    Return z-normalised E2 scalar tensor for `poi`.
    Normalisation stats computed from training POIs only.
    NaN (non-MD POIs, or missing features) → 0.0 after normalisation.
    """
    train_vals = e2_df.loc[e2_df.index.isin(train_pois), E2_COLS]
    means = train_vals.mean()
    stds  = train_vals.std().replace(0, 1.0)

    if poi in e2_df.index:
        raw = e2_df.loc[poi, E2_COLS].values.astype(float)
    else:
        raw = np.full(len(E2_COLS), np.nan)

    normed = (raw - means.values) / stds.values
    normed = np.nan_to_num(normed, nan=0.0)
    return torch.tensor(normed, dtype=torch.float32, device=device).unsqueeze(0)


# ── Training loop ─────────────────────────────────────────────────────────────

def train_step(model, opt, static_g, dynamic_g, e2_vec, label, device):
    model.train()
    opt.zero_grad()
    static_g = static_g.to(device)
    dynamic_g = dynamic_g.to(device) if dynamic_g is not None else None
    out = model(static_g, dynamic_g, e2_vec)
    tgt = torch.tensor([label], dtype=torch.float32, device=device)
    loss = F.binary_cross_entropy_with_logits(out, tgt)
    loss.backward()
    opt.step()
    return float(loss)


@torch.no_grad()
def predict(model, static_g, dynamic_g, e2_vec, device):
    model.eval()
    static_g = static_g.to(device)
    dynamic_g = dynamic_g.to(device) if dynamic_g is not None else None
    out = model(static_g, dynamic_g, e2_vec)
    return float(torch.sigmoid(out).cpu())


# ── LOPO-CV ───────────────────────────────────────────────────────────────────

def run_lopo(static_graphs, dynamic_map, e2_df, args):
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Device: {device}")
    print(f"E2 scalars: {E2_COLS}\n")

    clean    = [g for g in static_graphs if not g.artifact]
    md_pois  = set(dynamic_map.keys())
    n_e2     = len(E2_COLS)

    print(f"Clean POIs: {len(clean)}")
    print(f"POIs with MD:  {sorted(md_pois)}")
    print(f"POIs with E2:  {sorted(e2_df.dropna(how='all').index.tolist())}\n")

    pred_rows = []

    for fold_i, test_g in enumerate(clean):
        poi        = test_g.poi
        test_dyn   = dynamic_map.get(poi, None)
        train_pois = [g.poi for g in clean if g.poi != poi]

        model = E2DualEncoderGNN(
            n_e2_scalars=n_e2,
            hidden=args.hidden, heads=args.heads, dropout=args.dropout,
        ).to(device)
        optimizer = torch.optim.Adam(
            model.parameters(), lr=args.lr, weight_decay=1e-4)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
            optimizer, args.epochs)

        train_pairs = [(g, dynamic_map.get(g.poi, None)) for g in clean
                       if g.poi != poi]

        for epoch in range(args.epochs):
            np.random.shuffle(train_pairs)
            for sg, dg in train_pairs:
                e2v = get_e2_vec(sg.poi, e2_df, train_pois, device)
                train_step(model, optimizer, sg, dg, e2v,
                           float(sg.y.item()), device)
            scheduler.step()
            if args.verbose and epoch % 100 == 0:
                print(f"    fold {fold_i+1}/{len(clean)} epoch {epoch}")

        e2_test = get_e2_vec(poi, e2_df, train_pois, device)
        prob    = predict(model, test_g, test_dyn, e2_test, device)
        label   = int(test_g.y.item())
        has_md  = poi in md_pois
        has_e2  = (poi in e2_df.index and
                   not e2_df.loc[poi, E2_COLS].isna().all())

        # Show raw E2 values for this POI
        if has_e2:
            raw_e2 = e2_df.loc[poi, 'e2_mean_min_dist_A']
            e2_str = f"e2_dist={raw_e2:.1f}Å"
        else:
            e2_str = "e2=NaN"

        print(f"  [{poi:14}] label={label}  prob={prob:.4f}  "
              f"has_md={has_md}  {e2_str}  "
              f"pred={'DEG' if prob>=0.5 else 'NON'}  "
              f"{'✓' if (prob>=0.5)==bool(label) else '✗'}", flush=True)

        pred_rows.append({
            'fold': fold_i + 1, 'poi': poi,
            'label': label, 'prob_degrader': round(prob, 4),
            'pred': int(prob >= 0.5),
            'correct': int((prob >= 0.5) == bool(label)),
            'has_md': int(has_md),
            'has_e2': int(has_e2),
            'e2_mean_min_dist_A': (
                e2_df.loc[poi, 'e2_mean_min_dist_A']
                if poi in e2_df.index else float('nan')
            ),
        })

    pred_df = pd.DataFrame(pred_rows)
    labels  = pred_df['label'].values
    probs   = pred_df['prob_degrader'].values
    preds   = pred_df['pred'].values

    try:
        auroc = roc_auc_score(labels, probs)
    except Exception:
        auroc = float('nan')

    acc       = accuracy_score(labels, preds)
    n_correct = int(pred_df['correct'].sum())
    n_total   = len(pred_df)

    md_mask  = pred_df['has_md'] == 1
    nmd_mask = pred_df['has_md'] == 0
    try:
        auroc_md  = roc_auc_score(labels[md_mask],  probs[md_mask])
        auroc_nmd = roc_auc_score(labels[nmd_mask], probs[nmd_mask])
    except Exception:
        auroc_md = auroc_nmd = float('nan')

    print(f"\n{'='*60}")
    print(f"E2-SCALAR GNN LOPO-CV RESULTS")
    print(f"{'='*60}")
    print(f"  Overall:    accuracy={n_correct}/{n_total}={acc:.3f}  AUROC={auroc:.4f}")
    print(f"  MD POIs:    AUROC={auroc_md:.4f}  (n={int(md_mask.sum())})")
    print(f"  Non-MD POIs: AUROC={auroc_nmd:.4f}  (n={int(nmd_mask.sum())})")
    print(f"\n  Reference (04d dual-encoder, no E2 scalars): AUROC=0.850 ACC=14/17")
    print(f"\n  Per-POI summary (MD POIs only):")
    for _, r in pred_df[pred_df['has_md']==1].iterrows():
        e2v = r['e2_mean_min_dist_A']
        e2s = f"{e2v:.1f}Å" if not pd.isna(e2v) else "NaN"
        mark = '✓' if r['correct'] else '✗'
        lbl  = 'DEG' if r['label'] else 'NON'
        pred = 'DEG' if r['pred'] else 'NON'
        print(f"    {r['poi']:14}  {lbl}→{pred}  p={r['prob_degrader']:.4f}  "
              f"e2={e2s}  {mark}")

    summ_df = pd.DataFrame([{
        'split': 'E2-scalar GNN LOPO-CV',
        'n_folds': n_total,
        'n_correct': n_correct,
        'accuracy': round(acc, 4),
        'auroc_overall': round(auroc, 4),
        'auroc_md_pois': round(auroc_md, 4),
        'auroc_nonmd_pois': round(auroc_nmd, 4),
        'e2_cols': '+'.join(E2_COLS),
    }])

    return pred_df, summ_df


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    if not HAS_DEPS:
        print(f"[ERROR] Missing: {MISSING}")
        sys.exit(1)

    parser = argparse.ArgumentParser()
    parser.add_argument('--static',   default=str(STATIC_PT))
    parser.add_argument('--dynamic',  default=str(DYNAMIC_PT))
    parser.add_argument('--ml_table', default=str(ML_TABLE))
    parser.add_argument('--epochs',   type=int,   default=400)
    parser.add_argument('--hidden',   type=int,   default=64)
    parser.add_argument('--heads',    type=int,   default=4)
    parser.add_argument('--dropout',  type=float, default=0.3)
    parser.add_argument('--lr',       type=float, default=0.001)
    parser.add_argument('--verbose',  action='store_true')
    args = parser.parse_args()

    for p, name in [(args.static, 'static graphs'),
                    (args.dynamic, 'MD graphs'),
                    (args.ml_table, 'ml_table')]:
        if not Path(p).exists():
            print(f"[ERROR] {name} not found: {p}")
            sys.exit(1)

    static_graphs = torch.load(args.static, map_location='cpu', weights_only=False)
    md_graphs     = torch.load(args.dynamic, map_location='cpu', weights_only=False)
    dynamic_map   = {g.poi: g for g in md_graphs}
    e2_df         = load_e2_table(args.ml_table)

    print(f"Loaded {len(static_graphs)} static graphs")
    print(f"Loaded {len(md_graphs)} MD graphs: {sorted(dynamic_map.keys())}")
    print(f"E2 table: {len(e2_df)} POIs, {e2_df.notna().any(axis=1).sum()} with ≥1 E2 value")

    pred_df, summ_df = run_lopo(static_graphs, dynamic_map, e2_df, args)

    GNN_DIR.mkdir(parents=True, exist_ok=True)
    pred_df.to_csv(OUT_PRED, sep='\t', index=False)
    summ_df.to_csv(OUT_SUMM, sep='\t', index=False)
    print(f"\nSaved: {OUT_PRED}")
    print(f"Saved: {OUT_SUMM}")


if __name__ == '__main__':
    main()
