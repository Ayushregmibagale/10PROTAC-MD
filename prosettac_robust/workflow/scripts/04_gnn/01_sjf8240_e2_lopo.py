#!/usr/bin/env python3
"""
04g — Dual-encoder GNN + E2 reachability scalars injected at fusion stage.

Identical to 04d (dual-encoder) except the fusion MLP receives additional
graph-level E2 scalar features:
  - e2_mean_min_dist_A      : producer-defined per-pose nearest E2-point/lysine distance mean (Å)
  - e2_delta_reach_50A      : Δ E2 reach at 50 Å (bridging vs non-bridging)
  - e2_delta_reach_80A      : Δ E2 reach at 80 Å
  - lk_lys_nearest_e2_dist_A: rank-selected static nearest lysine NZ to its modeled E2 point (Å)

Missing descriptors remain NaN before fold normalization. Absence of MD and absence
of bridging are different conditions; exact scalar definitions/conditioning belong
in the reviewed input provenance. Historical all-frame means are not implicitly
redefined as bridging-conditioned means.
Imputation: per-fold z-normalise on training set; NaN → 0.0 after normalisation
Missing-value imputation does not eliminate missingness confounding; coverage
and missing-data sensitivity must be reviewed separately.

Hypothesis: ABL1/AXL/EPHA2 may be separable from DEGs via E2 accessibility
patterns not captured by graph topology alone.

Usage:
  python3 scripts/02_analysis/04g_e2_scalar_gnn.py
  python3 scripts/02_analysis/04g_e2_scalar_gnn.py --epochs 400 --hidden 64
"""
import argparse
import json,random
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from gnn_contract import verify_gnn_release,select_training_graphs,normalize_fold_scalars,fold_seed,validate_feature_graph
from protac_qc import ValidationError,sha256
from gnn_contract import configure_native_pyg
DISABLED_OPTIONAL_GNN_BACKENDS=configure_native_pyg()
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

E2_COLS = [
    'e2_mean_min_dist_A',
    'e2_delta_reach_50A',
    'e2_delta_reach_80A',
    'lk_lys_nearest_e2_dist_A',
]

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
    if 'feature_schema' not in df or not df['feature_schema'].eq('SJF8240-periodic-contacts-v1').all():raise ValidationError('Historical or incompatible scalar feature schema')
    if 'ligand_inchikey' not in df or not df['ligand_inchikey'].eq('SCUZPTNDZHEHLB-ZWUDBOSPSA-N').all():raise ValidationError('Scalar feature chemistry is unverified')
    if df['poi'].duplicated().any():raise ValidationError('Duplicate POI in scalar input')
    df = df[df['prostac_artifact'] == 0].copy()
    df = df.set_index('poi')
    return df[E2_COLS]


def get_e2_vec(poi, e2_df, train_pois, device):
    """
    Return z-normalised E2 scalar tensor for `poi`.
    Normalisation stats computed from training POIs only.
    NaN (non-MD POIs, or missing features) → 0.0 after normalisation.
    """
    normed=normalize_fold_scalars(poi,e2_df,train_pois,E2_COLS)
    return torch.tensor(normed,dtype=torch.float32,device=device).unsqueeze(0)


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
    return float(loss.detach())


@torch.no_grad()
def predict(model, static_g, dynamic_g, e2_vec, device):
    model.eval()
    static_g = static_g.to(device)
    dynamic_g = dynamic_g.to(device) if dynamic_g is not None else None
    out = model(static_g, dynamic_g, e2_vec)
    return float(torch.sigmoid(out).cpu())


# ── LOPO-CV ───────────────────────────────────────────────────────────────────

def save_fold_checkpoint(path,model,poi,train_pois,seed,e2_df,e2_test,args,prediction):
    """Retain inference state and exact training-only scalar reference rows."""
    path=Path(path)
    if any(path.resolve().is_relative_to(root) for root in [Path('/home/ayushregmibagale/protac-md'),Path('/mnt/d/PROTAC degradation Prediction')]):raise ValidationError('Checkpoint output points into original work')
    if path.exists():raise ValidationError('Fold checkpoint already exists')
    if not np.isfinite(prediction) or not 0<=prediction<=1:raise ValidationError('Invalid model score')
    if poi in train_pois:raise ValidationError('Held-out target leaked into checkpoint training cohort')
    reference=e2_df.loc[e2_df.index.isin(train_pois),E2_COLS]
    checkpoint={'format_version':1,'purpose':'EXPLORATORY_FOLD_INFERENCE_NOT_CALIBRATED_PROBABILITY',
        'held_out_poi':poi,'training_pois':list(train_pois),'seed':int(seed),
        'architecture':{'n_e2_scalars':len(E2_COLS),'hidden':args.hidden,'heads':args.heads,'dropout':args.dropout},
        'model_state_dict':{name:value.detach().cpu().clone() for name,value in model.state_dict().items()},
        'scalar_columns':list(E2_COLS),'scalar_reference_pois':list(reference.index),
        'scalar_reference_values':torch.tensor(reference.to_numpy(dtype=float),dtype=torch.float64),
        'scalar_normalization':'Training rows only; finite mean and sample SD (ddof=1); all-missing mean 0; singleton/zero scale 1; target NaN becomes normalized 0',
        'held_out_normalized_scalars':e2_test.detach().cpu().clone(),
        'held_out_score':float(prediction),'feature_schema':'SJF8240-periodic-contacts-v1',
        'source_sha256':sha256(Path(__file__)),'contract_sha256':sha256(Path(__file__).resolve().parents[1]/'lib/gnn_contract.py'),
        'torch_version':str(torch.__version__)}
    path.parent.mkdir(parents=True,exist_ok=True);torch.save(checkpoint,path)
    return {'path':str(path),'sha256':sha256(path)}


def run_lopo(static_graphs, dynamic_map, e2_df, args):
    device = 'cpu'
    print(f"Device: {device}")
    print(f"E2 scalars: {E2_COLS}\n")

    clean = static_graphs  # Already filtered against the reviewed exact cohort.
    md_pois  = set(dynamic_map.keys())
    n_e2     = len(E2_COLS)

    print(f"Clean POIs: {len(clean)}")
    print(f"POIs with MD:  {sorted(md_pois)}")
    print(f"POIs with E2:  {sorted(e2_df.dropna(how='all').index.tolist())}\n")

    pred_rows = []

    for fold_i, test_g in enumerate(clean):
        poi        = test_g.poi
        stage_seed=fold_seed(args.seed,poi)
        random.seed(stage_seed);np.random.seed(stage_seed);torch.manual_seed(stage_seed)
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
        checkpoint=save_fold_checkpoint(args.output/'fold_checkpoints'/f'fold_{fold_i+1:04d}.pt',model,poi,train_pois,stage_seed,e2_df,e2_test,args,prob)
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
            'label': label, 'prob_degrader': float(prob),
            'fold_checkpoint':str(Path(checkpoint['path']).relative_to(args.output)),
            'fold_checkpoint_sha256':checkpoint['sha256'],
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
    parser.add_argument('--static',required=True)
    parser.add_argument('--dynamic',required=True)
    parser.add_argument('--ml_table',required=True)
    parser.add_argument('--epochs',   type=int,   default=400)
    parser.add_argument('--hidden',   type=int,   default=64)
    parser.add_argument('--heads',    type=int,   default=4)
    parser.add_argument('--dropout',  type=float, default=0.3)
    parser.add_argument('--lr',       type=float, default=0.001)
    parser.add_argument('--verbose',  action='store_true')
    parser.add_argument('--seed',type=int,required=True)
    parser.add_argument('--release',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.epochs<1 or args.hidden<1 or args.heads<1 or not 0<=args.dropout<1 or (not np.isfinite(args.lr)) or args.lr<=0 or args.hidden%args.heads!=0:
        raise ValidationError('Invalid GNN training budget')
    old_roots=[Path('/home/ayushregmibagale/protac-md'),Path('/mnt/d/PROTAC degradation Prediction')]
    if any(args.output.resolve().is_relative_to(p) for p in old_roots):raise ValidationError('Output points into original work')
    if args.output.exists():raise ValidationError('GNN output directory already exists')
    release=verify_gnn_release(args.release,{'static':Path(args.static),'dynamic':Path(args.dynamic),'ml_table':Path(args.ml_table)},'SCUZPTNDZHEHLB-ZWUDBOSPSA-N')
    source_digest=sha256(Path(__file__));contract_digest=sha256(Path(__file__).resolve().parents[1]/'lib/gnn_contract.py');release_digest=sha256(args.release)
    torch.set_num_threads(1);torch.use_deterministic_algorithms(True)

    for p, name in [(args.static, 'static graphs'),
                    (args.dynamic, 'MD graphs'),
                    (args.ml_table, 'ml_table')]:
        if not Path(p).exists():
            print(f"[ERROR] {name} not found: {p}")
            sys.exit(1)

    static_graphs = torch.load(args.static, map_location='cpu', weights_only=False)
    md_graphs     = torch.load(args.dynamic, map_location='cpu', weights_only=False)
    static_graphs,excluded=select_training_graphs(static_graphs,release['experimental_labels'])
    for g in static_graphs:validate_feature_graph(g,9)
    dynamic_map={}
    for g in md_graphs:
        if g.poi not in release['experimental_labels']:continue
        if g.poi in dynamic_map or float(g.y.item())!=release['experimental_labels'][g.poi]:
            raise ValidationError('Duplicate/mislabeled dynamic graph: '+g.poi)
        validate_feature_graph(g,3)
        dynamic_map[g.poi]=g
    e2_df         = load_e2_table(args.ml_table)

    print(f"Loaded {len(static_graphs)} static graphs")
    print(f"Loaded {len(md_graphs)} MD graphs: {sorted(dynamic_map.keys())}")
    print(f"E2 table: {len(e2_df)} POIs, {e2_df.notna().any(axis=1).sum()} with ≥1 E2 value")

    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'training_intent.json').write_text(json.dumps({'status':'INTENT_ONLY_NOT_COMPLETION',
        'seed':args.seed,'source_sha256':source_digest,'contract_sha256':contract_digest,'release_sha256':release_digest,
        'inputs_sha256':{k:sha256(Path(getattr(args,k))) for k in ['static','dynamic','ml_table']},
        'settings':{k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
        'completion_requires':'receipt.json with EXPLORATORY_LOPO_COMPLETE'},indent=2)+'\n')
    pred_df, summ_df = run_lopo(static_graphs, dynamic_map, e2_df, args)

    if sha256(Path(__file__))!=source_digest or sha256(Path(__file__).resolve().parents[1]/'lib/gnn_contract.py')!=contract_digest or sha256(args.release)!=release_digest:raise ValidationError('GNN source/release changed during training')
    verify_gnn_release(args.release,{'static':Path(args.static),'dynamic':Path(args.dynamic),'ml_table':Path(args.ml_table)},'SCUZPTNDZHEHLB-ZWUDBOSPSA-N')
    pred_df.to_csv(args.output/'predictions.tsv',sep='\t',index=False)
    summ_df.to_csv(args.output/'summary.tsv',sep='\t',index=False)
    receipt={'status':'EXPLORATORY_LOPO_COMPLETE','seed':args.seed,'device':'cpu','deterministic_algorithms':True,
        'source_sha256':sha256(Path(__file__)),'contract_sha256':sha256(Path(__file__).resolve().parents[1]/'lib/gnn_contract.py'),
        'release_sha256':sha256(args.release),'cohort':sorted(release['experimental_labels']),'excluded_unknown_or_outside_cohort':excluded,
        'dynamic_graph_targets':sorted(dynamic_map),'missing_dynamic_targets':sorted(set(release['experimental_labels'])-set(dynamic_map)),
        'scalar_observed_counts':{c:int(e2_df.loc[e2_df.index.isin(release['experimental_labels']),c].notna().sum()) for c in E2_COLS},
        'settings':{k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
        'inputs_sha256':{k:sha256(Path(getattr(args,k))) for k in ['static','dynamic','ml_table']},
        'outputs_sha256':{str(p.relative_to(args.output)):sha256(p) for p in args.output.rglob('*') if p.is_file()},
        'fold_model_checkpoints_saved':True,'scores_rounded_before_metrics':False,
        'disabled_optional_GNN_backends':DISABLED_OPTIONAL_GNN_BACKENDS,'torch_version':torch.__version__,'numpy_version':np.__version__,
        'generalization_certified':False,'limitations':['LOPO on the development cohort is exploratory; independent held-out experimental validation remains required.','A model score is not a calibrated degradation probability.']}
    (args.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(f"Saved isolated exploratory results: {args.output}")


if __name__ == '__main__':
    main()
