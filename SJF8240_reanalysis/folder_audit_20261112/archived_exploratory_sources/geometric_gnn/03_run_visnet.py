#!/usr/bin/env python3
"""Quick ViSNet-only run + compile final comparison report."""
import sys, json, time, random, warnings
import numpy as np
import pandas as pd
from pathlib import Path

warnings.filterwarnings('ignore')

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import global_mean_pool
from torch_geometric.nn.models import ViSNet
from sklearn.metrics import roc_auc_score

SCRIPT_DIR  = Path(__file__).resolve().parent
GRAPH_DIR   = SCRIPT_DIR.parent / 'corrected_md_only_v1/inputs/elob'
POS_DIR     = SCRIPT_DIR / 'positions'
RESULTS_DIR = SCRIPT_DIR / 'results'
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

LABELS = {
    'DDR2': 1, 'ABL1': 0, 'AXL': 0, 'EPHA2': 1,
    'EPHB2': 1, 'MAPK14_full': 1, 'MET': 1, 'RIPK2_full': 1
}
POI_ORDER = ['DDR2', 'ABL1', 'AXL', 'EPHA2', 'EPHB2', 'MAPK14_full', 'MET', 'RIPK2_full']

EPOCHS = 300
LR     = 1e-3
HIDDEN = 64
N_RUNS = 5
SEED   = 42

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


def load_graphs():
    graphs = []
    for poi in POI_ORDER:
        g = np.load(GRAPH_DIR / f'{poi}.npz')
        p = np.load(POS_DIR / f'{poi}.npz')
        x = torch.tensor(g['x'], dtype=torch.float32)
        pos = np.concatenate([p['vhl_pos'], p['kin_pos']], axis=0)
        pos = torch.tensor(pos, dtype=torch.float32)
        pos = pos - pos.mean(dim=0, keepdim=True)
        n_vhl = len(p['vhl_resids'])
        n_kin = len(p['kin_resids'])
        z = torch.cat([torch.ones(n_vhl, dtype=torch.long),
                       torch.full((n_kin,), 2, dtype=torch.long)])
        graphs.append({'poi': poi, 'label': int(g['label'][0]),
                       'x': x, 'z': z, 'pos': pos})
    return graphs


class ViSNetClassifier(nn.Module):
    def __init__(self, hidden=HIDDEN, num_layers=3, max_z=3, cutoff=12.0, dropout=0.3):
        super().__init__()
        self.visnet = ViSNet(lmax=1, num_layers=num_layers, hidden_channels=hidden,
                             cutoff=cutoff, max_z=max_z, reduce_op='sum', derivative=False)
        self.clf = nn.Sequential(nn.Linear(hidden, hidden), nn.ReLU(),
                                 nn.Dropout(dropout), nn.Linear(hidden, 1))

    def forward(self, data):
        z   = data['z'].to(device)
        pos = data['pos'].to(device)
        batch = torch.zeros(z.size(0), dtype=torch.long, device=device)
        x, v = self.visnet.representation_model(z, pos, batch)
        v_norm = v.norm(dim=1)
        g = global_mean_pool(x + v_norm, batch)
        return self.clf(g).squeeze(-1)


def train_one_fold(train_graphs, test_graph, seed):
    torch.manual_seed(seed)
    model = ViSNetClassifier().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, EPOCHS, eta_min=1e-5)
    pos_weight = torch.tensor([3.0], device=device)

    model.train()
    for _ in range(EPOCHS):
        for g in train_graphs:
            optimizer.zero_grad()
            logit = model(g)
            label = torch.tensor([float(g['label'])], device=device)
            loss = F.binary_cross_entropy_with_logits(logit, label, pos_weight=pos_weight)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
        scheduler.step()

    model.eval()
    with torch.no_grad():
        return torch.sigmoid(model(test_graph)).item()


def main():
    print(f"ViSNet LOPO-CV  Device={device}  Epochs={EPOCHS}  Runs/fold={N_RUNS}")
    graphs = load_graphs()

    results = []
    t0 = time.time()
    for i, test_g in enumerate(graphs):
        train_gs = [g for j, g in enumerate(graphs) if j != i]
        scores = [train_one_fold(train_gs, test_g, SEED + r*100) for r in range(N_RUNS)]
        mean_sc = float(np.mean(scores))
        results.append({'poi': test_g['poi'], 'label': test_g['label'],
                        'score': mean_sc, 'model': 'ViSNet_SE3'})
        correct = round(mean_sc) == test_g['label']
        print(f"  {test_g['poi']:18s} y={test_g['label']}  score={mean_sc:.4f}  "
              f"{'✓' if correct else '✗'}")

    scores_arr = [r['score'] for r in results]
    labels_arr = [r['label'] for r in results]
    n_ok = sum(round(s) == l for s, l in zip(scores_arr, labels_arr))
    try:
        auroc = roc_auc_score(labels_arr, scores_arr)
    except Exception:
        auroc = float('nan')
    print(f"  AUROC={auroc:.4f}  Accuracy={n_ok}/8  ({time.time()-t0:.0f}s)")

    # Compile full comparison with previously saved results
    prev = {
        'Baseline_GAT': {'DDR2':0.8137,'ABL1':0.8131,'AXL':0.8881,'EPHA2':0.8219,
                         'EPHB2':0.7993,'MAPK14_full':0.8057,'MET':0.8052,'RIPK2_full':0.8256,
                         'auroc':0.25,'n_correct':6},
        'GeomGNN_SchNet': {'DDR2':0.9716,'ABL1':0.9818,'AXL':0.9828,'EPHA2':0.9616,
                           'EPHB2':0.9596,'MAPK14_full':0.9619,'MET':0.9645,'RIPK2_full':0.9620,
                           'auroc':0.0,'n_correct':6},
        'EGNN': {'DDR2':0.9842,'ABL1':0.9961,'AXL':0.9892,'EPHA2':0.9871,
                 'EPHB2':0.9792,'MAPK14_full':0.9841,'MET':0.9736,'RIPK2_full':0.9867,
                 'auroc':0.0,'n_correct':6},
    }
    visnet_scores = {r['poi']: r['score'] for r in results}

    # Build comparison table
    rows = []
    for poi in POI_ORDER:
        row = {'poi': poi, 'label': LABELS[poi]}
        for mname, mdata in prev.items():
            row[mname] = mdata[poi]
        row['ViSNet_SE3'] = visnet_scores.get(poi, float('nan'))
        rows.append(row)

    df = pd.DataFrame(rows).set_index('poi')
    df.to_csv(RESULTS_DIR / 'geometric_comparison.tsv', sep='\t', float_format='%.4f')

    summary = []
    for mname, mdata in prev.items():
        summary.append({'model': mname, 'auroc': mdata['auroc'],
                        'n_correct': mdata['n_correct'], 'n_total': 8})
    summary.append({'model': 'ViSNet_SE3', 'auroc': auroc, 'n_correct': n_ok, 'n_total': 8})
    pd.DataFrame(summary).to_csv(RESULTS_DIR / 'geometric_summary.tsv', sep='\t', index=False)

    print("\n" + "="*60)
    print("FULL COMPARISON")
    print(df.to_string())
    print("\nSUMMARY TABLE")
    print(pd.DataFrame(summary).to_string(index=False))
    print(f"\nReference (E2-GNN historical best): AUROC ~0.85, 7/8 correct")
    print(f"Results saved to {RESULTS_DIR}")


if __name__ == '__main__':
    main()
