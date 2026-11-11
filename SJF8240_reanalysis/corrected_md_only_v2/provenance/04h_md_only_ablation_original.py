#!/usr/bin/env python3
"""MD-only GNN ablation — the missing piece of the ablation table."""
import sys, numpy as np, pandas as pd
from pathlib import Path

try:
    import torch, torch.nn as nn, torch.nn.functional as F
    from torch_geometric.nn import GATConv, global_mean_pool, global_max_pool
    from sklearn.metrics import roc_auc_score, accuracy_score
except ImportError as e:
    print(f"[ERROR] {e}"); sys.exit(1)

GNN_DIR  = Path.home() / 'protac-md' / 'analysis' / 'gnn'
DYN_PT   = GNN_DIR / 'md_contact_graphs.pt'
STAT_PT  = GNN_DIR / 'ensemble_graphs.pt'
OUT      = GNN_DIR / 'md_only_ablation.tsv'

EPOCHS=300; HIDDEN=64; HEADS=4; DROP=0.3; LR=0.001

class MDEncoder(nn.Module):
    def __init__(self):
        super().__init__()
        self.ep = nn.Linear(3, HIDDEN)
        self.c1 = GATConv(3, HIDDEN//HEADS, heads=HEADS, edge_dim=HIDDEN, dropout=DROP, concat=True)
        self.b1 = nn.BatchNorm1d(HIDDEN)
        self.c2 = GATConv(HIDDEN, HIDDEN, heads=1, edge_dim=HIDDEN, dropout=DROP, concat=False)
        self.b2 = nn.BatchNorm1d(HIDDEN)
        self.drop = nn.Dropout(DROP)
        self.head = nn.Sequential(nn.Linear(HIDDEN*2, HIDDEN), nn.ReLU(), nn.Linear(HIDDEN, 1))

    def forward(self, d):
        x, ei, ea = d.x, d.edge_index, d.edge_attr
        b = getattr(d,'batch',None)
        if b is None: b = torch.zeros(x.size(0), dtype=torch.long, device=x.device)
        ea = F.relu(self.ep(ea))
        x = F.relu(self.b1(self.c1(x,ei,ea))); x = self.drop(x)
        x = F.relu(self.b2(self.c2(x,ei,ea)))
        emb = torch.cat([global_mean_pool(x,b), global_max_pool(x,b)], dim=1)
        return self.head(emb).squeeze(1)

def main():
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Device: {device}, Epochs: {EPOCHS}")

    all_md = torch.load(DYN_PT, map_location='cpu', weights_only=False)
    all_st = torch.load(STAT_PT, map_location='cpu', weights_only=False)
    clean_pois = {g.poi for g in all_st if not g.artifact}

    # Only MD POIs that are also clean (exclude TNIK artifact)
    md_clean = [g for g in all_md if g.poi in clean_pois]
    print(f"MD-clean POIs ({len(md_clean)}): {[g.poi for g in md_clean]}")
    print(f"Labels: {[(g.poi, int(g.y.item())) for g in md_clean]}\n")

    rows = []
    for i, test in enumerate(md_clean):
        train = [g for g in md_clean if g.poi != test.poi]
        model = MDEncoder().to(device)
        opt   = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=1e-4)
        sch   = torch.optim.lr_scheduler.CosineAnnealingLR(opt, EPOCHS)
        for ep in range(EPOCHS):
            model.train()
            np.random.shuffle(train)
            for g in train:
                opt.zero_grad()
                out = model(g.to(device))
                loss = F.binary_cross_entropy_with_logits(
                    out, torch.tensor([float(g.y.item())], device=device))
                loss.backward(); opt.step()
            sch.step()
        model.eval()
        with torch.no_grad():
            prob = float(torch.sigmoid(model(test.to(device))).cpu())
        lbl = int(test.y.item())
        ok  = int((prob>=0.5)==bool(lbl))
        print(f"  [{test.poi:14}] label={lbl}  prob={prob:.4f}  "
              f"pred={'DEG' if prob>=0.5 else 'NON'}  {'✓' if ok else '✗'}", flush=True)
        rows.append({'poi': test.poi, 'label': lbl, 'prob': round(prob,4), 'correct': ok})

    df = pd.DataFrame(rows)
    auroc = roc_auc_score(df.label, df.prob)
    acc   = df.correct.sum()
    print(f"\n{'='*50}")
    print(f"MD-only LOPO-CV ({len(df)} MD POIs)")
    print(f"  AUROC = {auroc:.4f}")
    print(f"  ACC   = {acc}/{len(df)}")
    print(f"\nKnown comparison:")
    print(f"  pRosettaC-only (17 POIs): AUROC = 0.185")
    print(f"  MD-only        (8 POIs):  AUROC = {auroc:.4f}  ← NEW")
    print(f"  Dual-encoder   (17 POIs): AUROC = 0.850")
    df.to_csv(OUT, sep='\t', index=False)
    print(f"\nSaved: {OUT}")

if __name__ == '__main__':
    main()
