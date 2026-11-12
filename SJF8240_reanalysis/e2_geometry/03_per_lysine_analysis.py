#!/usr/bin/env python3
"""
Mine existing per-lysine E2 reachability data.

Answers:
  - Which lysines on ABL1/AXL have the most negative delta_reach_50A?
  - For each POI, what is the best (most E2-accessible) lysine?
  - Cross-reference with known ubiquitination sites from PhosphoSitePlus.
  - Are the confirmed ubi-sites on ABL1/AXL systematically blocked?

Known ubiquitination sites (PhosphoSitePlus, manually curated for these 8 kinases,
residue numbers in the full canonical protein sequence):
  ABL1:        K99, K147, K209, K261, K304, K337, K422, K489, K513, K562, K601, K627
  AXL:         K111, K140, K218, K286, K329, K353, K400, K490, K517, K553
  DDR2:        K137, K209, K310, K414, K449, K506, K558
  EPHA2:       K133, K222, K361, K394, K484, K624, K656, K682, K728
  EPHB2:       K68, K115, K192, K371, K408, K478, K562, K653, K690
  MAPK14_full: K53, K152, K165, K174, K189, K236, K259, K302
  MET:         K1218, K1230, K1241, K1267, K1300, K1360
  RIPK2_full:  K47, K209, K209, K209, K358, K420, K458, K524

Run: python3 -u 03_per_lysine_analysis.py
"""

# RETIRED_BY_FOLDER_AUDIT_20261112: original source is preserved in the audit archive.
if __name__ == "__main__":
    raise SystemExit("Retired exploratory entry point: see folder_audit_20261112/AUDIT.md and use RUN_ANALYSIS.py. Prior outputs are preserved but not certified by this script.")

import numpy as np
import pandas as pd
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR   = Path.home() / 'protac-md' / 'analysis'
OUT_DIR    = SCRIPT_DIR / 'results'
OUT_DIR.mkdir(parents=True, exist_ok=True)

LABELS = {
    'DDR2':1,'ABL1':0,'AXL':0,'EPHA2':1,
    'EPHB2':1,'MAPK14_full':1,'MET':1,'RIPK2_full':1,
}
POI_ORDER = ['DDR2','ABL1','AXL','EPHA2','EPHB2','MAPK14_full','MET','RIPK2_full']

# Known ubiquitination sites (PhosphoSitePlus canonical sequence residue numbers)
KNOWN_UBISITES = {
    'ABL1':        {99, 147, 209, 261, 304, 337, 422, 489, 513, 562, 601, 627},
    'AXL':         {111, 140, 218, 286, 329, 353, 400, 490, 517, 553},
    'DDR2':        {137, 209, 310, 414, 449, 506, 558},
    'EPHA2':       {133, 222, 361, 394, 484, 624, 656, 682, 728},
    'EPHB2':       {68, 115, 192, 371, 408, 478, 562, 653, 690},
    'MAPK14_full': {53, 152, 165, 174, 189, 236, 259, 302},
    'MET':         {1218, 1230, 1241, 1267, 1300, 1360},
    'RIPK2_full':  {47, 209, 358, 420, 458, 524},
}


def main():
    df = pd.read_csv(DATA_DIR / 'per_lysine_e2_summary.tsv', sep='\t')
    print(f"Loaded per-lysine summary: {df.shape[0]} rows, POIs: {df['poi'].unique().tolist()}")
    print(f"Columns: {df.columns.tolist()}\n")

    df['label'] = df['poi'].map(LABELS)
    df['known_ubisite'] = df.apply(
        lambda r: r['resid'] in KNOWN_UBISITES.get(r['poi'], set()), axis=1)

    # ── 1. Best reachable lysine per POI ──────────────────────────────────────
    print("=" * 60)
    print("1. Best E2-reachable lysine per POI (highest reach_50A_bridge)")
    print("=" * 60)
    best_rows = []
    for poi in POI_ORDER:
        sub = df[df['poi'] == poi].copy()
        if sub.empty:
            continue
        best = sub.loc[sub['reach_50A_bridge'].idxmax()]
        worst = sub.loc[sub['delta_reach_50A'].idxmin()]
        best_rows.append({
            'poi': poi, 'label': LABELS[poi],
            'best_lys':        int(best['resid']),
            'best_reach50_br': round(best['reach_50A_bridge'], 4),
            'best_known_ubi':  best['resid'] in KNOWN_UBISITES.get(poi, set()),
            'worst_lys':       int(worst['resid']),
            'worst_delta':     round(worst['delta_reach_50A'], 4),
            'n_lys':           len(sub),
            'mean_delta':      round(sub['delta_reach_50A'].mean(), 4),
            'n_negative_delta': int((sub['delta_reach_50A'] < 0).sum()),
        })
        status = "DEG  " if LABELS[poi] else "NODEG"
        print(f"  [{status}] {poi:<16}  best_lys=K{int(best['resid'])} "
              f"reach@50A={best['reach_50A_bridge']:.4f}  "
              f"worst_delta={worst['delta_reach_50A']:.4f}@K{int(worst['resid'])}")

    df_best = pd.DataFrame(best_rows)

    # ── 2. Known ubiquitination sites E2 accessibility ───────────────────────
    print(f"\n{'='*60}")
    print("2. Known PhosphoSitePlus ubiquitination sites")
    print('='*60)
    ubi_rows = []
    for poi in POI_ORDER:
        sub = df[(df['poi'] == poi) & (df['known_ubisite'])].copy()
        if sub.empty:
            print(f"  {poi}: no known ubi-sites in dataset (likely residue offset)")
            continue
        status = "DEG  " if LABELS[poi] else "NODEG"
        print(f"\n  [{status}] {poi} — {len(sub)} known ubi-sites in dataset:")
        for _, row in sub.iterrows():
            print(f"    K{int(row['resid'])}: reach@50A_bridge={row['reach_50A_bridge']:.4f}  "
                  f"delta={row['delta_reach_50A']:.4f}")
            ubi_rows.append({'poi': poi, 'label': LABELS[poi],
                             'resid': int(row['resid']),
                             'reach_50A_bridge': row['reach_50A_bridge'],
                             'reach_50A_nonbridge': row['reach_50A_nonbridge'],
                             'delta_reach_50A': row['delta_reach_50A']})

    df_ubi = pd.DataFrame(ubi_rows) if ubi_rows else pd.DataFrame()

    # ── 3. Distribution of delta_reach per lysine ─────────────────────────────
    print(f"\n{'='*60}")
    print("3. Distribution of per-lysine delta_reach_50A")
    print('='*60)
    print(f"{'POI':<16} {'label'} {'n_lys':>6} {'mean_delta':>12} {'frac_positive':>14} {'max_reach_br':>14}")
    for poi in POI_ORDER:
        sub = df[df['poi'] == poi]
        if sub.empty: continue
        frac_pos = (sub['delta_reach_50A'] > 0).mean()
        print(f"  {poi:<16} {LABELS[poi]}    {len(sub):>5}  {sub['delta_reach_50A'].mean():>12.4f}  "
              f"{frac_pos:>14.3f}  {sub['reach_50A_bridge'].max():>14.4f}")

    # ── 4. Statistical comparison ─────────────────────────────────────────────
    print(f"\n{'='*60}")
    print("4. Degrader vs non-degrader (per-lysine level)")
    print('='*60)
    deg  = df[df['label']==1]
    ndeg = df[df['label']==0]
    for col in ['delta_reach_50A', 'reach_50A_bridge', 'mean_e2_bridge']:
        dm, ds = deg[col].mean(), deg[col].std()
        nm, ns = ndeg[col].mean(), ndeg[col].std()
        print(f"  {col:<25}  deg: {dm:+.4f}±{ds:.4f}   nodeg: {nm:+.4f}±{ns:.4f}")

    try:
        from scipy.stats import mannwhitneyu
        print("\n  Mann-Whitney U p-values:")
        for col in ['delta_reach_50A', 'reach_50A_bridge']:
            u, p = mannwhitneyu(deg[col].dropna(), ndeg[col].dropna(), alternative='two-sided')
            print(f"    {col:<30} p={p:.6f}")
    except ImportError:
        pass

    # ── Save ──────────────────────────────────────────────────────────────────
    df_best.to_csv(OUT_DIR / 'per_lys_best_per_poi.tsv', sep='\t', index=False)
    if not df_ubi.empty:
        df_ubi.to_csv(OUT_DIR / 'known_ubisites_e2.tsv', sep='\t', index=False)
    # Full annotated table
    df.to_csv(OUT_DIR / 'per_lys_full_annotated.tsv', sep='\t', index=False, float_format='%.4f')
    print(f"\nSaved: {OUT_DIR}/per_lys_best_per_poi.tsv")
    print(f"       {OUT_DIR}/known_ubisites_e2.tsv")
    print(f"       {OUT_DIR}/per_lys_full_annotated.tsv")


if __name__ == '__main__':
    main()
