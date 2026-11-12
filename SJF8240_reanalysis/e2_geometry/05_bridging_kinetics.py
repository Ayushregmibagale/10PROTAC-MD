#!/usr/bin/env python3
"""
Bridging event kinetics from per-frame metrics.csv files.

For each of the 141 MD poses, loads the per-frame distance data and identifies
contiguous bridging events (PROTAC-VHL <5Å AND PROTAC-kinase <5Å).

Measures:
  - Total bridging fraction
  - Number of distinct bridging events (contiguous runs)
  - Mean / max event duration (ns)
  - Mean time between events (ns)
  - Longest single bridging event (ns)

Hypothesis: ABL1/AXL might have high total bridging but short individual events —
transient contacts that don't give E2 enough time to position.

Run: python3 -u 05_bridging_kinetics.py
"""

# RETIRED_BY_FOLDER_AUDIT_20261112: original source is preserved in the audit archive.
if __name__ == "__main__":
    raise SystemExit("Retired exploratory entry point: see folder_audit_20261112/AUDIT.md and use RUN_ANALYSIS.py. Prior outputs are preserved but not certified by this script.")

import numpy as np
import pandas as pd
from pathlib import Path

SCRIPT_DIR  = Path(__file__).resolve().parent
ANA_DIR     = Path.home() / 'protac-md' / 'analysis' / 'random_screen'
OUT_DIR     = SCRIPT_DIR / 'results'
OUT_DIR.mkdir(parents=True, exist_ok=True)

LABELS = {
    'DDR2':1,'ABL1':0,'AXL':0,'EPHA2':1,
    'EPHB2':1,'MAPK14_full':1,'MET':1,'RIPK2_full':1,
}
POI_ORDER  = ['DDR2','ABL1','AXL','EPHA2','EPHB2','MAPK14_full','MET','RIPK2_full']
VHL_CUTOFF = 5.0
KIN_CUTOFF = 5.0


def run_length_encode(arr):
    """Return (values, starts, lengths) for run-length encoding."""
    n = len(arr)
    if n == 0:
        return [], [], []
    vals, starts, lengths = [], [], []
    cur = arr[0]; s = 0; l = 1
    for i in range(1, n):
        if arr[i] == cur:
            l += 1
        else:
            vals.append(cur); starts.append(s); lengths.append(l)
            cur = arr[i]; s = i; l = 1
    vals.append(cur); starts.append(s); lengths.append(l)
    return vals, starts, lengths


def analyse_pose(csv_path):
    df = pd.read_csv(csv_path, comment='#')
    if 'd_protac_vhl_A' not in df.columns or 'd_protac_kinase_A' not in df.columns:
        return None

    t   = df['time_ns'].values
    dvhl = df['d_protac_vhl_A'].values
    dkin = df['d_protac_kinase_A'].values

    bridging = ((dvhl < VHL_CUTOFF) & (dkin < KIN_CUTOFF)).astype(int)
    dt = float(np.median(np.diff(t))) if len(t) > 1 else 0.1  # ns per frame

    frac_bridging = bridging.mean()
    vals, starts, lengths = run_length_encode(bridging)

    bridge_lengths_ns = [l * dt for v, l in zip(vals, lengths) if v == 1]
    n_events = len(bridge_lengths_ns)

    if n_events == 0:
        return {
            'frac_bridging': frac_bridging,
            'n_events': 0,
            'mean_event_ns': 0.0,
            'max_event_ns': 0.0,
            'mean_gap_ns': float(t[-1] - t[0]) if len(t) > 1 else 0.0,
        }

    gap_lengths_ns = [l * dt for v, l in zip(vals, lengths) if v == 0]

    return {
        'frac_bridging':  float(frac_bridging),
        'n_events':       n_events,
        'mean_event_ns':  float(np.mean(bridge_lengths_ns)),
        'max_event_ns':   float(np.max(bridge_lengths_ns)),
        'mean_gap_ns':    float(np.mean(gap_lengths_ns)) if gap_lengths_ns else 0.0,
        'total_time_ns':  float(t[-1]) if len(t) > 0 else 20.0,
    }


def main():
    rows = []
    for poi in POI_ORDER:
        label = LABELS[poi]
        csvs  = sorted(ANA_DIR.glob(f'{poi}/pose*/metrics.csv'))
        print(f"\n[{poi}] label={label}  n_poses={len(csvs)}")

        for csv in csvs:
            pose_id = csv.parts[-2]
            res = analyse_pose(csv)
            if res is None:
                print(f"  {pose_id}: PARSE FAIL")
                continue
            rows.append({'poi': poi, 'label': label, 'pose_id': pose_id, **res})
            print(f"  {pose_id}: frac={res['frac_bridging']:.4f}  "
                  f"n_events={res['n_events']:3d}  "
                  f"mean_dur={res['mean_event_ns']:.3f}ns  "
                  f"max_dur={res['max_event_ns']:.3f}ns", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / 'bridging_kinetics_all_poses.tsv', sep='\t', index=False)

    print(f"\n{'='*70}")
    print("PER-POI SUMMARY")
    print('='*70)
    summary = []
    for poi in POI_ORDER:
        sub = df[df['poi']==poi]
        if sub.empty: continue
        r = {
            'poi': poi, 'label': LABELS[poi],
            'frac_br_mean':     sub['frac_bridging'].mean(),
            'n_events_mean':    sub['n_events'].mean(),
            'mean_event_mean':  sub['mean_event_ns'].mean(),
            'max_event_mean':   sub['max_event_ns'].mean(),
            'mean_gap_mean':    sub['mean_gap_ns'].mean(),
        }
        summary.append(r)
        status = "DEG  " if LABELS[poi] else "NODEG"
        print(f"  [{status}] {poi:<16}  "
              f"frac={r['frac_br_mean']:.4f}  "
              f"n_events={r['n_events_mean']:.1f}  "
              f"mean_dur={r['mean_event_mean']:.3f}ns  "
              f"max_dur={r['max_event_mean']:.3f}ns")

    df_sum = pd.DataFrame(summary)
    df_sum.to_csv(OUT_DIR / 'bridging_kinetics_summary.tsv', sep='\t', index=False, float_format='%.4f')

    # Degrader vs non-degrader
    deg  = df[df['label']==1]
    ndeg = df[df['label']==0]
    print(f"\n{'':20s}  {'degraders':>16}  {'non-degraders':>16}")
    for col in ['frac_bridging','n_events','mean_event_ns','max_event_ns']:
        dm, ds = deg[col].mean(), deg[col].std()
        nm, ns = ndeg[col].mean(), ndeg[col].std()
        print(f"  {col:<22}  {dm:>8.4f}±{ds:.4f}  {nm:>8.4f}±{ns:.4f}")

    try:
        from scipy.stats import mannwhitneyu
        print("\n  Mann-Whitney p-values:")
        for col in ['frac_bridging','n_events','mean_event_ns','max_event_ns']:
            _, p = mannwhitneyu(deg[col], ndeg[col], alternative='two-sided')
            print(f"    {col:<25} p={p:.4f}")
    except ImportError:
        pass

    print(f"\nSaved: {OUT_DIR}/bridging_kinetics_all_poses.tsv")


if __name__ == '__main__':
    main()
