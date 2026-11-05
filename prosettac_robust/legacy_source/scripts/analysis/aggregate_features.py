#!/usr/bin/env python3
"""
Aggregate per-POI pRosettaC features into a single cross-POI ML feature table.

Reads cluster_features.csv and premd_descriptors.csv from each POI's run
directory and combines them into one table for degradation prediction.

Usage:
    python aggregate_features.py              # use best MD-candidate cluster per POI
    python aggregate_features.py --all-clusters   # one row per cluster per POI

Output:
    analysis/degradation_prediction/ml_features.csv
    analysis/degradation_prediction/ml_features_best.csv   (best cluster per POI)
"""

import sys, csv, json, argparse
import numpy as np
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RUNS_ROOT    = PROJECT_ROOT / 'pRosettaC' / 'runs'
OUT_DIR      = PROJECT_ROOT / 'analysis' / 'degradation_prediction'
OUT_DIR.mkdir(parents=True, exist_ok=True)

# POI registry — must match run_poi_full_analysis.py POI_CONFIGS
POI_REGISTRY = {
    'MET':     {'category': 'high_affinity_degraders',      'degrader': True,  'affinity': 'high'},
    'DDR2':    {'category': 'high_affinity_degraders',      'degrader': True,  'affinity': 'high'},
    'RIPK2':   {'category': 'high_affinity_degraders',      'degrader': True,  'affinity': 'high'},
    'EPHB2':   {'category': 'high_affinity_degraders',      'degrader': True,  'affinity': 'high'},
    'MAPK14':  {'category': 'high_affinity_degraders',      'degrader': True,  'affinity': 'high'},
    'AXL':     {'category': 'high_affinity_no_degradation', 'degrader': False, 'affinity': 'high'},
    'SLK':     {'category': 'high_affinity_no_degradation', 'degrader': False, 'affinity': 'high'},
    'ABL1':    {'category': 'high_affinity_no_degradation', 'degrader': False, 'affinity': 'high'},
    'EPHA2':   {'category': 'high_affinity_no_degradation', 'degrader': False, 'affinity': 'high'},
    'MAP4K5':  {'category': 'high_affinity_no_degradation', 'degrader': False, 'affinity': 'high'},
    'TNIK':    {'category': 'low_affinity_degraders',       'degrader': True,  'affinity': 'low'},
    'PIP4K2C': {'category': 'low_affinity_degraders',       'degrader': True,  'affinity': 'low'},
    'CDK17':   {'category': 'low_affinity_degraders',       'degrader': True,  'affinity': 'low'},
}

# Numeric columns to carry from cluster_features.csv
CLUSTER_NUMERIC = [
    'n_poses', 'score_orig_REU', 'interface_contacts_5A', 'interface_area_est_A2',
    'poi_iface_atoms', 'vhl_iface_atoms', 'polar_hbond_pairs',
    'protac_disconnected', 'protac_overlaps', 'clash_interface', 'clash_protac_protein',
    'n_lys_within_25A_vhl', 'n_lys_e2_optimal',
]

# Numeric columns to carry from premd_descriptors.csv (MD candidates only)
PREMD_NUMERIC = [
    'score_REU', 'vhl_pct', 'poi_pct',
    'protac_n_atoms', 'protac_e2e_A',
    'pt0_atoms', 'pt0_bonds', 'pt0_charge_sum',
    'pt1_atoms', 'pt1_bonds', 'pt1_charge_sum',
    'clash_iface_severe_15A', 'clash_iface_bad_20A', 'clash_iface_warn_25A',
    'clash_protac_severe_15A', 'clash_protac_bad_20A',
    'interface_contacts_5A', 'closest_lys_A', 'n_lys_e2_optimal', 'best_e2_A',
]


def load_csv(path):
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def parse_best_e2_lys(s):
    """Parse 'K1232(37.4A);K1244(37.5A);K1110(51.3A)' → list of (resid, dist)."""
    if not isinstance(s, str) or not s.strip():
        return []
    pairs = []
    for tok in s.split(';'):
        tok = tok.strip()
        if '(' in tok and 'A)' in tok:
            try:
                resid = int(tok[1:tok.index('(')])
                dist  = float(tok[tok.index('(')+1:tok.index('A)')])
                pairs.append((resid, dist))
            except (ValueError, IndexError):
                pass
    return pairs


def cluster_summary(cf_df, poi):
    """Compute per-POI summary stats across all clusters."""
    if cf_df.empty:
        return {}
    cf = cf_df[cf_df['poi'] == poi].copy() if 'poi' in cf_df.columns else cf_df.copy()
    if cf.empty:
        return {}

    # Best cluster = highest n_poses (or md_candidate==1 if flagged)
    mc1 = cf[cf['md_candidate'] == 1] if 'md_candidate' in cf.columns else pd.DataFrame()
    best = mc1.iloc[0] if not mc1.empty else cf.sort_values('n_poses', ascending=False).iloc[0]

    summary = {
        'n_clusters': len(cf),
        'total_poses': cf['n_poses'].sum() if 'n_poses' in cf.columns else 0,
        'best_cluster': best.get('cluster', ''),
        'best_cluster_n_poses': int(best.get('n_poses', 0)),
    }
    for col in CLUSTER_NUMERIC:
        if col in cf.columns:
            summary[f'best_{col}'] = best.get(col, np.nan)
            summary[f'mean_{col}'] = cf[col].mean()
            summary[f'max_{col}']  = cf[col].max()

    # Parse top3 E2 Lys from best cluster
    top3 = parse_best_e2_lys(best.get('best_e2_lys_top3', ''))
    summary['best_e2_lys_dist_1'] = top3[0][1] if len(top3) > 0 else np.nan
    summary['best_e2_lys_dist_2'] = top3[1][1] if len(top3) > 1 else np.nan
    summary['best_e2_lys_dist_3'] = top3[2][1] if len(top3) > 2 else np.nan

    return summary


def premd_summary(pd_df, poi):
    """Take the primary MD candidate row (cluster rank 1)."""
    if pd_df.empty:
        return {}
    df = pd_df[pd_df['poi'] == poi].copy() if 'poi' in pd_df.columns else pd_df.copy()
    if df.empty:
        return {}

    row = df.iloc[0]  # first row = primary MD candidate
    summary = {}
    for col in PREMD_NUMERIC:
        if col in df.columns:
            summary[f'md1_{col}'] = row.get(col, np.nan)
    summary['md1_verdict']   = row.get('md_verdict', '')
    summary['md1_warnings']  = row.get('warnings', '')
    summary['md1_issues']    = row.get('issues', '')
    summary['elongB_present'] = row.get('elongB_present', np.nan)
    summary['elongC_present'] = row.get('elongC_present', np.nan)
    return summary


def aggregate_all():
    rows = []
    for poi, meta in POI_REGISTRY.items():
        run_dir = RUNS_ROOT / meta['category'] / poi
        cf_path = run_dir / 'cluster_features.csv'
        pd_path = run_dir / 'premd_descriptors.csv'

        cf_df = load_csv(cf_path)
        pd_df = load_csv(pd_path)

        if cf_df.empty:
            print(f'  [{poi}] No cluster_features.csv — skipping')
            continue

        row = {
            'poi':       poi,
            'category':  meta['category'],
            'degrader':  int(meta['degrader']),
            'affinity':  meta['affinity'],
        }
        row.update(cluster_summary(cf_df, poi))
        row.update(premd_summary(pd_df, poi))
        rows.append(row)
        print(f'  [{poi}] {len(cf_df)} clusters, md1_contacts={row.get("md1_interface_contacts_5A", "N/A")}')

    if not rows:
        print('No data found. Run run_poi_full_analysis.py --all first.')
        return

    df = pd.DataFrame(rows).sort_values(['affinity', 'degrader'], ascending=[False, False])
    out = OUT_DIR / 'ml_features_best.csv'
    df.to_csv(out, index=False)
    print(f'\nWrote {len(df)} rows → {out}')

    # Print key discriminative features
    print('\n── Key features (best cluster per POI) ──')
    cols = ['poi', 'degrader', 'best_interface_contacts_5A', 'best_n_lys_e2_optimal',
            'best_polar_hbond_pairs', 'best_n_poses']
    available = [c for c in cols if c in df.columns]
    print(df[available].to_string(index=False))

    return df


def aggregate_all_clusters():
    """One row per cluster per POI — for cluster-level ML."""
    rows = []
    for poi, meta in POI_REGISTRY.items():
        run_dir = RUNS_ROOT / meta['category'] / poi
        cf_df = load_csv(run_dir / 'cluster_features.csv')
        if cf_df.empty:
            continue
        cf = cf_df[cf_df['poi'] == poi].copy() if 'poi' in cf_df.columns else cf_df.copy()
        for _, row in cf.iterrows():
            r = {
                'poi':       poi,
                'category':  meta['category'],
                'degrader':  int(meta['degrader']),
                'affinity':  meta['affinity'],
            }
            for col in CLUSTER_NUMERIC + ['cluster', 'representative', 'best_e2_lys_top3']:
                if col in cf.columns:
                    r[col] = row.get(col, np.nan)
            # Parse E2 Lys distances
            top3 = parse_best_e2_lys(row.get('best_e2_lys_top3', ''))
            r['best_e2_lys_dist_1'] = top3[0][1] if len(top3) > 0 else np.nan
            r['best_e2_lys_dist_2'] = top3[1][1] if len(top3) > 1 else np.nan
            r['best_e2_lys_dist_3'] = top3[2][1] if len(top3) > 2 else np.nan
            rows.append(r)

    if not rows:
        print('No cluster data found.')
        return

    df = pd.DataFrame(rows)
    out = OUT_DIR / 'ml_features_all_clusters.csv'
    df.to_csv(out, index=False)
    print(f'Wrote {len(df)} cluster rows → {out}')
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--all-clusters', action='store_true',
                        help='Output one row per cluster (default: one row per POI)')
    args = parser.parse_args()

    print('=== Cross-POI Feature Aggregation ===')
    if args.all_clusters:
        aggregate_all_clusters()
    else:
        aggregate_all()


if __name__ == '__main__':
    main()
