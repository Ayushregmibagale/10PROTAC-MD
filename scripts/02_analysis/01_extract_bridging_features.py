#!/usr/bin/env python3
"""
01_extract_bridging_features.py — extract pRosettaC bridging features for all runs
====================================================================================
Rsyncs completed pRosettaC outputs from Ernie, extracts bridging fraction and
geometry features, and writes a master feature table ready for ML scoring.

Input:
  Stage 2 manifest (runs/manifest_stage2.tsv)
  pRosettaC score.sc + complex PDBs on Ernie

Output:
  analysis/bridging_features_all.tsv   — one row per (PROTAC, gene, klifs_id)
  analysis/bridging_features_{PROTAC}.tsv

Usage:
  python 01_extract_bridging_features.py --rsync         # pull outputs from Ernie
  python 01_extract_bridging_features.py --extract       # compute features
  python 01_extract_bridging_features.py --run           # rsync + extract
  python 01_extract_bridging_features.py --status
"""

import argparse
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path

ROOT         = Path(__file__).resolve().parents[2]
RUNS_DIR     = ROOT / 'runs'
ANA_DIR      = ROOT / 'analysis'
MANIFEST_S2  = RUNS_DIR / 'manifest_stage2.tsv'
OUT_ALL      = ANA_DIR / 'bridging_features_all.tsv'

ERNIE_HOST   = 'ernie'
ERNIE_ROOT   = '/home/ayush/10PROTAC-MD'

BRIDGE_DIST  = 5.0
VHL_MAX_RES  = 103
KIN_MIN_RES  = 200

# Geometry feature columns (mirrors protac-md feature set)
FEATURE_COLS = [
    'bridge_pct',
    'bridge_n_poses',
    'mean_vhl_dist',
    'mean_kin_dist',
    'mean_linker_rmsd',
    'rosetta_score_mean',
    'rosetta_score_bridge_mean',
]


# ── Rsync ──────────────────────────────────────────────────────────────────────

def rsync_completed(args):
    if not MANIFEST_S2.exists():
        print('No stage 2 manifest.')
        return
    df = pd.read_csv(MANIFEST_S2, sep='\t')
    complete = df[df['status'] == 'COMPLETE']
    print(f'Rsyncing {len(complete)} completed runs from Ernie ...')
    ok = fail = skip = 0
    for _, row in complete.iterrows():
        local_dir = RUNS_DIR / 'stage2' / row['protac'] / row['run_name']
        if local_dir.exists() and (local_dir / 'score.sc').exists():
            skip += 1
            continue
        local_dir.mkdir(parents=True, exist_ok=True)
        ernie_dir = row['ernie_outdir']
        result = subprocess.run(
            ['rsync', '-az', '--include=score.sc', '--include=*.pdb',
             '--exclude=*', f'{ERNIE_HOST}:{ernie_dir}/', str(local_dir) + '/'],
            capture_output=True
        )
        if result.returncode == 0:
            ok += 1
        else:
            fail += 1
            print(f'  rsync failed: {row["protac"]}/{row["gene"]}')
    print(f'rsync done — ok: {ok}, skip: {skip}, fail: {fail}')


# ── Feature extraction ─────────────────────────────────────────────────────────

def parse_score_sc(score_path: Path) -> pd.DataFrame:
    try:
        lines = [l for l in score_path.read_text().splitlines()
                 if not l.startswith('#') and l.strip()]
        if len(lines) < 2:
            return pd.DataFrame()
        header = lines[0].split()
        rows   = [l.split() for l in lines[1:] if len(l.split()) == len(header)]
        df     = pd.DataFrame(rows, columns=header)
        for col in df.columns:
            try:
                df[col] = pd.to_numeric(df[col])
            except Exception:
                pass
        return df
    except Exception:
        return pd.DataFrame()


def extract_features_one(run_dir: Path, protac: str, gene: str) -> dict:
    score_path = run_dir / 'score.sc'
    if not score_path.exists():
        return {}

    sc = parse_score_sc(score_path)
    if sc.empty:
        return {}

    # Bridging columns: vhl_dist, kin_dist — adapt to actual pRosettaC column names
    vhl_col = next((c for c in sc.columns if 'vhl' in c.lower() and 'dist' in c.lower()), None)
    kin_col = next((c for c in sc.columns if 'kin' in c.lower() and 'dist' in c.lower()), None)
    score_col = next((c for c in sc.columns if 'total' in c.lower() or 'score' in c.lower()), None)

    feats = {'protac': protac, 'gene': gene, 'n_poses': len(sc)}

    if vhl_col and kin_col:
        bridging = (sc[vhl_col] <= BRIDGE_DIST) & (sc[kin_col] <= BRIDGE_DIST)
        feats['bridge_n_poses'] = int(bridging.sum())
        feats['bridge_pct']     = round(bridging.mean() * 100, 2)
        feats['mean_vhl_dist']  = round(sc[vhl_col].mean(), 3)
        feats['mean_kin_dist']  = round(sc[kin_col].mean(), 3)
        if score_col:
            feats['rosetta_score_mean']        = round(sc[score_col].mean(), 3)
            feats['rosetta_score_bridge_mean'] = (
                round(sc.loc[bridging, score_col].mean(), 3) if bridging.any() else float('nan')
            )
    else:
        feats['bridge_n_poses'] = -1
        feats['bridge_pct']     = -1.0

    return feats


def cmd_extract(args):
    if not MANIFEST_S2.exists():
        print('No stage 2 manifest.')
        return

    df = pd.read_csv(MANIFEST_S2, sep='\t')
    complete = df[df['status'] == 'COMPLETE']
    print(f'Extracting features for {len(complete)} runs ...')

    records = []
    for _, row in complete.iterrows():
        local_dir = RUNS_DIR / 'stage2' / row['protac'] / row['run_name']
        feats = extract_features_one(local_dir, row['protac'], row['gene'])
        if feats:
            records.append(feats)
        else:
            print(f'  WARN: no features for {row["protac"]}/{row["gene"]}')

    if not records:
        print('No features extracted.')
        return

    out_df = pd.DataFrame(records)
    ANA_DIR.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(OUT_ALL, sep='\t', index=False)
    print(f'Saved {len(out_df)} rows → {OUT_ALL}')

    # Per-PROTAC files
    for protac in out_df['protac'].unique():
        sub = out_df[out_df['protac'] == protac]
        sub.to_csv(ANA_DIR / f'bridging_features_{protac}.tsv', sep='\t', index=False)
        n_bridging = int((sub['bridge_pct'] > 0).sum())
        print(f'  {protac:12s}  {len(sub)} runs, {n_bridging} with any bridging')


def cmd_status(args):
    if OUT_ALL.exists():
        df = pd.read_csv(OUT_ALL, sep='\t')
        print(f'Feature table: {len(df)} rows')
        for protac in df['protac'].unique():
            sub = df[df['protac'] == protac]
            top = sub.nlargest(5, 'bridge_pct')[['gene', 'bridge_pct']].to_string(index=False)
            print(f'\n  {protac} — top 5 by bridge_pct:\n{top}')
    else:
        print('No feature table yet.')


def main():
    p = argparse.ArgumentParser()
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--rsync',   action='store_true')
    g.add_argument('--extract', action='store_true')
    g.add_argument('--run',     action='store_true', help='rsync + extract')
    g.add_argument('--status',  action='store_true')
    args = p.parse_args()

    if args.rsync:
        rsync_completed(args)
    elif args.extract:
        cmd_extract(args)
    elif args.run:
        rsync_completed(args)
        cmd_extract(args)
    elif args.status:
        cmd_status(args)


if __name__ == '__main__':
    main()
