#!/usr/bin/env python3
"""
02_batch_prosettac.py — batch pRosettaC submission for all kinase × PROTAC pairs
==================================================================================
Single stage: N_DOCK=200 → ~10,000 Rosetta-refined poses per kinase.
Matches the protac-md 19-POI training campaign exactly.

pRosettaC interface (Ernie):
  python3 RUNNER Protac_params.txt --n_dock 200 --n_workers W

Runs on Ernie via SSH subprocess.  Fully resumable — checks for score.sc to detect
completion.  Submits up to MAX_PARALLEL jobs simultaneously.

Prerequisites:
  - ~/10PROTAC-MD/prosettac_inputs/{PROTAC}/{GENE}_{UID}/Protac_params.txt
  - pRosettaC runner on Ernie at ERNIE_RUNNER
  - SSH key-based access to Ernie (host alias 'ernie')

Usage:
  python 02_batch_prosettac.py --plan                     # list all runs + ETA
  python 02_batch_prosettac.py --run                      # submit all
  python 02_batch_prosettac.py --run --protac SJF8240     # one PROTAC
  python 02_batch_prosettac.py --run --gene DDR2          # one kinase
  python 02_batch_prosettac.py --status                   # progress report

State machine per run (stored in manifest):
  PENDING → RUNNING → COMPLETE | FAILED
"""

import argparse
import subprocess
import time
import sys
import pandas as pd
from pathlib import Path
from datetime import datetime

ROOT         = Path(__file__).resolve().parents[2]
INPUT_DIR    = ROOT / 'prosettac_inputs'
RUNS_DIR     = ROOT / 'runs'
MANIFEST     = ROOT / 'runs' / 'manifest.tsv'

# ── Ernie config ───────────────────────────────────────────────────────────────
ERNIE_HOST    = 'ernie'
ERNIE_ROOT    = '/home/s1638308/10PROTAC-MD'
ERNIE_PYTHON  = '/home/s1638308/miniconda3/envs/protac/bin/python3'
ERNIE_RUNNER  = '/home/s1638308/protac/scripts/02_rosetta_utils/03_run_local_prosettac.py'
MAX_PARALLEL  = 10    # concurrent kinases (10 × 2 workers = 20 CPUs)
N_WORKERS_PER = 2     # parallel Rosetta workers per job

# ── pRosettaC params ───────────────────────────────────────────────────────────
N_DOCK = 200   # → ~10,000 Rosetta-refined poses per kinase (matches protac-md)


# ── Helpers ────────────────────────────────────────────────────────────────────

def ernie_run(cmd: str, check=True) -> subprocess.CompletedProcess:
    return subprocess.run(['ssh', ERNIE_HOST, cmd], capture_output=True, text=True, check=check)


def ernie_exists(path: str) -> bool:
    return ernie_run(f'test -f {path}', check=False).returncode == 0


def get_all_inputs(protac_filter=None, gene_filter=None):
    pairs = []
    for p_dir in sorted(INPUT_DIR.iterdir()):
        if not p_dir.is_dir():
            continue
        protac = p_dir.name
        if protac_filter and protac != protac_filter:
            continue
        for run_dir in sorted(p_dir.iterdir()):
            if not run_dir.is_dir():
                continue
            if not (run_dir / 'Protac_params.txt').exists():
                continue
            parts = run_dir.name.rsplit('_', 1)
            gene  = parts[0] if len(parts) == 2 else run_dir.name
            uid   = parts[1] if len(parts) == 2 else run_dir.name
            if gene_filter and gene != gene_filter:
                continue
            pairs.append({
                'protac':    protac,
                'gene':      gene,
                'uid':       uid,
                'run_name':  run_dir.name,
                'local_dir': str(run_dir),
            })
    return pairs


def build_manifest(protac_filter=None, gene_filter=None) -> pd.DataFrame:
    pairs = get_all_inputs(protac_filter, gene_filter)
    if not pairs:
        print('No input pairs found.')
        sys.exit(1)

    df = pd.DataFrame(pairs)
    df['status']       = 'PENDING'
    df['bridge_poses'] = -1
    df['bridge_pct']   = -1.0
    df['submit_time']  = ''
    df['done_time']    = ''
    df['ernie_outdir'] = df.apply(
        lambda r: f'{ERNIE_ROOT}/prosettac_inputs/{r["protac"]}/{r["run_name"]}', axis=1
    )

    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    if MANIFEST.exists():
        existing = pd.read_csv(MANIFEST, sep='\t')
        key = ['protac', 'gene', 'uid']
        df = df.set_index(key)
        existing = existing.set_index(key)
        for col in ['status', 'bridge_poses', 'bridge_pct', 'submit_time', 'done_time']:
            if col in existing.columns:
                df[col] = existing[col].reindex(df.index).fillna(df[col])
        df = df.reset_index()
    df.to_csv(MANIFEST, sep='\t', index=False)
    return df


def submit_run(row: pd.Series) -> str:
    local_dir = Path(row['local_dir'])
    out_dir   = row['ernie_outdir']

    ernie_run(f'mkdir -p {out_dir}')
    subprocess.run(
        ['rsync', '-az', '--copy-links', str(local_dir) + '/', f'{ERNIE_HOST}:{out_dir}/'],
        check=True, capture_output=True
    )

    run_cmd = (
        f'export PATH="$HOME/miniconda3/envs/protac/bin:$PATH"; '
        f'export OMP_NUM_THREADS=1; '
        f'cd {out_dir} && '
        f'nohup {ERNIE_PYTHON} {ERNIE_RUNNER} '
        f'Protac_params.txt '
        f'--n_dock {N_DOCK} --n_workers {N_WORKERS_PER} '
        f'> {out_dir}/pRosettaC.log 2>&1 & echo $!'
    )
    r = ernie_run(run_cmd)
    return r.stdout.strip()


def check_complete(row: pd.Series) -> tuple:
    score_sc = f'{row["ernie_outdir"]}/score.sc'
    if not ernie_exists(score_sc):
        return False, -1, -1.0
    local_score = RUNS_DIR / 'scores' / row['protac'] / f'{row["run_name"]}_score.sc'
    local_score.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(['rsync', '-az', f'{ERNIE_HOST}:{score_sc}', str(local_score)],
                   capture_output=True)
    return True, *_count_bridging(local_score)


def _count_bridging(score_sc: Path):
    if not score_sc.exists():
        return 0, 0.0
    try:
        df = pd.read_csv(score_sc, sep=r'\s+', comment='#', header=0)
        cols = [c for c in df.columns if 'bridg' in c.lower()]
        if cols:
            n   = int((df[cols[0]] > 0).sum())
            pct = round(n / len(df) * 100, 2)
            return n, pct
        return 0, 0.0
    except Exception:
        return 0, 0.0


# ── Commands ───────────────────────────────────────────────────────────────────

def cmd_plan(args):
    pairs = get_all_inputs(args.protac, args.gene)
    print(f'Total input pairs: {len(pairs)}  (N_DOCK={N_DOCK}, ~10k poses each)')
    by_protac = {}
    for p in pairs:
        by_protac.setdefault(p['protac'], 0)
        by_protac[p['protac']] += 1
    for k, v in sorted(by_protac.items()):
        wall = v * 32 / 60 / MAX_PARALLEL   # ~32 min/kinase at N_DOCK=200
        print(f'  {k:12s}  {v:4d} pairs  ~{wall:.0f} hrs wall-time')


def cmd_run(args):
    df = build_manifest(args.protac, args.gene)
    pending = df[df['status'] == 'PENDING']
    print(f'{len(pending)} PENDING, '
          f'{(df["status"]=="COMPLETE").sum()} COMPLETE, '
          f'{(df["status"]=="FAILED").sum()} FAILED')

    while True:
        df = pd.read_csv(MANIFEST, sep='\t')

        for idx, row in df[df['status'] == 'RUNNING'].iterrows():
            done, n_bridge, pct = check_complete(row)
            if done:
                df.at[idx, 'status']       = 'COMPLETE'
                df.at[idx, 'bridge_poses'] = n_bridge
                df.at[idx, 'bridge_pct']   = pct
                df.at[idx, 'done_time']    = datetime.now().isoformat()
                print(f'  COMPLETE: {row["protac"]}/{row["gene"]} '
                      f'bridge={n_bridge} ({pct:.1f}%)')

        n_running = int((df['status'] == 'RUNNING').sum())
        for idx, row in df[df['status'] == 'PENDING'].head(MAX_PARALLEL - n_running).iterrows():
            try:
                pid = submit_run(row)
                df.at[idx, 'status']      = 'RUNNING'
                df.at[idx, 'submit_time'] = datetime.now().isoformat()
                print(f'  SUBMITTED: {row["protac"]}/{row["gene"]}  pid={pid}')
            except Exception as e:
                df.at[idx, 'status'] = 'FAILED'
                print(f'  FAILED: {row["protac"]}/{row["gene"]}: {e}')

        df.to_csv(MANIFEST, sep='\t', index=False)

        n_remaining = int((df['status'].isin(['PENDING', 'RUNNING'])).sum())
        n_done      = int((df['status'] == 'COMPLETE').sum())
        print(f'  [{datetime.now().strftime("%H:%M:%S")}] {n_done}/{len(df)} complete, '
              f'{n_remaining} remaining')

        if n_remaining == 0:
            print('All done.')
            break
        time.sleep(60)


def cmd_status(args):
    if not MANIFEST.exists():
        print('No manifest — run --plan first.')
        return
    df = pd.read_csv(MANIFEST, sep='\t')
    print(f'{len(df)} total runs — '
          + ', '.join(f'{s}={int((df["status"]==s).sum())}'
                      for s in ['PENDING', 'RUNNING', 'COMPLETE', 'FAILED']))
    print('\nBy PROTAC:')
    for protac in sorted(df['protac'].unique()):
        sub = df[df['protac'] == protac]
        done = int((sub['status'] == 'COMPLETE').sum())
        bridging = sub[sub['bridge_poses'] > 0]
        print(f'  {protac:12s}  {done}/{len(sub)} complete  '
              f'{len(bridging)} bridging hits')


def main():
    p = argparse.ArgumentParser()
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--plan',   action='store_true')
    g.add_argument('--run',    action='store_true')
    g.add_argument('--status', action='store_true')
    p.add_argument('--protac', default=None)
    p.add_argument('--gene',   default=None)
    args = p.parse_args()

    if args.plan:     cmd_plan(args)
    elif args.run:    cmd_run(args)
    elif args.status: cmd_status(args)


if __name__ == '__main__':
    main()
