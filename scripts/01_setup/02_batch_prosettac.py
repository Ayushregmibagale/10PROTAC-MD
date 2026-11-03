#!/usr/bin/env python3
"""
02_batch_prosettac.py — batch pRosettaC submission for all kinase × PROTAC pairs
==================================================================================
Two-stage screen:
  Stage 1 (fast):  n_poses=50  → identify any bridging poses (bridge_pct > 0)
  Stage 2 (deep):  n_poses=500 → top bridging hits from Stage 1 only

Runs on Ernie via SSH subprocess. Fully resumable — checks for score.sc to detect
completion. Submits up to MAX_PARALLEL jobs simultaneously.

Prerequisites:
  - ~/10PROTAC-MD/prosettac_inputs/{PROTAC}/{GENE}_{ID}/complex_input.pdb
  - pRosettaC installed on Ernie at ERNIE_PROSETTAC_DIR
  - SSH access to Ernie without password (key-based)

Usage:
  python 02_batch_prosettac.py --plan                         # list all runs
  python 02_batch_prosettac.py --stage1                       # run stage 1 (n=50)
  python 02_batch_prosettac.py --stage2                       # run stage 2 (top hits, n=500)
  python 02_batch_prosettac.py --stage1 --protac SJF8240      # one PROTAC
  python 02_batch_prosettac.py --stage1 --gene DDR2           # one kinase
  python 02_batch_prosettac.py --status                       # progress report
  python 02_batch_prosettac.py --filter-stage1                # identify stage2 candidates

State machine per run (stored in manifest):
  PENDING → SUBMITTED → RUNNING → COMPLETE | FAILED
"""

import argparse
import subprocess
import time
import sys
import json
import pandas as pd
from pathlib import Path
from datetime import datetime

ROOT         = Path(__file__).resolve().parents[2]
INPUT_DIR    = ROOT / 'prosettac_inputs'
RUNS_DIR     = ROOT / 'runs'
MANIFEST_S1  = ROOT / 'runs' / 'manifest_stage1.tsv'
MANIFEST_S2  = ROOT / 'runs' / 'manifest_stage2.tsv'

# ── Ernie config ───────────────────────────────────────────────────────────────
ERNIE_HOST          = 'ernie'                          # SSH host alias in ~/.ssh/config
ERNIE_ROOT          = '/home/ayush/10PROTAC-MD'        # project root on Ernie
ERNIE_PROSETTAC     = '/home/ayush/pRosettaC'          # pRosettaC install on Ernie
ERNIE_PYTHON        = '/home/ayush/.conda/envs/pRosettaC/bin/python3'
MAX_PARALLEL        = 12                               # max concurrent pRosettaC jobs on Ernie

# ── pRosettaC params ───────────────────────────────────────────────────────────
STAGE1_POSES = 50
STAGE2_POSES = 500
BRIDGE_DIST  = 5.0     # Å, VHL-warhead and kinase-warhead simultaneously
VHL_MAX_RES  = 103     # last VHL residue in combined PDB
KIN_MIN_RES  = 200     # first kinase residue in combined PDB (renumbered)

# Stage 2 selection criterion: at least 1 bridging pose out of 50
STAGE2_MIN_BRIDGE_POSES = 1


# ── Helpers ────────────────────────────────────────────────────────────────────

def ernie_run(cmd: str, check=True) -> subprocess.CompletedProcess:
    full = ['ssh', ERNIE_HOST, cmd]
    return subprocess.run(full, capture_output=True, text=True, check=check)


def ernie_exists(path: str) -> bool:
    r = ernie_run(f'test -f {path}', check=False)
    return r.returncode == 0


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
            complex_pdb = run_dir / 'complex_input.pdb'
            if not complex_pdb.exists():
                continue
            gene = run_dir.name.rsplit('_', 1)[0]
            klifs_id = run_dir.name.rsplit('_', 1)[-1]
            if gene_filter and gene != gene_filter:
                continue
            pairs.append({
                'protac':    protac,
                'gene':      gene,
                'klifs_id':  klifs_id,
                'run_name':  run_dir.name,
                'input_pdb': str(complex_pdb),
            })
    return pairs


def build_manifest(stage: int, protac_filter=None, gene_filter=None) -> pd.DataFrame:
    manifest_path = MANIFEST_S1 if stage == 1 else MANIFEST_S2
    pairs = get_all_inputs(protac_filter, gene_filter)

    if stage == 2:
        # Filter to stage 1 successes with bridging poses
        s1 = pd.read_csv(MANIFEST_S1, sep='\t') if MANIFEST_S1.exists() else pd.DataFrame()
        if s1.empty:
            print('No stage 1 manifest — run stage 1 first.')
            sys.exit(1)
        passing = set(
            zip(s1[s1['bridge_poses'] >= STAGE2_MIN_BRIDGE_POSES]['protac'],
                s1[s1['bridge_poses'] >= STAGE2_MIN_BRIDGE_POSES]['gene'])
        )
        pairs = [p for p in pairs if (p['protac'], p['gene']) in passing]
        print(f'Stage 2 candidates (bridge_poses >= {STAGE2_MIN_BRIDGE_POSES}): {len(pairs)}')

    df = pd.DataFrame(pairs)
    df['status']       = 'PENDING'
    df['n_poses']      = STAGE1_POSES if stage == 1 else STAGE2_POSES
    df['bridge_poses'] = -1
    df['bridge_pct']   = -1.0
    df['submit_time']  = ''
    df['done_time']    = ''
    df['ernie_outdir'] = df.apply(
        lambda r: f'{ERNIE_ROOT}/runs/stage{stage}/{r["protac"]}/{r["run_name"]}', axis=1
    )

    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    if manifest_path.exists():
        existing = pd.read_csv(manifest_path, sep='\t')
        # Preserve status of already-submitted runs
        key = ['protac', 'gene', 'klifs_id']
        df = df.set_index(key)
        existing = existing.set_index(key)
        for col in ['status', 'bridge_poses', 'bridge_pct', 'submit_time', 'done_time']:
            if col in existing.columns:
                df[col] = existing[col].reindex(df.index).fillna(df[col])
        df = df.reset_index()
    df.to_csv(manifest_path, sep='\t', index=False)
    return df


def submit_run(row: pd.Series, stage: int) -> str:
    """Submit one pRosettaC job on Ernie via nohup. Returns job pid string."""
    n_poses  = row['n_poses']
    out_dir  = row['ernie_outdir']
    # rsync input PDB to Ernie
    ernie_input = f'{out_dir}/complex_input.pdb'
    rsync_cmd = [
        'rsync', '-az', row['input_pdb'],
        f'{ERNIE_HOST}:{ernie_input}'
    ]
    subprocess.run(['ssh', ERNIE_HOST, f'mkdir -p {out_dir}'], check=True,
                   capture_output=True)
    subprocess.run(rsync_cmd, check=True, capture_output=True)

    run_cmd = (
        f'cd {ERNIE_PROSETTAC} && '
        f'nohup {ERNIE_PYTHON} pRosettaC.py '
        f'--receptor {ernie_input} '
        f'--n_dock {n_poses} '
        f'--out_dir {out_dir} '
        f'--vhl_chain A --kin_chain D '
        f'> {out_dir}/pRosettaC.log 2>&1 & echo $!'
    )
    r = ernie_run(run_cmd)
    pid = r.stdout.strip()
    return pid


def check_complete(row: pd.Series) -> tuple:
    """Return (is_complete, bridge_poses, bridge_pct)."""
    score_sc = f'{row["ernie_outdir"]}/score.sc'
    if not ernie_exists(score_sc):
        return False, -1, -1.0
    # rsync score.sc back and count bridging poses
    local_score = RUNS_DIR / 'stage_scores' / row['protac'] / f'{row["run_name"]}_score.sc'
    local_score.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(['rsync', '-az', f'{ERNIE_HOST}:{score_sc}', str(local_score)],
                   capture_output=True)
    bridge_poses, bridge_pct = _count_bridging(local_score, row)
    return True, bridge_poses, bridge_pct


def _count_bridging(score_sc: Path, row: pd.Series):
    """Parse score.sc to count poses with PROTAC bridging both VHL and kinase."""
    if not score_sc.exists():
        return 0, 0.0
    try:
        df = pd.read_csv(score_sc, sep='\s+', comment='#', header=0)
        # bridging columns depend on pRosettaC version — adapt as needed
        bridge_col = [c for c in df.columns if 'bridg' in c.lower() or 'bridge' in c.lower()]
        if bridge_col:
            n = int((df[bridge_col[0]] > 0).sum())
            pct = n / len(df) * 100
        else:
            n, pct = 0, 0.0
        return n, round(pct, 2)
    except Exception:
        return 0, 0.0


# ── Commands ───────────────────────────────────────────────────────────────────

def cmd_plan(args):
    pairs = get_all_inputs(args.protac, args.gene)
    print(f'Total input pairs ready for pRosettaC: {len(pairs)}')
    by_protac = {}
    for p in pairs:
        by_protac.setdefault(p['protac'], 0)
        by_protac[p['protac']] += 1
    for k, v in sorted(by_protac.items()):
        s1_hrs = v * (STAGE1_POSES / 50 * 8) / 60   # ~8 min per 50-pose run
        print(f'  {k:12s}  {v:4d} pairs  ~{s1_hrs:.0f} CPU-hrs (stage1)')


def cmd_stage(args, stage: int):
    df = build_manifest(stage, args.protac, args.gene)
    pending = df[df['status'] == 'PENDING']
    print(f'Stage {stage}: {len(pending)} PENDING, '
          f'{(df["status"]=="COMPLETE").sum()} COMPLETE, '
          f'{(df["status"]=="FAILED").sum()} FAILED')

    manifest_path = MANIFEST_S1 if stage == 1 else MANIFEST_S2

    while True:
        df = pd.read_csv(manifest_path, sep='\t')

        # Check running jobs for completion
        running = df[df['status'] == 'RUNNING']
        for idx, row in running.iterrows():
            done, n_bridge, pct = check_complete(row)
            if done:
                df.at[idx, 'status']       = 'COMPLETE'
                df.at[idx, 'bridge_poses'] = n_bridge
                df.at[idx, 'bridge_pct']   = pct
                df.at[idx, 'done_time']    = datetime.now().isoformat()
                print(f'  COMPLETE: {row["protac"]}/{row["gene"]} '
                      f'bridge={n_bridge} ({pct:.1f}%)')

        n_running = int((df['status'] == 'RUNNING').sum())
        n_slots   = MAX_PARALLEL - n_running

        # Submit pending jobs into empty slots
        pending = df[df['status'] == 'PENDING']
        for idx, row in pending.head(n_slots).iterrows():
            try:
                pid = submit_run(row, stage)
                df.at[idx, 'status']      = 'RUNNING'
                df.at[idx, 'submit_time'] = datetime.now().isoformat()
                print(f'  SUBMITTED: {row["protac"]}/{row["gene"]}  pid={pid}')
            except Exception as e:
                df.at[idx, 'status'] = 'FAILED'
                print(f'  FAILED submit: {row["protac"]}/{row["gene"]}: {e}')

        df.to_csv(manifest_path, sep='\t', index=False)

        n_remaining = int((df['status'].isin(['PENDING', 'RUNNING'])).sum())
        n_total     = len(df)
        n_done      = int((df['status'] == 'COMPLETE').sum())
        print(f'  [{datetime.now().strftime("%H:%M:%S")}] '
              f'{n_done}/{n_total} complete, {n_remaining} remaining')

        if n_remaining == 0:
            print('All done.')
            break
        time.sleep(60)


def cmd_filter_stage1(args):
    if not MANIFEST_S1.exists():
        print('No stage 1 manifest.')
        return
    df = pd.read_csv(MANIFEST_S1, sep='\t')
    complete = df[df['status'] == 'COMPLETE']
    bridging = complete[complete['bridge_poses'] >= STAGE2_MIN_BRIDGE_POSES]
    print(f'Stage 1 complete: {len(complete)}/{len(df)}')
    print(f'Passing bridge filter (>={STAGE2_MIN_BRIDGE_POSES} bridging poses): {len(bridging)}')
    print(f'\nBy PROTAC:')
    for protac in bridging['protac'].unique():
        n = (bridging['protac'] == protac).sum()
        print(f'  {protac:12s}  {n} hits')


def cmd_status(args):
    for stage, path in [(1, MANIFEST_S1), (2, MANIFEST_S2)]:
        if not path.exists():
            print(f'Stage {stage}: no manifest')
            continue
        df = pd.read_csv(path, sep='\t')
        print(f'Stage {stage}: {len(df)} runs — '
              + ', '.join(f'{s}={int((df["status"]==s).sum())}'
                          for s in ['PENDING','RUNNING','COMPLETE','FAILED']))


def main():
    p = argparse.ArgumentParser()
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--plan',          action='store_true')
    g.add_argument('--stage1',        action='store_true')
    g.add_argument('--stage2',        action='store_true')
    g.add_argument('--filter-stage1', action='store_true')
    g.add_argument('--status',        action='store_true')
    p.add_argument('--protac', default=None)
    p.add_argument('--gene',   default=None)
    args = p.parse_args()

    if args.plan:
        cmd_plan(args)
    elif args.stage1:
        cmd_stage(args, stage=1)
    elif args.stage2:
        cmd_stage(args, stage=2)
    elif args.filter_stage1:
        cmd_filter_stage1(args)
    elif args.status:
        cmd_status(args)


if __name__ == '__main__':
    main()
