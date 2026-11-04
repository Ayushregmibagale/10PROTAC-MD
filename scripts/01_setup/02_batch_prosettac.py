#!/usr/bin/env python3
"""
02_batch_prosettac.py — batch pRosettaC submission for all kinase × PROTAC pairs
==================================================================================
Two-stage screen:
  Stage 1 (screen):  N_DOCK=50  → ~200 Rosetta-refined poses per kinase
  Stage 2 (deep):    N_DOCK=1000, Full: True → ~20,000 poses for top bridging hits

Runs on Ernie via SSH subprocess.  Fully resumable — checks for score.sc to detect
completion.  Submits up to MAX_PARALLEL jobs simultaneously.

pRosettaC interface (Ernie):
  python3 RUNNER Protac_params.txt --n_dock N --n_workers W

Stage 2 adds "Full: True" and "Anchor atoms: 6 1" lines to Protac_params.txt
(written as Protac_params_full.txt alongside the original).

Prerequisites:
  - ~/10PROTAC-MD/prosettac_inputs/{PROTAC}/{GENE}_{UID}/Protac_params.txt
  - pRosettaC runner on Ernie at ERNIE_RUNNER
  - SSH key-based access to Ernie (host alias 'ernie')

Usage:
  python 02_batch_prosettac.py --plan                         # list all runs
  python 02_batch_prosettac.py --stage1                       # screen (N_DOCK=50)
  python 02_batch_prosettac.py --stage2                       # deep (N_DOCK=1000, Full)
  python 02_batch_prosettac.py --stage1 --protac SJF8240      # one PROTAC
  python 02_batch_prosettac.py --stage1 --gene DDR2           # one kinase
  python 02_batch_prosettac.py --status                       # progress report
  python 02_batch_prosettac.py --filter-stage1                # list stage2 candidates

State machine per run (stored in manifest):
  PENDING → SUBMITTED → RUNNING → COMPLETE | FAILED
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
MANIFEST_S1  = ROOT / 'runs' / 'manifest_stage1.tsv'
MANIFEST_S2  = ROOT / 'runs' / 'manifest_stage2.tsv'

# ── Ernie config ───────────────────────────────────────────────────────────────
ERNIE_HOST    = 'ernie'
ERNIE_ROOT    = '/home/s1638308/10PROTAC-MD'
ERNIE_PYTHON  = '/home/s1638308/miniconda3/envs/protac/bin/python3'
ERNIE_RUNNER  = '/home/s1638308/protac/scripts/02_rosetta_utils/03_run_local_prosettac.py'
MAX_PARALLEL  = 10    # concurrent kinases on Ernie (10 × 2 workers = 20 CPUs ~20%)
N_WORKERS_PER = 2     # parallel Rosetta workers per pRosettaC job

# ── pRosettaC params ───────────────────────────────────────────────────────────
STAGE1_NDOCK = 200    # 200 global PatchDock poses → screen
STAGE2_NDOCK = 1000   # + Full: True → ~20,000 poses per kinase
BRIDGE_DIST  = 5.0    # Å, used in bridging pose counting
VHL_MAX_RES  = 103    # last VHL residue in combined PDB
KIN_MIN_RES  = 200    # first kinase residue in combined PDB (renumbered)

# Stage 2 selection: at least 1 bridging pose out of Stage 1
STAGE2_MIN_BRIDGE_POSES = 1


# ── Helpers ────────────────────────────────────────────────────────────────────

def ernie_run(cmd: str, check=True) -> subprocess.CompletedProcess:
    full = ['ssh', ERNIE_HOST, cmd]
    return subprocess.run(full, capture_output=True, text=True, check=check)


def ernie_exists(path: str) -> bool:
    r = ernie_run(f'test -f {path}', check=False)
    return r.returncode == 0


def get_all_inputs(protac_filter=None, gene_filter=None):
    """Return list of dicts for all prepared pRosettaC input dirs (have Protac_params.txt)."""
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
            params_file = run_dir / 'Protac_params.txt'
            if not params_file.exists():
                continue
            # run_dir.name format: {GENE}_{UNIPROT_ID}
            parts   = run_dir.name.rsplit('_', 1)
            gene    = parts[0] if len(parts) == 2 else run_dir.name
            uid     = parts[1] if len(parts) == 2 else run_dir.name
            if gene_filter and gene != gene_filter:
                continue
            pairs.append({
                'protac':      protac,
                'gene':        gene,
                'uid':         uid,
                'run_name':    run_dir.name,
                'local_dir':   str(run_dir),
                'params_file': str(params_file),
            })
    return pairs


def build_manifest(stage: int, protac_filter=None, gene_filter=None) -> pd.DataFrame:
    manifest_path = MANIFEST_S1 if stage == 1 else MANIFEST_S2
    pairs = get_all_inputs(protac_filter, gene_filter)

    if stage == 2:
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
    if df.empty:
        print('No input pairs found.')
        sys.exit(1)

    df['status']        = 'PENDING'
    df['n_dock']        = STAGE1_NDOCK if stage == 1 else STAGE2_NDOCK
    df['bridge_poses']  = -1
    df['bridge_pct']    = -1.0
    df['submit_time']   = ''
    df['done_time']     = ''
    df['ernie_outdir']  = df.apply(
        lambda r: f'{ERNIE_ROOT}/runs/stage{stage}/{r["protac"]}/{r["run_name"]}', axis=1
    )

    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    if manifest_path.exists():
        existing = pd.read_csv(manifest_path, sep='\t')
        key = ['protac', 'gene', 'uid']
        df = df.set_index(key)
        existing = existing.set_index(key)
        for col in ['status', 'bridge_poses', 'bridge_pct', 'submit_time', 'done_time']:
            if col in existing.columns:
                df[col] = existing[col].reindex(df.index).fillna(df[col])
        df = df.reset_index()
    df.to_csv(manifest_path, sep='\t', index=False)
    return df


def _make_full_params(local_dir: Path) -> Path:
    """Write Protac_params_full.txt alongside Protac_params.txt with Full:True added."""
    src = local_dir / 'Protac_params.txt'
    dst = local_dir / 'Protac_params_full.txt'
    text = src.read_text()
    # Remove any existing Full: or Anchor atoms: lines then append
    lines = [l for l in text.splitlines()
             if not l.startswith('Full:') and not l.startswith('Anchor atoms:')]
    lines.append('Full: True')
    lines.append('Anchor atoms: 6 1')
    dst.write_text('\n'.join(lines) + '\n')
    return dst


def submit_run(row: pd.Series, stage: int) -> str:
    """Rsync input dir to Ernie and launch pRosettaC. Returns pid string."""
    local_dir  = Path(row['local_dir'])
    out_dir    = row['ernie_outdir']
    n_dock     = row['n_dock']

    # For stage 2, inject Full: True into a separate params file
    if stage == 2:
        params_file = _make_full_params(local_dir)
        ernie_params = f'{out_dir}/Protac_params_full.txt'
    else:
        params_file = local_dir / 'Protac_params.txt'
        ernie_params = f'{out_dir}/Protac_params.txt'

    # Create remote dir and rsync input files
    ernie_run(f'mkdir -p {out_dir}')
    subprocess.run(
        ['rsync', '-az', str(local_dir) + '/', f'{ERNIE_HOST}:{out_dir}/'],
        check=True, capture_output=True
    )

    run_cmd = (
        f'export PATH="$HOME/miniconda3/envs/protac/bin:$PATH"; '
        f'export OMP_NUM_THREADS=1; '
        f'cd {out_dir} && '
        f'nohup {ERNIE_PYTHON} {ERNIE_RUNNER} '
        f'{ernie_params} '
        f'--n_dock {n_dock} --n_workers {N_WORKERS_PER} '
        f'> {out_dir}/pRosettaC.log 2>&1 & echo $!'
    )
    r = ernie_run(run_cmd)
    return r.stdout.strip()


def check_complete(row: pd.Series) -> tuple:
    """Return (is_complete, bridge_poses, bridge_pct)."""
    score_sc = f'{row["ernie_outdir"]}/score.sc'
    if not ernie_exists(score_sc):
        return False, -1, -1.0
    local_score = RUNS_DIR / 'stage_scores' / row['protac'] / f'{row["run_name"]}_score.sc'
    local_score.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(['rsync', '-az', f'{ERNIE_HOST}:{score_sc}', str(local_score)],
                   capture_output=True)
    bridge_poses, bridge_pct = _count_bridging(local_score, row)
    return True, bridge_poses, bridge_pct


def _count_bridging(score_sc: Path, row: pd.Series):
    if not score_sc.exists():
        return 0, 0.0
    try:
        df = pd.read_csv(score_sc, sep=r'\s+', comment='#', header=0)
        bridge_col = [c for c in df.columns if 'bridg' in c.lower()]
        if bridge_col:
            n   = int((df[bridge_col[0]] > 0).sum())
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
    # ~32 min per kinase at N_DOCK=200 (4× the kinome_scan 50-pose ~8 min estimate)
    for k, v in sorted(by_protac.items()):
        s1_hrs = v * 32 / 60
        print(f'  {k:12s}  {v:4d} pairs  ~{s1_hrs:.0f} CPU-hrs (stage1, N_DOCK={STAGE1_NDOCK})')
    total = sum(by_protac.values())
    print(f'\nTotal stage1 wall-time at MAX_PARALLEL={MAX_PARALLEL}: '
          f'~{total * 32 / 60 / MAX_PARALLEL:.0f} hrs')


def cmd_stage(args, stage: int):
    df = build_manifest(stage, args.protac, args.gene)
    pending = df[df['status'] == 'PENDING']
    print(f'Stage {stage}: {len(pending)} PENDING, '
          f'{(df["status"]=="COMPLETE").sum()} COMPLETE, '
          f'{(df["status"]=="FAILED").sum()} FAILED')

    manifest_path = MANIFEST_S1 if stage == 1 else MANIFEST_S2

    while True:
        df = pd.read_csv(manifest_path, sep='\t')

        running = df[df['status'] == 'RUNNING']
        for idx, row in running.iterrows():
            done, n_bridge, pct = check_complete(row)
            if done:
                df.at[idx, 'status']        = 'COMPLETE'
                df.at[idx, 'bridge_poses']  = n_bridge
                df.at[idx, 'bridge_pct']    = pct
                df.at[idx, 'done_time']     = datetime.now().isoformat()
                print(f'  COMPLETE: {row["protac"]}/{row["gene"]} '
                      f'bridge={n_bridge} ({pct:.1f}%)')

        n_running = int((df['status'] == 'RUNNING').sum())
        n_slots   = MAX_PARALLEL - n_running

        pending = df[df['status'] == 'PENDING']
        for idx, row in pending.head(n_slots).iterrows():
            try:
                pid = submit_run(row, stage)
                df.at[idx, 'status']       = 'RUNNING'
                df.at[idx, 'submit_time']  = datetime.now().isoformat()
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
    for protac in sorted(bridging['protac'].unique()):
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
                          for s in ['PENDING', 'RUNNING', 'COMPLETE', 'FAILED']))


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
