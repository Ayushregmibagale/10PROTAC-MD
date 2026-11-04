#!/usr/bin/env python3
"""
00_download_klifs_structures.py — bulk download human kinase structures from KLIFS
===================================================================================
Uses opencadd-klifs session API to fetch all human kinase structures, filter by
DFG conformation state, deduplicate to best structure per (gene, DFG state), and
download PDB files. Fully resumable.

Usage:
  python 00_download_klifs_structures.py --query              # fetch manifest only
  python 00_download_klifs_structures.py --download           # download all PDBs
  python 00_download_klifs_structures.py --download --dfg out # DFG-out only
  python 00_download_klifs_structures.py --download --gene EGFR
  python 00_download_klifs_structures.py --status

KLIFS DFG values: 'in', 'out', 'out-like', 'inter', 'na'

Env: ~/.conda/envs/klifs_screen  (conda-forge opencadd=1.0.1)
"""

import argparse
import sys
import time
import tempfile
import json
import hashlib
import shutil
import pandas as pd
from pathlib import Path
from tqdm import tqdm

ROOT         = Path(__file__).resolve().parents[2]
STRUCT_DIR   = ROOT / 'structures'
MANIFEST_TSV = STRUCT_DIR / 'klifs_manifest.tsv'

# DFG states each warhead type needs (KLIFS notation: 'in','out','out-like','inter')
DFG_FOR_WARHEAD = {
    'foretinib': ['out', 'out-like'],
    'lapatinib': ['in'],
    'VS4718':    ['in', 'out', 'out-like', 'inter'],
    'JHI25':     ['in'],
    'MK2206':    ['in', 'out', 'out-like', 'inter'],
}

# All DFG states we care about across all PROTACs
ALL_DFG_NEEDED = sorted({s for states in DFG_FOR_WARHEAD.values() for s in states})


def get_session():
    from opencadd.databases.klifs import setup_remote
    return setup_remote()


def cmd_query(args):
    if MANIFEST_TSV.exists():
        raise FileExistsError('Preserving existing KLIFS manifest; use a new project snapshot for a new query')
    print('Connecting to KLIFS and fetching all structures ...')
    session = get_session()
    df = session.structures.all_structures()
    print(f'  Total structures in KLIFS: {len(df)}')

    # Filter to Human only
    df = df[df['species.klifs'] == 'Human'].copy()
    print(f'  Human structures: {len(df)}')

    # Keep only DFG states we need
    df = df[df['structure.dfg'].isin(ALL_DFG_NEEDED)].copy()
    print(f'  After DFG filter ({ALL_DFG_NEEDED}): {len(df)} structures')

    # Rename to friendlier columns
    df = df.rename(columns={
        'structure.klifs_id':        'klifs_id',
        'structure.pdb_id':          'pdb_id',
        'structure.chain':           'chain',
        'kinase.klifs_name':         'gene',
        'structure.dfg':             'dfg_state',
        'structure.resolution':      'resolution',
        'structure.missing_residues':'missing_residues',
        'structure.qualityscore':    'quality_score',
        'structure.curation_flag':   'curated',
    })

    keep = ['klifs_id', 'gene', 'pdb_id', 'chain', 'dfg_state',
            'resolution', 'missing_residues', 'quality_score', 'curated']
    df = df[[c for c in keep if c in df.columns]].copy()
    df['resolution'] = pd.to_numeric(df['resolution'], errors='coerce')

    STRUCT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(MANIFEST_TSV, sep='\t', index=False)
    print(f'\nManifest saved → {MANIFEST_TSV}')
    print(f'Unique kinases: {df["gene"].nunique()}')
    print(f'DFG breakdown:\n{df["dfg_state"].value_counts().to_string()}')


def cmd_download(args):
    if not MANIFEST_TSV.exists():
        print('No manifest — run --query first.')
        sys.exit(1)

    df = pd.read_csv(MANIFEST_TSV, sep='\t')

    # Optional filters
    if args.dfg:
        df = df[df['dfg_state'] == args.dfg]
        print(f'DFG filter "{args.dfg}": {len(df)} structures')
    if args.gene:
        df = df[df['gene'] == args.gene]
        print(f'Gene filter "{args.gene}": {len(df)} structures')

    # Deduplicate: best structure per (gene, dfg_state)
    # Prefer: curated=True > lowest missing_residues > lowest resolution
    df = df.sort_values(
        ['curated', 'missing_residues', 'resolution'],
        ascending=[False, True, True]
    )
    df_best = df.drop_duplicates(subset=['gene', 'dfg_state'], keep='first').reset_index(drop=True)
    print(f'After dedup (best per gene×DFG): {len(df_best)} structures, '
          f'{df_best["gene"].nunique()} kinases')

    session = get_session()
    ok = fail = skip = 0

    for _, row in tqdm(df_best.iterrows(), total=len(df_best), desc='Downloading'):
        gene     = str(row['gene'])
        klifs_id = int(row['klifs_id'])
        dfg      = str(row['dfg_state'])
        out_dir  = STRUCT_DIR / gene
        out_dir.mkdir(parents=True, exist_ok=True)
        out_file = out_dir / f'{klifs_id}_{dfg}.pdb'

        provenance=out_file.with_suffix('.provenance.json')
        if out_file.exists():
            if not provenance.exists(): raise ValueError('Existing unverified download: '+str(out_file))
            meta=json.loads(provenance.read_text())
            if meta['klifs_id']!=klifs_id or meta['sha256']!=hashlib.sha256(out_file.read_bytes()).hexdigest():
                raise ValueError('Downloaded input provenance mismatch: '+str(out_file))
            skip += 1
            continue

        try:
            # OpenCADD filenames vary; isolate one API request instead of guessing
            # the saved name or counting an absent file as a successful download.
            with tempfile.TemporaryDirectory(dir=out_dir,prefix='.download-') as temporary:
                session.coordinates.to_pdb(structure_klifs_id=klifs_id,output_path=Path(temporary),entity='complex')
                candidates=list(Path(temporary).rglob('*.pdb'))
                if len(candidates)!=1 or not any(l.startswith('ATOM  ') for l in candidates[0].read_text().splitlines()):
                    raise ValueError('API did not produce exactly one nonempty protein PDB')
                shutil.copy2(candidates[0],out_file)
            meta={'klifs_id':klifs_id,'gene':gene,'pdb_id':str(row['pdb_id']),
                  'chain':str(row['chain']),'dfg_state':dfg,
                  'sha256':hashlib.sha256(out_file.read_bytes()).hexdigest()}
            provenance.write_text(json.dumps(meta,indent=2)+'\n')
            ok += 1
            time.sleep(0.3)
        except Exception as e:
            print(f'  FAIL {gene} {klifs_id}: {e}')
            fail += 1

    print(f'\nDone — downloaded: {ok}, skipped (exists): {skip}, failed: {fail}')


def cmd_status(args):
    if not MANIFEST_TSV.exists():
        print('No manifest.')
        return
    df = pd.read_csv(MANIFEST_TSV, sep='\t')
    n_on_disk = sum(1 for g in STRUCT_DIR.iterdir()
                    if g.is_dir() for _ in g.glob('*.pdb'))
    print(f'Manifest : {len(df)} entries, {df["gene"].nunique()} genes')
    print(f'PDB files: {n_on_disk}')
    print(f'\nDFG breakdown:\n{df["dfg_state"].value_counts().to_string()}')
    print(f'\nTop 10 genes by structure count:\n'
          f'{df["gene"].value_counts().head(10).to_string()}')


def main():
    p = argparse.ArgumentParser()
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--query',    action='store_true', help='Fetch manifest (no download)')
    g.add_argument('--download', action='store_true', help='Download PDB files (resumable)')
    g.add_argument('--status',   action='store_true', help='Show progress')
    p.add_argument('--dfg',  default=None, help='DFG state filter: in, out, out-like, inter')
    p.add_argument('--gene', default=None, help='Download one kinase only')
    args = p.parse_args()

    if args.query:
        cmd_query(args)
    elif args.download:
        cmd_download(args)
    elif args.status:
        cmd_status(args)


if __name__ == '__main__':
    main()
