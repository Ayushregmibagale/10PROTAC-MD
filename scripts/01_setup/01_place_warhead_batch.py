#!/usr/bin/env python3
"""
01_place_warhead_batch.py — automated warhead placement for all kinase × PROTAC pairs
=======================================================================================
Places each PROTAC's warhead into every kinase structure frame using Cα superposition
onto a reference co-crystal, then transfers warhead coordinates.

Foretinib PROTACs (SJF8240, SJFδ) are already staged from kinome_scan via
00b_build_merged_manifest.py — this script handles the remaining 3 PROTACs only.

Strategy per warhead:
  lapatinib (SJF1521):  Cα superposition vs EGFR in 1XKK (Type I, DFG-in)
  VS-4718   (GSK215):   Cα superposition vs FAK  in 5JC7 (FERM-domain kinase)
  JH-I-25   (JNJ-1013): Cα superposition vs IRAK1 in 6EGA (Type I, DFG-in)

Input structure selection from merged_manifest.tsv:
  SJF1521 / JNJ-1013: KLIFS DFG-in first, AF2 fallback
  GSK215: KLIFS any DFG, AF2 fallback

Prerequisites:
  ~/10PROTAC-MD/warheads/{lapatinib,VS4718,JHI25}_warhead.sdf
  ~/10PROTAC-MD/reference_structures/{1XKK,5JC7,6EGA}_kinase.pdb
  ~/10PROTAC-MD/structures/merged_manifest.tsv
  Conda env: ~/.conda/envs/klifs_screen (has rdkit via conda-forge opencadd)

Usage:
  python 01_place_warhead_batch.py --plan                   # count pairs, check prereqs
  python 01_place_warhead_batch.py --run                    # place all (resumable)
  python 01_place_warhead_batch.py --run --protac SJF1521   # one PROTAC
  python 01_place_warhead_batch.py --run --gene EGFR        # one gene
  python 01_place_warhead_batch.py --status

Output per pair:
  ~/10PROTAC-MD/prosettac_inputs/{PROTAC}/{GENE}_{UID}/
      kinase.pdb        — kinase chain (chain D, residues 200+)
      warhead.sdf       — warhead in kinase frame
      Protac_params.txt — pRosettaC config (ernie paths)
"""

import argparse, sys, re
import numpy as np
import pandas as pd
from pathlib import Path

ROOT        = Path(__file__).resolve().parents[2]
WARHEAD_DIR = ROOT / 'warheads'
REF_DIR     = ROOT / 'reference_structures'
INPUT_DIR   = ROOT / 'prosettac_inputs'
MANIFEST    = ROOT / 'structures' / 'merged_manifest.tsv'

VHL_PDB     = Path.home() / 'protac' / 'data' / 'structures' / 'VHL' / '4W9H_VBC_clean.pdb'
VHL_WH_SDF  = Path.home() / 'protac' / 'data' / 'structures' / 'VHL' / 'SJF8240_vhll_placed.sdf'
ERNIE_USER  = 's1638308'
ERNIE_HOME  = '/home/s1638308'

PROTAC_CONFIG = {
    'SJF8240':   {'skip': True,  'skip_reason': 'Already staged from kinome_scan'},
    'SJF_delta': {'skip': True,  'skip_reason': 'Already staged from kinome_scan'},
    'SJF1521': {
        'skip':        False,
        'warhead_sdf': WARHEAD_DIR / 'lapatinib_warhead.sdf',
        'ref_kinase':  REF_DIR / '1XKK_kinase.pdb',
        'struct_col':  'SJF1521_path',
        'src_col':     'SJF1521_source',
    },
    'GSK215': {
        'skip':        False,
        'warhead_sdf': WARHEAD_DIR / 'VS4718_warhead.sdf',
        'ref_kinase':  REF_DIR / '5JC7_kinase.pdb',
        'struct_col':  'GSK215_path',
        'src_col':     'GSK215_source',
    },
    'JNJ1013': {
        'skip':        False,
        'warhead_sdf': WARHEAD_DIR / 'JHI25_warhead.sdf',
        'ref_kinase':  REF_DIR / '6EGA_kinase.pdb',
        'struct_col':  'JNJ1013_path',
        'src_col':     'JNJ1013_source',
    },
    'MS143': {'skip': True, 'skip_reason': 'Allosteric AKT — excluded from kinome screen'},
}


# ── Kabsch alignment ──────────────────────────────────────────────────────────

def parse_ca(pdb_path: Path) -> np.ndarray:
    """Return (N,3) Cα coordinate array in residue order."""
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if line[:4] == 'ATOM' and line[12:16].strip() == 'CA':
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(coords) if coords else np.zeros((0, 3))


def kabsch(mobile: np.ndarray, ref: np.ndarray):
    """Return (R, t) minimising RMSD: mobile @ R.T + t ≈ ref."""
    mc, rc = mobile.mean(0), ref.mean(0)
    H  = (mobile - mc).T @ (ref - rc)
    U, _, Vt = np.linalg.svd(H)
    d  = np.linalg.det(Vt.T @ U.T)
    R  = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, rc - mc @ R.T


def common_ca_align(kin_pdb: Path, ref_pdb: Path):
    """
    Align up to min(N_kin, N_ref) Cα atoms (N-terminal truncation to common length).
    Returns (R, t) or (None, None) if too few atoms.
    """
    kin_ca = parse_ca(kin_pdb)
    ref_ca = parse_ca(ref_pdb)
    n = min(len(kin_ca), len(ref_ca))
    if n < 30:
        return None, None
    R, t = kabsch(kin_ca[:n], ref_ca[:n])
    rmsd = float(np.sqrt(((kin_ca[:n] @ R.T + t - ref_ca[:n]) ** 2).mean()))
    return R, t


# ── Warhead transform ─────────────────────────────────────────────────────────

def transform_warhead(sdf_in: Path, R: np.ndarray, t: np.ndarray, sdf_out: Path) -> bool:
    """Apply inverse rotation (ref→kinase frame) to warhead atoms via RDKit."""
    try:
        from rdkit import Chem
        mol = Chem.MolFromMolFile(str(sdf_in), removeHs=False)
        if mol is None:
            return False
        conf = mol.GetConformer()
        pos  = conf.GetPositions()
        # ref frame → kinase frame: inverse of (mob @ R.T + t = ref)
        new_pos = (pos - t) @ R
        for i, xyz in enumerate(new_pos):
            conf.SetAtomPosition(i, xyz.tolist())
        w = Chem.SDWriter(str(sdf_out))
        w.write(mol); w.close()
        return True
    except Exception as e:
        print(f'    RDKit error: {e}')
        return False


# ── Kinase chain writer ───────────────────────────────────────────────────────

def write_kinase_chain(src: Path, dst: Path, chain='D', offset=200):
    """Write kinase PDB with chain relabelled to D and residues offset to 200+."""
    lines = []
    seen_ter = False
    with open(src) as f:
        for line in f:
            if line.startswith('TER') and not seen_ter:
                seen_ter = True
                lines.append('TER\n')
            elif line[:4] in ('ATOM', 'HETA') and not seen_ter:
                try:
                    rnum = int(line[22:26]) + offset
                except ValueError:
                    continue
                line = line[:21] + chain + f'{rnum:4d}' + line[26:]
                lines.append(line)
    with open(dst, 'w') as f:
        f.writelines(lines)


# ── Protac_params writer ──────────────────────────────────────────────────────

def write_protac_params(out_dir: Path, protac: str, uid: str):
    ernie_dir  = f'{ERNIE_HOME}/10PROTAC-MD/prosettac_inputs/{protac}/{out_dir.name}'
    # VHL SMILES from SJF8240 — same VHL ligand for all PROTACs
    smi_path   = f'{ERNIE_HOME}/protac/data/structures/PROTAC/SJF8240/SJF8240.smi'
    params = (
        f'Structures: {ERNIE_HOME}/protac/data/structures/VHL/4W9H_VBC_clean.pdb '
        f'{ernie_dir}/protein.pdb\n'
        f'Chains: ABC A\n'
        f'Heads: {ERNIE_HOME}/protac/data/structures/VHL/SJF8240_vhll_placed.sdf '
        f'{ernie_dir}/warhead.sdf\n'
        f'Smiles: {smi_path}\n'
        f'OutDir: {ernie_dir}\n'
    )
    (out_dir / 'Protac_params.txt').write_text(params)


# ── Placement core ────────────────────────────────────────────────────────────

def place_one(gene: str, uid: str, struct_path: str, protac: str) -> bool:
    cfg     = PROTAC_CONFIG[protac]
    out_dir = INPUT_DIR / protac / f'{gene}_{uid}'
    done    = out_dir / 'warhead.sdf'
    if done.exists():
        return True

    kin_pdb = Path(struct_path)
    if not kin_pdb.exists():
        print(f'    Missing structure: {kin_pdb}')
        return False

    R, t = common_ca_align(kin_pdb, cfg['ref_kinase'])
    if R is None:
        print(f'    Too few Cα atoms for alignment: {gene}')
        return False

    out_dir.mkdir(parents=True, exist_ok=True)

    if not transform_warhead(cfg['warhead_sdf'], R, t, out_dir / 'warhead.sdf'):
        return False

    write_kinase_chain(kin_pdb, out_dir / 'protein.pdb')
    write_protac_params(out_dir, protac, uid)
    return True


# ── Commands ─────────────────────────────────────────────────────────────────

def load_manifest():
    if not MANIFEST.exists():
        print('No merged manifest — run 00b_build_merged_manifest.py first.')
        sys.exit(1)
    return pd.read_csv(MANIFEST, sep='\t')


def cmd_plan(args):
    df = load_manifest()
    print('Warhead placement plan:')
    for protac, cfg in PROTAC_CONFIG.items():
        if cfg.get('skip'):
            col = f'{protac}_path'
            if col in df.columns:
                n = df[col].notna().sum()
                print(f'  {protac:12s}  SKIP ({cfg["skip_reason"]}) — {n} already staged')
            else:
                print(f'  {protac:12s}  SKIP ({cfg["skip_reason"]})')
            continue
        col = cfg['struct_col']
        if col not in df.columns:
            print(f'  {protac:12s}  no column in manifest')
            continue
        n_total = df[col].notna().sum()
        n_done  = sum(1 for _, row in df.iterrows()
                      if pd.notna(row.get(col)) and
                      (INPUT_DIR / protac / f'{row.get("hgnc_gene","?")}_{row.get("uniprot_id","?")}' / 'warhead.sdf').exists())
        print(f'  {protac:12s}  {n_total} pairs  ({n_done} already done)')
        # check prereqs
        for req in [cfg['warhead_sdf'], cfg['ref_kinase']]:
            if not req.exists():
                print(f'    MISSING: {req}')


def cmd_run(args):
    df = load_manifest()
    for protac, cfg in PROTAC_CONFIG.items():
        if cfg.get('skip'):
            continue
        if args.protac and protac != args.protac:
            continue
        col = cfg.get('struct_col')
        if not col or col not in df.columns:
            continue
        subset = df[df[col].notna()].copy()
        if args.gene:
            subset = subset[subset['hgnc_gene'] == args.gene]
        print(f'\n{protac}: {len(subset)} kinases to place ...')
        ok = skip = fail = 0
        for _, row in subset.iterrows():
            gene = str(row.get('hgnc_gene', row.get('uniprot_id', '?')))
            uid  = str(row.get('uniprot_id', '?'))
            path = str(row[col])
            result = place_one(gene, uid, path, protac)
            if result is True:
                # check if it was already done or just placed
                out = INPUT_DIR / protac / f'{gene}_{uid}' / 'warhead.sdf'
                ok += 1
            else:
                fail += 1
        print(f'  placed: {ok}, failed: {fail}')


def cmd_status(args):
    df = load_manifest() if MANIFEST.exists() else None
    print('Placement status:')
    for protac in PROTAC_CONFIG:
        d = INPUT_DIR / protac
        n = sum(1 for x in d.glob('*/warhead.sdf')) if d.exists() else 0
        cfg = PROTAC_CONFIG[protac]
        if cfg.get('skip'):
            print(f'  {protac:12s}  {n} staged (pre-existing kinome_scan)')
        else:
            col = cfg.get('struct_col', '')
            expected = df[col].notna().sum() if df is not None and col in df.columns else '?'
            print(f'  {protac:12s}  {n}/{expected} done')


def main():
    p = argparse.ArgumentParser()
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--plan',   action='store_true')
    g.add_argument('--run',    action='store_true')
    g.add_argument('--status', action='store_true')
    p.add_argument('--protac', default=None)
    p.add_argument('--gene',   default=None)
    args = p.parse_args()
    if args.plan:   cmd_plan(args)
    elif args.run:  cmd_run(args)
    elif args.status: cmd_status(args)


if __name__ == '__main__':
    main()
