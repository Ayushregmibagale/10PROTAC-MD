#!/usr/bin/env python3
"""
01_place_warhead_batch.py — automated warhead placement for all kinase × PROTAC pairs
=======================================================================================
For each (PROTAC, kinase structure), superimposes the kinase onto the PROTAC's
reference co-crystal structure via Cα MCS alignment, then transfers the warhead
coordinates into the kinase frame. Output is a combined PDB (VHL + kinase + warhead)
ready to feed into pRosettaC.

Strategy per warhead type:
  foretinib (SJF8240, SJFδ):  MCS vs foretinib in 3LQ8 (MET) → existing approach
  lapatinib (SJF1521):         MCS vs lapatinib in 1XKK (EGFR) → Type I alignment
  VS-4718   (GSK215):          MCS vs VS-4718 in 5JC7 (FAK)    → FERM domain kinase
  JH-I-25   (JNJ-1013):        MCS vs JH-I-25 in IRAK1 co-crystal
  MK-2206   (MS143):           SKIP — allosteric AKT, non-ATP-site, not portable

Prerequisites:
  - ~/10PROTAC-MD/structures/{GENE}/{KLIFS_ID}.pdb  (from 00_download_klifs_structures.py)
  - ~/10PROTAC-MD/warheads/{PROTAC}_warhead.sdf     (placed manually once)
  - ~/10PROTAC-MD/reference_structures/{PROTAC}_ref.pdb  (reference co-crystal)
  - RDKit, numpy, pandas (conda: drMD or ml_env)

Usage:
  python 01_place_warhead_batch.py --plan                   # list all pairs, check refs
  python 01_place_warhead_batch.py --run                    # place all (resumable)
  python 01_place_warhead_batch.py --run --protac SJF8240   # one PROTAC only
  python 01_place_warhead_batch.py --run --gene DDR2        # one kinase only
  python 01_place_warhead_batch.py --status                 # show progress

Output:
  ~/10PROTAC-MD/prosettac_inputs/{PROTAC}/{GENE}_{KLIFS_ID}/
      kinase.pdb        — cleaned kinase chain
      warhead.sdf       — warhead placed in kinase frame
      complex_input.pdb — VHL + kinase + warhead (pRosettaC input)
"""

import argparse
import sys
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path

ROOT           = Path(__file__).resolve().parents[2]
STRUCT_DIR     = ROOT / 'structures'
WARHEAD_DIR    = ROOT / 'warheads'
REF_DIR        = ROOT / 'reference_structures'
INPUT_DIR      = ROOT / 'prosettac_inputs'
MANIFEST_TSV   = ROOT / 'structures' / 'klifs_manifest.tsv'
PLACEMENT_TSV  = ROOT / 'prosettac_inputs' / 'placement_manifest.tsv'

VHL_PDB        = Path.home() / 'protac-md' / 'reference_structures' / 'vhl_elongin.pdb'

# ── PROTAC config ──────────────────────────────────────────────────────────────
PROTAC_CONFIG = {
    'SJF8240': {
        'warhead':    'foretinib',
        'ref_pdb':    '3LQ8',
        'dfg_states': ['DFG-out', 'DFG-out-like'],
        'skip':       False,
    },
    'SJF_delta': {
        'warhead':    'foretinib',
        'ref_pdb':    '3LQ8',
        'dfg_states': ['DFG-out', 'DFG-out-like'],
        'skip':       False,
    },
    'GSK215': {
        'warhead':    'VS4718',
        'ref_pdb':    '5JC7',
        'dfg_states': ['DFG-in', 'DFG-out', 'DFG-out-like', 'DFG-inter'],
        'skip':       False,
    },
    'SJF1521': {
        'warhead':    'lapatinib',
        'ref_pdb':    '1XKK',
        'dfg_states': ['DFG-in'],
        'skip':       False,
    },
    'JNJ1013': {
        'warhead':    'JHI25',
        'ref_pdb':    None,          # need IRAK1 co-crystal — fill in when found
        'dfg_states': ['DFG-in'],
        'skip':       True,          # skip until reference PDB confirmed
        'skip_reason': 'No IRAK1 co-crystal with JH-I-25 confirmed yet',
    },
    'MS143': {
        'warhead':    'MK2206',
        'ref_pdb':    None,
        'dfg_states': [],
        'skip':       True,
        'skip_reason': 'Allosteric AKT PH-domain binding — not portable to other kinases',
    },
}


def kabsch_align(mobile: np.ndarray, ref: np.ndarray):
    """Return rotation R and translation t such that mobile @ R.T + t ≈ ref."""
    mob_c = mobile.mean(axis=0)
    ref_c = ref.mean(axis=0)
    mob_centered = mobile - mob_c
    ref_centered = ref - ref_c
    H = mob_centered.T @ ref_centered
    U, _, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = ref_c - mob_c @ R.T
    return R, t


def parse_ca_coords(pdb_path: Path):
    """Extract {residue_number: CA_xyz} from a PDB file."""
    coords = {}
    with open(pdb_path) as f:
        for line in f:
            if line[:4] in ('ATOM', 'HETA') and line[12:16].strip() == 'CA':
                rnum = int(line[22:26])
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                coords[rnum] = np.array([x, y, z])
    return coords


def place_warhead(kinase_pdb: Path, protac: str, out_dir: Path) -> bool:
    """
    Place warhead into kinase frame by Cα superposition onto reference co-crystal.
    Writes kinase.pdb, warhead.sdf, complex_input.pdb into out_dir.
    Returns True on success.
    """
    cfg = PROTAC_CONFIG[protac]
    ref_pdb_path = REF_DIR / f'{cfg["ref_pdb"]}_kinase.pdb'
    warhead_sdf  = WARHEAD_DIR / f'{cfg["warhead"]}_warhead.sdf'

    for required in [ref_pdb_path, warhead_sdf, VHL_PDB]:
        if not required.exists():
            print(f'  MISSING: {required}')
            return False

    out_dir.mkdir(parents=True, exist_ok=True)

    # Load Cα coords from kinase and reference
    kin_ca  = parse_ca_coords(kinase_pdb)
    ref_ca  = parse_ca_coords(ref_pdb_path)

    # Common residues (pocket-aligned by KLIFS — use intersection)
    common = sorted(set(kin_ca) & set(ref_ca))
    if len(common) < 20:
        print(f'  Too few shared CA ({len(common)}) for {kinase_pdb.name}')
        return False

    mob_pts = np.array([kin_ca[r] for r in common])
    ref_pts = np.array([ref_ca[r] for r in common])
    R, t = kabsch_align(mob_pts, ref_pts)

    # Apply inverse transform (ref→kinase frame) to warhead coordinates via RDKit
    import_ok = _transform_warhead_rdkit(warhead_sdf, R, t, out_dir / 'warhead.sdf',
                                         direction='ref_to_kinase')
    if not import_ok:
        return False

    # Write kinase.pdb (copy with chain reassigned to D)
    _write_kinase_chain(kinase_pdb, out_dir / 'kinase.pdb', chain='D', resnum_offset=200)

    # Assemble complex_input.pdb = VHL + kinase + warhead (as HETATM chain X)
    _assemble_complex(VHL_PDB, out_dir / 'kinase.pdb', out_dir / 'warhead.sdf',
                      out_dir / 'complex_input.pdb')
    return True


def _transform_warhead_rdkit(sdf_in: Path, R: np.ndarray, t: np.ndarray,
                              sdf_out: Path, direction: str) -> bool:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromMolFile(str(sdf_in), removeHs=False)
        if mol is None:
            return False
        conf = mol.GetConformer()
        pos = conf.GetPositions()
        if direction == 'ref_to_kinase':
            # inverse: kinase = ref @ R.T + t  →  ref = (kinase - t) @ R
            new_pos = (pos - t) @ R
        else:
            new_pos = pos @ R.T + t
        for i, xyz in enumerate(new_pos):
            conf.SetAtomPosition(i, xyz)
        writer = Chem.SDWriter(str(sdf_out))
        writer.write(mol)
        writer.close()
        return True
    except Exception as e:
        print(f'  RDKit transform failed: {e}')
        return False


def _write_kinase_chain(pdb_in: Path, pdb_out: Path, chain: str, resnum_offset: int):
    lines = []
    with open(pdb_in) as f:
        for line in f:
            if line[:4] in ('ATOM', 'HETA'):
                rnum = int(line[22:26]) + resnum_offset
                line = line[:21] + chain + f'{rnum:4d}' + line[26:]
            lines.append(line)
    with open(pdb_out, 'w') as f:
        f.writelines(lines)


def _assemble_complex(vhl_pdb: Path, kin_pdb: Path, warhead_sdf: Path, out_pdb: Path):
    from rdkit import Chem
    lines = []
    with open(vhl_pdb) as f:
        for l in f:
            if l.startswith(('ATOM', 'HETA', 'TER')):
                lines.append(l)
    lines.append('TER\n')
    with open(kin_pdb) as f:
        for l in f:
            if l.startswith(('ATOM', 'HETA', 'TER')):
                lines.append(l)
    lines.append('TER\n')
    # warhead from SDF as HETATM chain X
    mol = Chem.MolFromMolFile(str(warhead_sdf), removeHs=False)
    if mol:
        conf = mol.GetConformer()
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            sym = atom.GetSymbol()
            lines.append(
                f'HETATM{i+1:5d}  {sym:<3s} LIG X   1    '
                f'{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}  1.00  0.00          {sym:>2s}\n'
            )
    lines.append('END\n')
    with open(out_pdb, 'w') as f:
        f.writelines(lines)


# ── Commands ───────────────────────────────────────────────────────────────────

def build_pairs(protac_filter=None, gene_filter=None):
    """Return list of (protac, gene, klifs_id, kinase_pdb_path) tuples."""
    if not MANIFEST_TSV.exists():
        print('No structure manifest — run 00_download_klifs_structures.py --query first')
        sys.exit(1)
    df = pd.read_csv(MANIFEST_TSV, sep='\t')
    pairs = []
    for protac, cfg in PROTAC_CONFIG.items():
        if cfg['skip']:
            continue
        if protac_filter and protac != protac_filter:
            continue
        subset = df[df['dfg_state'].isin(cfg['dfg_states'])].copy()
        if gene_filter:
            subset = subset[subset['gene'] == gene_filter]
        for _, row in subset.iterrows():
            gene     = str(row['gene'])
            klifs_id = int(row['klifs_id'])
            pdb_path = STRUCT_DIR / gene / f'{klifs_id}_{str(row["dfg_state"]).replace(" ", "_")}.pdb'
            pairs.append((protac, gene, klifs_id, pdb_path))
    return pairs


def cmd_plan(args):
    pairs = build_pairs(args.protac, args.gene)
    print(f'Total pairs: {len(pairs)}')
    for protac in PROTAC_CONFIG:
        cfg = PROTAC_CONFIG[protac]
        n = sum(1 for p, *_ in pairs if p == protac)
        if cfg['skip']:
            print(f'  {protac:12s}  SKIP — {cfg.get("skip_reason", "")}')
        else:
            print(f'  {protac:12s}  {n} pairs')
    print('\nSkipped PROTACs:')
    for protac, cfg in PROTAC_CONFIG.items():
        if cfg['skip']:
            print(f'  {protac}: {cfg.get("skip_reason", "no reason given")}')


def cmd_run(args):
    pairs = build_pairs(args.protac, args.gene)
    print(f'Running placement for {len(pairs)} pairs ...')
    ok = fail = skip = 0
    for protac, gene, klifs_id, pdb_path in pairs:
        out_dir = INPUT_DIR / protac / f'{gene}_{klifs_id}'
        done_flag = out_dir / 'complex_input.pdb'
        if done_flag.exists():
            skip += 1
            continue
        if not pdb_path.exists():
            print(f'  MISSING PDB: {pdb_path}')
            fail += 1
            continue
        success = place_warhead(pdb_path, protac, out_dir)
        if success:
            ok += 1
        else:
            fail += 1
    print(f'Done — placed: {ok}, skipped (exists): {skip}, failed: {fail}')


def cmd_status(args):
    if not MANIFEST_TSV.exists():
        print('No manifest.')
        return
    df = pd.read_csv(MANIFEST_TSV, sep='\t')
    for protac in PROTAC_CONFIG:
        d = INPUT_DIR / protac
        n_done = len(list(d.glob('*/complex_input.pdb'))) if d.exists() else 0
        cfg = PROTAC_CONFIG[protac]
        expected = len(df[df['dfg_state'].isin(cfg['dfg_states'])]) if not cfg['skip'] else 0
        tag = 'SKIP' if cfg['skip'] else f'{n_done}/{expected}'
        print(f'  {protac:12s}  {tag}')


def main():
    p = argparse.ArgumentParser()
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--plan',   action='store_true')
    g.add_argument('--run',    action='store_true')
    g.add_argument('--status', action='store_true')
    p.add_argument('--protac', default=None)
    p.add_argument('--gene',   default=None)
    args = p.parse_args()

    if args.plan:
        cmd_plan(args)
    elif args.run:
        cmd_run(args)
    elif args.status:
        cmd_status(args)


if __name__ == '__main__':
    main()
