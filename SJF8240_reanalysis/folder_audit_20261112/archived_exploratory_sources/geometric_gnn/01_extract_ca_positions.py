#!/usr/bin/env python3
"""
Extract mean Cα positions from pRosettaC rank-1 ternary complex PDB files.

For each POI, parses the complex_protein.pdb to get Cα (x, y, z) coordinates
for VHL residues (chain A, 1-103) and kinase residues (chain A/D, resid>103).

Outputs:
  geometric_gnn/positions/<POI>.npz
    - vhl_resids: (n_vhl,) int
    - vhl_pos:    (n_vhl, 3) float32  Angstroms
    - kin_resids: (n_kin,) int
    - kin_pos:    (n_kin, 3) float32  Angstroms
"""
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SIMUL = Path.home() / 'protac-md' / 'simulations'
RSCREEN = Path.home() / 'protac-md' / 'random_screen'
OUT_DIR = Path(__file__).resolve().parent / 'positions'

LABELS = {
    'DDR2':        1,
    'ABL1':        0,
    'AXL':         0,
    'EPHA2':       1,
    'EPHB2':       1,
    'MAPK14_full': 1,
    'MET':         1,
    'RIPK2_full':  1,
}

PDB_PATHS = {
    'DDR2':        SIMUL / 'DDR2'        / 'structure/rank1/complex_protein.pdb',
    'ABL1':        SIMUL / 'ABL1'        / 'structure/rank1/complex_protein.pdb',
    'AXL':         SIMUL / 'AXL'         / 'structure/rank1/complex_protein.pdb',
    'EPHA2':       SIMUL / 'EPHA2'       / 'structure/rank1/complex_protein.pdb',
    'EPHB2':       RSCREEN / 'EPHB2/pose001/structure/complex_protein.pdb',
    'MAPK14_full': SIMUL / 'MAPK14_full' / 'structure/rank1/complex_protein.pdb',
    'MET':         SIMUL / 'MET'         / 'structure/rank1/complex_protein.pdb',
    'RIPK2_full':  SIMUL / 'RIPK2_full'  / 'structure/rank1/complex_protein.pdb',
}


def parse_ca_positions(pdb_path: Path) -> dict:
    """
    Parse Cα positions from PDB. Returns dict of residue_num → (x, y, z).
    VHL = chain A, resid 1-103.
    Kinase = chain A residues >103 OR chain D (EPHB2 convention).
    """
    vhl = {}   # resid → (x, y, z)
    kin = {}

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            atom_name = line[12:16].strip()
            if atom_name != 'CA':
                continue
            chain = line[21]
            try:
                resid = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            if chain == 'A':
                if 1 <= resid <= 103:
                    vhl[resid] = (x, y, z)
                elif resid > 103:
                    kin[resid] = (x, y, z)
            elif chain in ('D', 'B') and chain != 'B':
                # chain D = kinase in EPHB2 pRosettaC convention
                kin[resid] = (x, y, z)
            elif chain == 'D':
                kin[resid] = (x, y, z)

    return vhl, kin


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"{'POI':<16} {'label':<6} {'n_vhl':>6} {'n_kin':>6}  pdb")
    print('-' * 70)

    for poi, pdb_path in PDB_PATHS.items():
        if not pdb_path.exists():
            print(f"{poi:<16} MISSING: {pdb_path}")
            continue

        vhl_dict, kin_dict = parse_ca_positions(pdb_path)

        vhl_resids = np.array(sorted(vhl_dict.keys()), dtype=np.int32)
        kin_resids = np.array(sorted(kin_dict.keys()), dtype=np.int32)
        vhl_pos    = np.array([vhl_dict[r] for r in vhl_resids], dtype=np.float32)
        kin_pos    = np.array([kin_dict[r] for r in kin_resids], dtype=np.float32)

        out_path = OUT_DIR / f'{poi}.npz'
        np.savez(out_path,
                 vhl_resids=vhl_resids, vhl_pos=vhl_pos,
                 kin_resids=kin_resids, kin_pos=kin_pos,
                 label=np.array([LABELS[poi]], dtype=np.int64))

        print(f"{poi:<16} {LABELS[poi]:<6} {len(vhl_resids):>6} {len(kin_resids):>6}  {pdb_path.name}")

    print(f"\nSaved {len(PDB_PATHS)} position files to {OUT_DIR}")


if __name__ == '__main__':
    main()
