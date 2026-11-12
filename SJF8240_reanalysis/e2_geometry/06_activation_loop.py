#!/usr/bin/env python3
"""
Activation loop conformation analysis.

For each POI, identifies the DFG motif in the kinase sequence (from PDB),
defines the activation loop (DFG to APE motif), then measures:
  - Activation loop RMSD (bridging frames vs non-bridging frames)
  - DFG-in vs DFG-out classification per frame (distance from DFG-Phe to αC-helix Glu)
  - Whether bridging events preferentially occur in DFG-in or DFG-out state

Hypothesis: ABL1/AXL might adopt DFG-out or intermediate conformations during
bridging that position the activation loop to block E2 access.

Run: ~/.conda/envs/gromacs/bin/python3 -u 06_activation_loop.py
"""

# RETIRED_BY_FOLDER_AUDIT_20261112: original source is preserved in the audit archive.
if __name__ == "__main__":
    raise SystemExit("Retired exploratory entry point: see folder_audit_20261112/AUDIT.md and use RUN_ANALYSIS.py. Prior outputs are preserved but not certified by this script.")

import numpy as np
import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

try:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import capped_distance
except ImportError:
    print("[ERROR] MDAnalysis not found. Run: ~/.conda/envs/gromacs/bin/python3 ...")
    import sys; sys.exit(1)

SCRIPT_DIR  = Path(__file__).resolve().parent
RS_SIM_DIR  = Path.home() / 'protac-md' / 'random_screen'
RS_ANA_DIR  = Path.home() / 'protac-md' / 'analysis' / 'random_screen'
OUT_DIR     = SCRIPT_DIR / 'results'
OUT_DIR.mkdir(parents=True, exist_ok=True)

LABELS = {
    'DDR2':1,'ABL1':0,'AXL':0,'EPHA2':1,
    'EPHB2':1,'MAPK14_full':1,'MET':1,'RIPK2_full':1,
}
POI_ORDER  = ['DDR2','ABL1','AXL','EPHA2','EPHB2','MAPK14_full','MET','RIPK2_full']
STRIDE     = 10
VHL_CUT    = 5.0
KIN_CUT    = 5.0


def _detect_atom_ranges(gro_path):
    """Detect VHL and kinase atom ranges from GRO file."""
    lines = open(gro_path).readlines()
    n_atoms = int(lines[1])
    atoms = lines[2:2+n_atoms]
    boundaries = [0]
    prev_res = None
    for i, l in enumerate(atoms):
        try:
            resnum = int(l[0:5])
        except ValueError:
            continue
        if prev_res is not None and abs(resnum - prev_res) > 30:
            boundaries.append(i)
        prev_res = resnum
    boundaries.append(n_atoms)
    chains = []
    for j in range(len(boundaries)-1):
        s, e = boundaries[j], boundaries[j+1]
        chains.append({'start': s+1, 'end': e,
                       'first_res': int(atoms[s][0:5]),
                       'last_res':  int(atoms[e-1][0:5])})
    try:
        vhl = next(c for c in chains if c['first_res'] == 1)
        kin = next(c for c in chains if c['first_res'] > 103)
    except StopIteration:
        return None
    return (vhl['start'], vhl['end']), (kin['start'], kin['end'])


def find_dfg_in_pdb(pdb_path, kin_start_resid):
    """
    Scan PDB for D-F-G triplet in kinase chain (chain D or chain A resid>103).
    Returns resid of D in DFG motif (the Asp), or None.
    """
    residues = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            ch = line[21]
            try:
                resid = int(line[22:26])
                resname = line[17:20].strip()
            except ValueError:
                continue
            if ch == 'D' or (ch == 'A' and resid > 103):
                residues[resid] = resname

    sorted_res = sorted(residues.keys())
    for i in range(len(sorted_res)-2):
        r1, r2, r3 = sorted_res[i], sorted_res[i+1], sorted_res[i+2]
        if (r2 == r1+1 and r3 == r2+1 and
            residues[r1] == 'ASP' and
            residues[r2] == 'PHE' and
            residues[r3] == 'GLY'):
            return r1  # resid of DFG-Asp
    return None


def find_ape_in_pdb(pdb_path, dfg_resid):
    """Find APE motif (DFG+~20 residues onwards). Returns APE-Ala resid."""
    residues = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            ch = line[21]
            try:
                resid = int(line[22:26])
                resname = line[17:20].strip()
            except ValueError:
                continue
            if ch == 'D' or (ch == 'A' and resid > 103):
                residues[resid] = resname

    sorted_res = sorted(residues.keys())
    # Look for APE after DFG (typically +20-30 residues)
    for i, r1 in enumerate(sorted_res):
        if r1 <= dfg_resid + 15:
            continue
        if r1 > dfg_resid + 50:
            break
        r2 = r1 + 1
        r3 = r1 + 2
        if (r2 in residues and r3 in residues and
            residues[r1] == 'ALA' and
            residues[r2] == 'PRO' and
            residues[r3] == 'GLU'):
            return r1
    return None


def process_pose(poi, pose_id, dfg_resid, ape_resid):
    """
    Measure activation loop RMSD and DFG conformation for each frame.
    Returns DataFrame with per-frame metrics.
    """
    pdir   = RS_SIM_DIR / poi / pose_id
    tpr    = pdir / 'prod' / 'prod.tpr'
    xtc    = RS_ANA_DIR / poi / pose_id / 'prod_nojump.xtc'
    gro    = pdir / 'topol' / 'complex_protein.gro'

    for f in [tpr, xtc, gro]:
        if not Path(f).exists():
            return None

    atom_ranges = _detect_atom_ranges(str(gro))
    if atom_ranges is None:
        return None
    (va0, va1), (ka0, ka1) = atom_ranges

    try:
        u = mda.Universe(str(tpr), str(xtc))
    except Exception:
        return None

    protac  = u.select_atoms("resname UNL and not name H*")
    vhl_ca  = u.select_atoms(f"bynum {va0}:{va1} and name CA")
    kin_ca  = u.select_atoms(f"bynum {ka0}:{ka1} and name CA")

    if len(protac) == 0 or len(kin_ca) == 0:
        return None

    # Activation loop = DFG to APE
    if dfg_resid and ape_resid:
        aloop_sel = f"bynum {ka0}:{ka1} and name CA and resid {dfg_resid}:{ape_resid}"
        aloop = u.select_atoms(aloop_sel)
        if len(aloop) < 5:
            aloop = None
    else:
        aloop = None

    rows = []
    ref_pos = None

    for ts in u.trajectory[::STRIDE]:
        pp = protac.positions
        vp = vhl_ca.positions
        kp = kin_ca.positions

        dv = np.sqrt(((pp[:, None] - vp[None]) ** 2).sum(2)).min(0)
        dk = np.sqrt(((pp[:, None] - kp[None]) ** 2).sum(2)).min(0)
        is_bridging = int(dv.min() < VHL_CUT and dk.min() < KIN_CUT)

        aloop_rmsd = float('nan')
        if aloop is not None and len(aloop) > 0:
            ap = aloop.positions
            if ref_pos is None:
                ref_pos = ap.copy()
            aloop_rmsd = float(np.sqrt(((ap - ref_pos) ** 2).sum(1).mean()))

        rows.append({
            'time_ns':     float(ts.time / 1000),
            'is_bridging': is_bridging,
            'aloop_rmsd':  aloop_rmsd,
        })

    return pd.DataFrame(rows)


def main():
    rows = []
    for poi in POI_ORDER:
        label = LABELS[poi]
        pdbs  = sorted((RS_SIM_DIR / poi).glob('pose*/structure/complex_protein.pdb'))
        if not pdbs:
            continue

        # Find DFG/APE from first available PDB
        dfg_resid = ape_resid = None
        for pdb in pdbs:
            dr = find_dfg_in_pdb(pdb, 103)
            if dr:
                dfg_resid = dr
                ape_resid = find_ape_in_pdb(pdb, dr)
                break

        print(f"\n[{poi}] label={label}  DFG-Asp={dfg_resid}  APE-Ala={ape_resid}")

        for pdb in pdbs:
            pose_id = pdb.parts[-3]
            df_pose = process_pose(poi, pose_id, dfg_resid, ape_resid)
            if df_pose is None:
                print(f"  {pose_id}: FAIL")
                continue

            br  = df_pose[df_pose['is_bridging'] == 1]
            nbr = df_pose[df_pose['is_bridging'] == 0]
            frac_br = df_pose['is_bridging'].mean()

            aloop_br  = br['aloop_rmsd'].mean()  if len(br) > 0  else float('nan')
            aloop_nbr = nbr['aloop_rmsd'].mean() if len(nbr) > 0 else float('nan')
            aloop_delta = aloop_br - aloop_nbr if not (np.isnan(aloop_br) or np.isnan(aloop_nbr)) else float('nan')

            rows.append({
                'poi': poi, 'label': label, 'pose_id': pose_id,
                'dfg_resid': dfg_resid,
                'frac_bridging': round(frac_br, 4),
                'aloop_rmsd_bridging': round(aloop_br, 3) if not np.isnan(aloop_br) else float('nan'),
                'aloop_rmsd_nonbr': round(aloop_nbr, 3) if not np.isnan(aloop_nbr) else float('nan'),
                'aloop_delta_rmsd': round(aloop_delta, 3) if not np.isnan(aloop_delta) else float('nan'),
                'n_br_frames': len(br),
                'n_frames': len(df_pose),
            })
            print(f"  {pose_id}: frac_br={frac_br:.4f}  "
                  f"aloop_rmsd: br={aloop_br:.2f}  nbr={aloop_nbr:.2f}  "
                  f"delta={aloop_delta:.2f}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / 'activation_loop_all_poses.tsv', sep='\t', index=False)

    print(f"\n{'='*70}")
    print("SUMMARY — activation loop RMSD: bridging vs non-bridging frames")
    print('='*70)
    for poi in POI_ORDER:
        sub = df[df['poi']==poi].dropna(subset=['aloop_delta_rmsd'])
        if sub.empty: continue
        status = "DEG  " if LABELS[poi] else "NODEG"
        print(f"  [{status}] {poi:<16}  "
              f"aloop_delta={sub['aloop_delta_rmsd'].mean():.3f}±{sub['aloop_delta_rmsd'].std():.3f}  "
              f"n_poses_with_br={sub[sub['n_br_frames']>0].shape[0]}")

    # Degrader vs non-degrader
    deg  = df[df['label']==1].dropna(subset=['aloop_delta_rmsd'])
    ndeg = df[df['label']==0].dropna(subset=['aloop_delta_rmsd'])
    if len(deg) > 0 and len(ndeg) > 0:
        print(f"\n  aloop_delta_rmsd: degraders={deg['aloop_delta_rmsd'].mean():.3f}  "
              f"non-degraders={ndeg['aloop_delta_rmsd'].mean():.3f}")
        try:
            from scipy.stats import mannwhitneyu
            _, p = mannwhitneyu(deg['aloop_delta_rmsd'], ndeg['aloop_delta_rmsd'], alternative='two-sided')
            print(f"  Mann-Whitney p={p:.4f}")
        except ImportError:
            pass

    print(f"\nSaved: {OUT_DIR}/activation_loop_all_poses.tsv")


if __name__ == '__main__':
    main()
