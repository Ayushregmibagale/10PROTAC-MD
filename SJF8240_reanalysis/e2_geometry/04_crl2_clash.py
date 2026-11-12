#!/usr/bin/env python3
"""
CRL2 clash analysis: does the kinase in the ternary complex sterically occlude
the CUL2 scaffold or RBX1 (E2 docking platform)?

Method:
  1. Load 5N4W (CUL2-A, VHL-V, RBX1-R, ELOB-B, ELOC-C)
  2. For each of the 20 pRosettaC poses per POI, superimpose the pose onto 5N4W
     by Kabsch-aligning ELOB+ELOC (chains B+C = same chains in both structures)
  3. After alignment, place the kinase (chain D) in the CRL2 frame
  4. Measure:
     - Min distance from kinase CAs to CUL2 (chain A)
     - Min distance from kinase CAs to RBX1 (chain R) — the E2-docking ring
     - Number of kinase CAs within 10/15/20 Å of RBX1
     - Number of kinase CAs within 10/15/20 Å of CUL2

Run: ~/.conda/envs/gromacs/bin/python3 -u 04_crl2_clash.py
"""

# RETIRED_BY_FOLDER_AUDIT_20261112: original source is preserved in the audit archive.
if __name__ == "__main__":
    raise SystemExit("Retired exploratory entry point: see folder_audit_20261112/AUDIT.md and use RUN_ANALYSIS.py. Prior outputs are preserved but not certified by this script.")

import numpy as np
import pandas as pd
from pathlib import Path

SCRIPT_DIR  = Path(__file__).resolve().parent
RSCREEN_DIR = Path.home() / 'protac-md' / 'random_screen'
REFDIR      = Path.home() / 'protac-md' / 'reference_structures'
CRL2_PDB    = REFDIR / '5N4W.pdb'
OUT_DIR     = SCRIPT_DIR / 'results'
OUT_DIR.mkdir(parents=True, exist_ok=True)

LABELS = {
    'DDR2':1,'ABL1':0,'AXL':0,'EPHA2':1,
    'EPHB2':1,'MAPK14_full':1,'MET':1,'RIPK2_full':1,
}
POI_ORDER = ['DDR2','ABL1','AXL','EPHA2','EPHB2','MAPK14_full','MET','RIPK2_full']
CUTOFFS = [10, 15, 20, 30]


def parse_ca(pdb_path, chains):
    """Parse Cα atoms for given chains. Returns dict chain → (n, 3) array."""
    out = {c: {} for c in chains}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            if line[12:16].strip() != 'CA':
                continue
            ch = line[21]
            if ch not in chains:
                continue
            try:
                resid = int(line[22:26])
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            out[ch][resid] = (x, y, z)
    return {c: np.array(list(d.values()), dtype=np.float32)
            for c, d in out.items() if d}


def parse_ca_pose(pdb_path):
    """
    Parse pRosettaC pose PDB.
    Chain A resid 1-103 → VHL. Chain B → ELOB. Chain C → ELOC. Chain D → kinase.
    Returns vhl, elob, eloc, kinase arrays.
    """
    out = {'vhl': {}, 'elob': {}, 'eloc': {}, 'kin': {}}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            if line[12:16].strip() != 'CA':
                continue
            ch = line[21]
            try:
                resid = int(line[22:26])
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            pos = (x, y, z)
            if ch == 'A' and 1 <= resid <= 103:
                out['vhl'][resid] = pos
            elif ch == 'B':
                out['elob'][resid] = pos
            elif ch == 'C':
                out['eloc'][resid] = pos
            elif ch == 'D':
                out['kin'][resid] = pos
            elif ch == 'A' and resid > 103:
                out['kin'][resid] = pos
    return {k: np.array(list(v.values()), dtype=np.float32)
            for k, v in out.items() if v}


def kabsch(mobile, fixed):
    """
    Kabsch rotation: aligns mobile onto fixed (both already mean-centred).
    Returns R such that mobile @ R.T ≈ fixed.
    """
    H = mobile.T @ fixed
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1., 1., d]) @ U.T


def align_and_transform(pose_chains, ref_chains):
    """
    Align pose ELOB+ELOC onto reference ELOB+ELOC (5N4W chains B+C).
    Applies the same transform to kinase positions.
    Returns aligned kinase positions in the 5N4W frame.
    """
    # Build matched reference array: concatenate ELOB+ELOC
    ref_b = ref_chains.get('B', np.empty((0, 3)))
    ref_c = ref_chains.get('C', np.empty((0, 3)))
    pos_b = pose_chains.get('elob', np.empty((0, 3)))
    pos_c = pose_chains.get('eloc', np.empty((0, 3)))

    # Trim to same length (take min)
    nb = min(len(ref_b), len(pos_b))
    nc = min(len(ref_c), len(pos_c))

    if nb + nc < 10:
        return None  # too few atoms for reliable alignment

    ref_m  = np.vstack([ref_b[:nb], ref_c[:nc]])
    pose_m = np.vstack([pos_b[:nb], pos_c[:nc]])

    ref_cen  = ref_m.mean(0)
    pose_cen = pose_m.mean(0)

    R = kabsch(pose_m - pose_cen, ref_m - ref_cen)

    # Transform kinase
    kin = pose_chains.get('kin')
    if kin is None or len(kin) == 0:
        return None
    kin_aligned = (kin - pose_cen) @ R.T + ref_cen
    return kin_aligned


def min_dist(pts_a, pts_b, chunk=500):
    md = np.inf
    for i in range(0, len(pts_a), chunk):
        d = np.sqrt(((pts_a[i:i+chunk, None] - pts_b[None]) ** 2).sum(2))
        md = min(md, float(d.min()))
    return md


def frac_within(pts, cloud, cutoff, chunk=200):
    cnt = 0
    for i in range(0, len(pts), chunk):
        d = np.sqrt(((pts[i:i+chunk, None] - cloud[None]) ** 2).sum(2))
        cnt += int((d.min(1) < cutoff).sum())
    return cnt / max(len(pts), 1)


def main():
    # Load 5N4W reference
    ref = parse_ca(CRL2_PDB, chains={'A', 'B', 'C', 'R', 'V'})
    cul2 = ref.get('A', np.empty((0, 3)))
    rbx1 = ref.get('R', np.empty((0, 3)))
    elob_ref = ref.get('B', np.empty((0, 3)))
    eloc_ref = ref.get('C', np.empty((0, 3)))
    ref_chains = {'B': elob_ref, 'C': eloc_ref}

    print(f"5N4W loaded: CUL2={len(cul2)} CA, RBX1={len(rbx1)} CA, "
          f"ELOB={len(elob_ref)} CA, ELOC={len(eloc_ref)} CA")

    rows = []
    for poi in POI_ORDER:
        label = LABELS[poi]
        pdbs  = sorted((RSCREEN_DIR / poi).glob('pose*/structure/complex_protein.pdb'))
        print(f"\n[{poi}] label={label}  n_poses={len(pdbs)}")

        for pdb in pdbs:
            pose_id = pdb.parts[-3]
            pose_chains = parse_ca_pose(pdb)

            kin_aligned = align_and_transform(pose_chains, ref_chains)
            if kin_aligned is None:
                print(f"  {pose_id}: ALIGN FAIL")
                continue

            md_cul2 = min_dist(kin_aligned, cul2)
            md_rbx1 = min_dist(kin_aligned, rbx1)

            frac_cul2 = {c: frac_within(kin_aligned, cul2, c) for c in CUTOFFS}
            frac_rbx1 = {c: frac_within(kin_aligned, rbx1, c) for c in CUTOFFS}

            rows.append({
                'poi': poi, 'label': label, 'pose_id': pose_id,
                'min_d_cul2': round(md_cul2, 1),
                'min_d_rbx1': round(md_rbx1, 1),
                **{f'frac_cul2_{c}A': round(frac_cul2[c], 4) for c in CUTOFFS},
                **{f'frac_rbx1_{c}A': round(frac_rbx1[c], 4) for c in CUTOFFS},
            })
            print(f"  {pose_id}: min_d_CUL2={md_cul2:.1f} Å  min_d_RBX1={md_rbx1:.1f} Å  "
                  f"kin@15A_RBX1={frac_rbx1[15]:.3f}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / 'crl2_clash_all_poses.tsv', sep='\t', index=False)

    print(f"\n{'='*70}")
    print("SUMMARY (mean ± std across 20 poses)")
    print('='*70)
    summary = []
    for poi in POI_ORDER:
        sub = df[df['poi']==poi]
        if sub.empty: continue
        r = {'poi': poi, 'label': LABELS[poi]}
        for col in ['min_d_cul2','min_d_rbx1','frac_cul2_15A','frac_rbx1_15A','frac_rbx1_20A']:
            r[f'{col}_mean'] = sub[col].mean()
            r[f'{col}_std']  = sub[col].std()
        summary.append(r)
        status = "DEG  " if LABELS[poi] else "NODEG"
        print(f"  [{status}] {poi:<16}  "
              f"CUL2={r['min_d_cul2_mean']:.1f}±{r['min_d_cul2_std']:.1f} Å  "
              f"RBX1={r['min_d_rbx1_mean']:.1f}±{r['min_d_rbx1_std']:.1f} Å  "
              f"kin@15A_RBX1={r['frac_rbx1_15A_mean']:.3f}")

    df_sum = pd.DataFrame(summary)
    df_sum.to_csv(OUT_DIR / 'crl2_clash_summary.tsv', sep='\t', index=False, float_format='%.3f')

    # Degrader vs non-degrader
    deg  = df[df['label']==1]
    ndeg = df[df['label']==0]
    print(f"\n{'':20s}  {'degraders':>14}  {'non-degraders':>14}")
    for col in ['min_d_cul2','min_d_rbx1','frac_rbx1_15A','frac_rbx1_20A']:
        print(f"  {col:<22}  {deg[col].mean():>8.2f}±{deg[col].std():.2f}  "
              f"{ndeg[col].mean():>8.2f}±{ndeg[col].std():.2f}")

    try:
        from scipy.stats import mannwhitneyu
        print("\n  Mann-Whitney p-values:")
        for col in ['min_d_cul2','min_d_rbx1','frac_rbx1_15A']:
            _, p = mannwhitneyu(deg[col], ndeg[col], alternative='two-sided')
            print(f"    {col:<25} p={p:.4f}")
    except ImportError:
        pass

    print(f"\nSaved: {OUT_DIR}/crl2_clash_all_poses.tsv")


if __name__ == '__main__':
    main()
