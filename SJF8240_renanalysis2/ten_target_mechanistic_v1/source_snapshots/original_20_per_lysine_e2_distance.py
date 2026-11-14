#!/usr/bin/env python3
"""
Per-lysine E2 reachability across all 8 random-screen POIs.

Unlike 16_all_poi_e2_reachability.py (which takes min over ALL kinase lysines),
this script tracks each individual lysine NZ separately, giving:
  - Which specific lysine is the primary E2 target per POI
  - Per-lysine Δ80A: does bridging help THIS lysine reach E2?
  - Mechanistic insight for hard non-degraders (AXL, EPHA2, ABL1)

Key question: for AXL/EPHA2, is there even ONE lysine where bridging
significantly improves E2 access, or is the geometry globally bad?

Output:
  analysis/per_lysine_e2_all_pois.tsv   — one row per (poi, pose_id, resid)
  analysis/per_lysine_e2_summary.tsv    — one row per (poi, resid), avg across poses
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import MDAnalysis as mda
from MDAnalysis import transformations
import warnings
warnings.filterwarnings('ignore')

ROOT       = Path.home() / 'protac-md'
REFDIR     = ROOT / 'reference_structures'
SIMDIR     = ROOT / 'simulations'
E2_ENS     = REFDIR / 'e2_c85_ensemble.json'
OUT_DIR    = ROOT / 'analysis'

BRIDGE_CUTOFF = 5.0
E2_THRESH_80  = 80.0
E2_THRESH_50  = 50.0
N_WORKERS     = 20

POI_CONFIG = {
    'ABL1':        {'label': 0},
    'AXL':         {'label': 0},
    'DDR2':        {'label': 1},
    'EPHA2':       {'label': 0},
    'EPHB2':       {'label': 1, 'use_pose_pdb': True},
    'MAPK14_full': {'label': 1},
    'MET':         {'label': 1},
    'RIPK2_full':  {'label': 1},
}


def load_e2_ensemble():
    with open(E2_ENS) as f:
        data = json.load(f)
    if isinstance(data, list):
        return np.array(data)
    if isinstance(data, dict):
        all_pos = []
        for v in data.values():
            all_pos.extend(v)
        return np.array(all_pos)
    raise ValueError(f"Unexpected format in {E2_ENS}")


def kabsch_rotation(mobile_c, fixed_c):
    H = mobile_c.T @ fixed_c
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1.0, 1.0, d]) @ U.T


def get_elongc_pdb(poi, cfg):
    main_pdb = SIMDIR / poi / 'structure' / 'rank1' / 'complex_protein.pdb'
    if main_pdb.exists() and not cfg.get('use_pose_pdb'):
        u   = mda.Universe(str(main_pdb))
        elc = u.select_atoms('chainID B and name CA')
        if len(elc) == 85:
            return elc.positions.copy()
    for pdb in sorted((ROOT / 'random_screen' / poi).glob('pose*/structure/complex_protein.pdb')):
        u   = mda.Universe(str(pdb))
        elc = u.select_atoms('chainID B and name CA')
        if len(elc) == 85:
            return elc.positions.copy()
    raise FileNotFoundError(f"No valid ElonginC PDB for {poi}")


def _worker(args):
    """Per-pose worker: returns per-lysine E2 distance stats."""
    poi, pose_id, label, e2_pdb_ens, elongc_pdb_coords = args

    rs_sim = ROOT / 'random_screen' / poi
    rs_ana = ROOT / 'analysis' / 'random_screen' / poi

    tpr = rs_sim / pose_id / 'prod' / 'prod.tpr'
    xtc = rs_ana / pose_id / 'prod_nojump.xtc'
    csv = rs_ana / pose_id / 'metrics.csv'

    if not all(p.exists() for p in [tpr, xtc, csv]):
        return None

    # Load bridging mask from metrics.csv
    rows = []
    for line in open(csv):
        if line.startswith(('#', 'time')):
            continue
        parts = line.strip().split(',')
        if len(parts) >= 3:
            try:
                rows.append((float(parts[0]), float(parts[1]), float(parts[2])))
            except ValueError:
                continue
    if not rows:
        return None

    arr      = np.array(rows)
    bridging = (arr[:, 1] < BRIDGE_CUTOFF) & (arr[:, 2] < BRIDGE_CUTOFF)
    bridge_pct = bridging.mean() * 100

    try:
        warnings.filterwarnings('ignore')
        u = mda.Universe(str(tpr), str(xtc))
        u.trajectory.add_transformations(
            transformations.unwrap(u.select_atoms('protein')))
    except Exception:
        return None

    elongc_traj = u.select_atoms('segid *Protein_chain_B and name CA')
    lys_ag      = u.select_atoms('segid *Protein_chain_D and resname LYS and name NZ')

    if len(elongc_traj) != 85 or len(lys_ag) == 0:
        return None

    # Get per-lysine residue IDs and unique residues
    lys_resids = lys_ag.resids
    unique_resids = np.unique(lys_resids)  # some may have multiple NZ if model artifact
    n_lys    = len(unique_resids)
    pdb_c    = elongc_pdb_coords.mean(0)
    n_frames = min(u.trajectory.n_frames, len(bridging))

    # Per-frame, per-lysine E2 distance arrays
    per_lys_e2   = np.zeros((n_frames, n_lys))  # min E2 dist per lysine per frame

    for fi, _ts in enumerate(u.trajectory[:n_frames]):
        mob   = elongc_traj.positions
        mob_c = mob.mean(0)
        R     = kabsch_rotation(elongc_pdb_coords - pdb_c, mob - mob_c)
        e2_sim = (e2_pdb_ens - pdb_c) @ R.T + mob_c  # (n_e2, 3)

        nz_pos = lys_ag.positions  # (n_nz_atoms, 3)

        for li, resid in enumerate(unique_resids):
            # Get NZ atoms for this residue
            mask = lys_resids == resid
            nz_r = nz_pos[mask]  # (1 or more, 3)
            # Min distance from this lysine's NZ to any E2 position
            d = np.sqrt(((nz_r[:, None] - e2_sim[None]) ** 2).sum(2))  # (nz, n_e2)
            per_lys_e2[fi, li] = d.min()

    b  = bridging[:n_frames].astype(bool)
    nb = ~b

    def sm(arr2d, mask):
        if mask.sum() == 0:
            return np.full(arr2d.shape[1], float('nan'))
        return arr2d[mask].mean(0)

    mean_bridge    = sm(per_lys_e2, b)
    mean_nonbridge = sm(per_lys_e2, nb)
    delta          = mean_bridge - mean_nonbridge  # negative = bridging brings E2 closer

    reach_80_b  = sm((per_lys_e2 < E2_THRESH_80).astype(float), b)
    reach_80_nb = sm((per_lys_e2 < E2_THRESH_80).astype(float), nb)
    reach_50_b  = sm((per_lys_e2 < E2_THRESH_50).astype(float), b)
    reach_50_nb = sm((per_lys_e2 < E2_THRESH_50).astype(float), nb)

    result_rows = []
    for li, resid in enumerate(unique_resids):
        result_rows.append({
            'poi':                  poi,
            'pose_id':              pose_id,
            'label':                label,
            'resid':                int(resid),
            'bridge_pct':           bridge_pct,
            'mean_e2_bridge':       float(mean_bridge[li]),
            'mean_e2_nonbridge':    float(mean_nonbridge[li]),
            'delta_e2':             float(delta[li]),        # negative = bridging helps
            'reach_80A_bridge':     float(reach_80_b[li]),
            'reach_80A_nonbridge':  float(reach_80_nb[li]),
            'delta_reach_80A':      float(reach_80_b[li] - reach_80_nb[li]),
            'reach_50A_bridge':     float(reach_50_b[li]),
            'reach_50A_nonbridge':  float(reach_50_nb[li]),
            'delta_reach_50A':      float(reach_50_b[li] - reach_50_nb[li]),
        })

    return result_rows


def main():
    print("=== Per-Lysine E2 Reachability — All 8 POIs (parallel) ===\n")

    e2_pdb_ens = load_e2_ensemble()
    print(f"E2 ensemble: {len(e2_pdb_ens)} positions\n")

    tasks = []
    for poi, cfg in POI_CONFIG.items():
        try:
            elc = get_elongc_pdb(poi, cfg)
        except FileNotFoundError as e:
            print(f"SKIP {poi}: {e}")
            continue
        rs_ana = ROOT / 'analysis' / 'random_screen' / poi
        for pd_dir in sorted(rs_ana.glob('pose*')):
            tasks.append((poi, pd_dir.name, cfg['label'], e2_pdb_ens, elc))

    print(f"Submitting {len(tasks)} pose-jobs across {N_WORKERS} workers...\n")

    all_rows = []
    completed = 0
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futs = {ex.submit(_worker, t): t for t in tasks}
        for fut in as_completed(futs):
            poi, pose_id, *_ = futs[fut]
            completed += 1
            try:
                r = fut.result()
                if r:
                    all_rows.extend(r)
                    print(f"  [{completed:3d}/{len(tasks)}] {poi}/{pose_id}  "
                          f"{len(r)} lysines", flush=True)
            except Exception as e:
                print(f"  [{completed:3d}/{len(tasks)}] {poi}/{pose_id} ERROR: {e}",
                      flush=True)

    if not all_rows:
        print("No data.")
        return

    df = pd.DataFrame(all_rows)
    pose_out = OUT_DIR / 'per_lysine_e2_all_pois.tsv'
    df.to_csv(pose_out, sep='\t', index=False, float_format='%.4f')
    print(f"\nPer-pose-per-lysine saved → {pose_out}")

    # Per-POI per-lysine summary (average across poses)
    summary = []
    for (poi, resid), grp in df.groupby(['poi', 'resid']):
        label = int(grp['label'].iloc[0])
        summary.append({
            'poi':                  poi,
            'label':                label,
            'resid':                int(resid),
            'n_poses':              len(grp),
            'mean_e2_bridge':       grp['mean_e2_bridge'].mean(),
            'mean_e2_nonbridge':    grp['mean_e2_nonbridge'].mean(),
            'delta_e2':             grp['delta_e2'].mean(),
            'reach_80A_bridge':     grp['reach_80A_bridge'].mean(),
            'reach_80A_nonbridge':  grp['reach_80A_nonbridge'].mean(),
            'delta_reach_80A':      grp['delta_reach_80A'].mean(),
            'reach_50A_bridge':     grp['reach_50A_bridge'].mean(),
            'reach_50A_nonbridge':  grp['reach_50A_nonbridge'].mean(),
            'delta_reach_50A':      grp['delta_reach_50A'].mean(),
        })

    sum_df = pd.DataFrame(summary)
    sum_out = OUT_DIR / 'per_lysine_e2_summary.tsv'
    sum_df.to_csv(sum_out, sep='\t', index=False, float_format='%.4f')
    print(f"Per-POI-per-lysine summary → {sum_out}\n")

    # Print best lysine per POI (lowest mean E2 dist in bridging frames)
    print(f"{'='*72}")
    print(f"Best lysine per POI (lowest mean E2 dist when bridging):\n")
    print(f"{'POI':<15} {'Label':<5} {'Resid':>6} {'E2_b(Å)':>9} "
          f"{'E2_nb(Å)':>9} {'ΔE2':>8} {'Δ80A':>8}")
    print('─' * 72)

    for poi in POI_CONFIG:
        sub = sum_df[sum_df['poi'] == poi]
        if sub.empty:
            continue
        best = sub.loc[sub['mean_e2_bridge'].idxmin()]
        lbl  = 'DEG' if best['label'] else 'non'
        print(f"{poi:<15} {lbl:<5} K{int(best['resid']):>5} "
              f"{best['mean_e2_bridge']:>9.1f} {best['mean_e2_nonbridge']:>9.1f} "
              f"{best['delta_e2']:>+8.2f} {best['delta_reach_80A']:>+8.4f}")

    # AXL and EPHA2 deep dive — all lysines sorted by E2 dist
    for focus_poi in ['AXL', 'EPHA2', 'ABL1']:
        sub = sum_df[sum_df['poi'] == focus_poi].sort_values('mean_e2_bridge')
        if sub.empty:
            continue
        print(f"\n  {focus_poi} — all lysines sorted by E2 dist (bridging):")
        print(f"  {'Resid':>6} {'E2_b(Å)':>9} {'E2_nb(Å)':>9} "
              f"{'ΔE2':>8} {'R80_b':>7} {'R80_nb':>7} {'Δ80A':>7}")
        for _, r in sub.iterrows():
            print(f"  K{int(r['resid']):>5} {r['mean_e2_bridge']:>9.1f} "
                  f"{r['mean_e2_nonbridge']:>9.1f} {r['delta_e2']:>+8.2f} "
                  f"{r['reach_80A_bridge']:>7.3f} {r['reach_80A_nonbridge']:>7.3f} "
                  f"{r['delta_reach_80A']:>+7.4f}")

    print(f"\nDone. Outputs: {pose_out.name}, {sum_out.name}")


if __name__ == '__main__':
    main()
