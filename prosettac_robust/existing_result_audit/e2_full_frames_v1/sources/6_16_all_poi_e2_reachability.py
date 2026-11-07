#!/usr/bin/env python3
"""
E2 reachability across all 8 random-screen POIs — parallel version.

All ~133 poses processed in parallel using ProcessPoolExecutor (20 workers).
Expected runtime: ~5-10 min on a 24-core machine vs ~67 min single-threaded.

Key question: does ternary-complex bridging INCREASE E2 C85 reachability
to kinase lysines? Compare Δ(reach_80A) = bridging − non-bridging across
degraders vs non-degraders.

Degraders (label=1):     DDR2, EPHB2, MAPK14_full, MET, RIPK2_full
Non-degraders (label=0): ABL1, AXL, EPHA2
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
E2_THRESHOLDS = [50, 80, 100]
N_WORKERS     = 20   # leave 4 cores for system

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


# ── worker function (runs in subprocess) ─────────────────────────────────────

def _worker(args):
    """Top-level picklable worker — processes one pose."""
    poi, pose_id, label, e2_pdb_ens, elongc_pdb_coords = args

    rs_sim = ROOT / 'random_screen' / poi
    rs_ana = ROOT / 'analysis' / 'random_screen' / poi

    tpr = rs_sim / pose_id / 'prod' / 'prod.tpr'
    xtc = rs_ana / pose_id / 'prod_nojump.xtc'
    csv = rs_ana / pose_id / 'metrics.csv'

    if not all(p.exists() for p in [tpr, xtc, csv]):
        return None

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

    arr        = np.array(rows)
    bridging   = (arr[:, 1] < BRIDGE_CUTOFF) & (arr[:, 2] < BRIDGE_CUTOFF)
    bridge_pct = bridging.mean() * 100

    try:
        import warnings
        warnings.filterwarnings('ignore')
        u = mda.Universe(str(tpr), str(xtc))
        u.trajectory.add_transformations(
            transformations.unwrap(u.select_atoms('protein')))
    except Exception as e:
        return None

    elongc_traj = u.select_atoms('segid *Protein_chain_B and name CA')
    lys_nz      = u.select_atoms('segid *Protein_chain_D and resname LYS and name NZ')

    if len(elongc_traj) != 85 or len(lys_nz) == 0:
        return None

    pdb_c    = elongc_pdb_coords.mean(0)
    n_frames = min(u.trajectory.n_frames, len(bridging))
    min_e2   = np.zeros(n_frames)
    reach    = {thr: np.zeros(n_frames) for thr in E2_THRESHOLDS}

    for fi, _ts in enumerate(u.trajectory[:n_frames]):
        mob    = elongc_traj.positions
        mob_c  = mob.mean(0)
        R      = kabsch_rotation(elongc_pdb_coords - pdb_c, mob - mob_c)
        e2_sim = (e2_pdb_ens - pdb_c) @ R.T + mob_c
        nz     = lys_nz.positions
        dists  = np.sqrt(((nz[:, None] - e2_sim[None]) ** 2).sum(2))
        min_e2[fi] = dists.min()
        for thr in E2_THRESHOLDS:
            reach[thr][fi] = (dists.min(0) < thr).mean()

    b  = bridging[:n_frames].astype(bool)
    nb = ~b

    def sm(a, m):
        return float(a[m].mean()) if m.sum() > 0 else float('nan')

    return {
        'poi': poi, 'pose_id': pose_id, 'label': label,
        'bridge_pct':      bridge_pct,
        'n_kinase_lys':    len(lys_nz),
        'min_e2_dist_mean': float(min_e2.mean()),
        **{f'reach_{thr}A_bridge':    sm(reach[thr], b)  for thr in E2_THRESHOLDS},
        **{f'reach_{thr}A_nonbridge': sm(reach[thr], nb) for thr in E2_THRESHOLDS},
    }


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    print("=== E2 Reachability — All 8 POIs (parallel) ===\n")

    e2_pdb_ens = load_e2_ensemble()
    print(f"E2 ensemble: {len(e2_pdb_ens)} positions, centroid {e2_pdb_ens.mean(0).round(1)}\n")

    # Build task list
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
                    all_rows.append(r)
                    print(f"  [{completed:3d}/{len(tasks)}] {poi}/{pose_id}  "
                          f"bridge={r['bridge_pct']:.0f}%  "
                          f"min_e2={r['min_e2_dist_mean']:.0f} Å  "
                          f"Δ80={r['reach_80A_bridge']-r['reach_80A_nonbridge']:+.3f}",
                          flush=True)
            except Exception as e:
                print(f"  [{completed:3d}/{len(tasks)}] {poi}/{pose_id} ERROR: {e}",
                      flush=True)

    if not all_rows:
        print("No data collected.")
        return

    df = pd.DataFrame(all_rows)
    pose_out = OUT_DIR / 'e2_reachability_all_pois.tsv'
    df.to_csv(pose_out, sep='\t', index=False, float_format='%.4f')

    # Per-POI summary
    summary = []
    for poi, grp in df.groupby('poi'):
        b80 = grp['reach_80A_bridge'].mean(skipna=True)
        n80 = grp['reach_80A_nonbridge'].mean(skipna=True)
        b100 = grp['reach_100A_bridge'].mean(skipna=True)
        n100 = grp['reach_100A_nonbridge'].mean(skipna=True)
        summary.append({
            'poi':                  poi,
            'label':                int(grp['label'].iloc[0]),
            'n_poses':              len(grp),
            'mean_bridge_pct':      grp['bridge_pct'].mean(),
            'mean_min_e2_dist_A':   grp['min_e2_dist_mean'].mean(),
            'reach_80A_bridge':     b80,
            'reach_80A_nonbridge':  n80,
            'delta_reach_80A':      b80 - n80,
            'reach_100A_bridge':    b100,
            'reach_100A_nonbridge': n100,
            'delta_reach_100A':     b100 - n100,
        })

    sum_df = pd.DataFrame(summary).sort_values(['label', 'delta_reach_80A'],
                                               ascending=[False, False])
    sum_out = OUT_DIR / 'e2_reachability_poi_summary.tsv'
    sum_df.to_csv(sum_out, sep='\t', index=False, float_format='%.4f')

    # Final table
    print(f"\n{'='*72}")
    print(f"{'POI':<15} {'Label':<5} {'n':>3} {'Bridge%':>8} {'MinE2Å':>8} "
          f"{'R80_b':>7} {'R80_nb':>7} {'Δ80':>7}")
    print('─' * 72)
    for _, r in sum_df.iterrows():
        lbl = 'DEG' if r['label'] else 'non'
        print(f"{r['poi']:<15} {lbl:<5} {r['n_poses']:>3.0f} "
              f"{r['mean_bridge_pct']:>8.1f} {r['mean_min_e2_dist_A']:>8.1f} "
              f"{r['reach_80A_bridge']:>7.3f} {r['reach_80A_nonbridge']:>7.3f} "
              f"{r['delta_reach_80A']:>+7.3f}")

    deg    = sum_df[sum_df['label'] == 1]
    nondeg = sum_df[sum_df['label'] == 0]
    print('─' * 72)
    print(f"\n  Degraders     (n={len(deg)}): mean Δ80A = {deg['delta_reach_80A'].mean():+.4f}, "
          f"mean min_E2 = {deg['mean_min_e2_dist_A'].mean():.1f} Å")
    print(f"  Non-degraders (n={len(nondeg)}): mean Δ80A = {nondeg['delta_reach_80A'].mean():+.4f}, "
          f"mean min_E2 = {nondeg['mean_min_e2_dist_A'].mean():.1f} Å")
    print(f"\nSaved → {pose_out.name}, {sum_out.name}")


if __name__ == '__main__':
    main()
