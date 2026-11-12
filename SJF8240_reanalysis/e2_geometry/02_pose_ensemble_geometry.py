#!/usr/bin/env python3
"""
Compute H1/H2/H3 geometry metrics across all 20 random_screen poses per POI.

Provides distributions rather than single rank-1 values.
Tests whether ABL1/AXL systematically show larger VHL→kin / VHL→E2 angles
across all 20 starting poses.

Run: ~/.conda/envs/gromacs/bin/python3 -u 02_pose_ensemble_geometry.py
"""

# RETIRED_BY_FOLDER_AUDIT_20261112: original source is preserved in the audit archive.
if __name__ == "__main__":
    raise SystemExit("Retired exploratory entry point: see folder_audit_20261112/AUDIT.md and use RUN_ANALYSIS.py. Prior outputs are preserved but not certified by this script.")

import json
import numpy as np
import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

SCRIPT_DIR  = Path(__file__).resolve().parent
RSCREEN_DIR = Path.home() / 'protac-md' / 'random_screen'
REFDIR      = Path.home() / 'protac-md' / 'reference_structures'
E2_ENS      = REFDIR / 'e2_c85_ensemble.json'
OUT_DIR     = SCRIPT_DIR / 'results'
OUT_DIR.mkdir(parents=True, exist_ok=True)

LABELS = {
    'DDR2':1,'ABL1':0,'AXL':0,'EPHA2':1,
    'EPHB2':1,'MAPK14_full':1,'MET':1,'RIPK2_full':1,
}
POI_ORDER = ['DDR2','ABL1','AXL','EPHA2','EPHB2','MAPK14_full','MET','RIPK2_full']


def load_e2_ensemble():
    with open(E2_ENS) as f:
        data = json.load(f)
    if isinstance(data, list):
        return np.array(data, dtype=np.float32)
    return np.array([p for v in data.values() for p in v], dtype=np.float32)


def parse_ca_from_pdb(pdb_path):
    """Parse VHL (chain A, resid 1-103) and kinase (chain A resid>103 or chain D) CAs."""
    vhl, kin = {}, {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            if line[12:16].strip() != 'CA':
                continue
            chain = line[21]
            try:
                resid = int(line[22:26])
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            if chain == 'A':
                if 1 <= resid <= 103:
                    vhl[resid] = (x, y, z)
                elif resid > 103:
                    kin[resid] = (x, y, z)
            elif chain == 'D':
                kin[resid] = (x, y, z)
    if not vhl or not kin:
        return None, None
    vhl_pos = np.array([vhl[r] for r in sorted(vhl)], dtype=np.float32)
    kin_pos = np.array([kin[r] for r in sorted(kin)], dtype=np.float32)
    return vhl_pos, kin_pos


def angle_deg(v1, v2):
    c = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9)
    return float(np.degrees(np.arccos(np.clip(c, -1, 1))))


def min_dist_to_cloud(points, cloud, chunk=500):
    min_d = np.inf
    for i in range(0, len(points), chunk):
        d = np.sqrt(((points[i:i+chunk, None] - cloud[None]) ** 2).sum(2))
        min_d = min(min_d, float(d.min()))
    return min_d


def fraction_within(points, cloud, cutoff, chunk=200):
    count = 0
    for i in range(0, len(points), chunk):
        d = np.sqrt(((points[i:i+chunk, None] - cloud[None]) ** 2).sum(2))
        count += int((d.min(1) < cutoff).sum())
    return count / max(len(points), 1)


def main():
    e2  = load_e2_ensemble()
    e2c = e2.mean(0)
    print(f"E2 ensemble: {len(e2)} pts  centroid={e2c.round(1)}")

    rows = []

    for poi in POI_ORDER:
        label = LABELS[poi]
        pose_dirs = sorted((RSCREEN_DIR / poi).glob('pose*/structure/complex_protein.pdb'))
        print(f"\n[{poi}] label={label}  n_poses={len(pose_dirs)}")

        for pdb in pose_dirs:
            pose_id = pdb.parts[-3]  # pose001, pose002, ...
            vhl_pos, kin_pos = parse_ca_from_pdb(pdb)
            if vhl_pos is None:
                print(f"  {pose_id}: PARSE FAIL")
                continue

            vhl_cen = vhl_pos.mean(0)
            kin_cen = kin_pos.mean(0)

            # H1
            v_e2  = e2c - vhl_cen
            v_kin = kin_cen - vhl_cen
            ang   = angle_deg(v_kin, v_e2)
            dist_cen = float(np.linalg.norm(kin_cen - e2c))

            # H2
            min_d  = min_dist_to_cloud(kin_pos, e2)
            f50    = fraction_within(kin_pos, e2, 50)
            f80    = fraction_within(kin_pos, e2, 80)
            rg     = float(np.sqrt(((kin_pos - kin_cen)**2).sum(1).mean()))

            # H3: N-lobe vs C-lobe (40/60 split)
            split  = int(len(kin_pos) * 0.4)
            nl_cen_d = float(np.linalg.norm(kin_pos[:split].mean(0) - e2c))
            cl_cen_d = float(np.linalg.norm(kin_pos[split:].mean(0) - e2c))

            rows.append({
                'poi': poi, 'label': label, 'pose_id': pose_id,
                'ang_kin_e2': round(ang, 1),
                'dist_kin_e2_cen': round(dist_cen, 1),
                'min_d_kin_e2': round(min_d, 1),
                'frac_e2_50A': round(f50, 3),
                'frac_e2_80A': round(f80, 3),
                'rg_kin': round(rg, 1),
                'nlobe_cen_d': round(nl_cen_d, 1),
                'clobe_cen_d': round(cl_cen_d, 1),
            })
            print(f"  {pose_id}  ang={ang:.1f}°  min_d={min_d:.1f}  f50={f50:.3f}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / 'e2_geometry_all_poses.tsv', sep='\t', index=False)

    # ── Summary stats per POI ─────────────────────────────────────────────────
    print(f"\n{'='*70}")
    print("PER-POI SUMMARY (mean ± std across 20 poses)")
    print('='*70)
    summary_rows = []
    for poi in POI_ORDER:
        sub = df[df['poi']==poi]
        if sub.empty:
            continue
        r = {
            'poi': poi, 'label': LABELS[poi],
            'ang_mean': sub['ang_kin_e2'].mean(),
            'ang_std':  sub['ang_kin_e2'].std(),
            'ang_min':  sub['ang_kin_e2'].min(),
            'ang_max':  sub['ang_kin_e2'].max(),
            'dist_mean': sub['dist_kin_e2_cen'].mean(),
            'min_d_mean': sub['min_d_kin_e2'].mean(),
            'frac50_mean': sub['frac_e2_50A'].mean(),
            'frac80_mean': sub['frac_e2_80A'].mean(),
        }
        summary_rows.append(r)
        status = "DEG  " if LABELS[poi] else "NODEG"
        print(f"  [{status}] {poi:<16}  ang={r['ang_mean']:.1f}±{r['ang_std']:.1f}°  "
              f"({r['ang_min']:.1f}–{r['ang_max']:.1f})  "
              f"dist={r['dist_mean']:.1f}  f50={r['frac50_mean']:.3f}  f80={r['frac80_mean']:.3f}")

    df_sum = pd.DataFrame(summary_rows)
    df_sum.to_csv(OUT_DIR / 'e2_geometry_poi_summary.tsv', sep='\t', index=False, float_format='%.2f')

    # Degrader vs non-degrader comparison
    print(f"\n{'':20s}  {'degraders (n=6)':>18}  {'non-degraders (n=2)':>20}")
    deg  = df[df['label']==1]
    ndeg = df[df['label']==0]
    for col, name in [('ang_kin_e2','angle(°)'),
                      ('dist_kin_e2_cen','dist_cen(Å)'),
                      ('min_d_kin_e2','min_d(Å)'),
                      ('frac_e2_50A','frac@50A'),
                      ('frac_e2_80A','frac@80A')]:
        dm  = deg[col].mean();  ds = deg[col].std()
        nm  = ndeg[col].mean(); ns = ndeg[col].std()
        print(f"  {name:<20s}  {dm:>8.2f} ± {ds:>5.2f}   {nm:>8.2f} ± {ns:>5.2f}")

    # Mann-Whitney U test
    try:
        from scipy.stats import mannwhitneyu
        print(f"\n  Mann-Whitney U p-values (degrader vs non-degrader):")
        for col in ['ang_kin_e2','dist_kin_e2_cen','min_d_kin_e2','frac_e2_50A','frac_e2_80A']:
            u, p = mannwhitneyu(deg[col], ndeg[col], alternative='two-sided')
            print(f"    {col:<25s}  p={p:.4f}")
    except ImportError:
        pass

    print(f"\nSaved:\n  {OUT_DIR}/e2_geometry_all_poses.tsv\n  {OUT_DIR}/e2_geometry_poi_summary.tsv")


if __name__ == '__main__':
    main()
