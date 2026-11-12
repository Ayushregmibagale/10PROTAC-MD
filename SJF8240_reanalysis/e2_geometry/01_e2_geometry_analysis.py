#!/usr/bin/env python3
"""
Structural analysis: why do ABL1/AXL have negative e2_delta_reach_50A?

Three hypotheses tested on pRosettaC rank-1 ternary complex structures:

H1 (Orientation):   angle between VHL→kinase vector and VHL→E2 centroid vector.
                    If kinase points AWAY from E2, E2 can't reach kinase lysines.

H2 (Bulk/shadow):   how many kinase CA atoms fall within 50/80Å of any E2 point?
                    Radius of gyration as a size proxy.
                    Does the kinase cloud physically overlap the E2 cloud?

H3 (Which face):    which part of the kinase (N-lobe / C-lobe / activation loop region)
                    faces the E2 cloud? Fraction of N-lobe vs C-lobe CAs within
                    cutoff of E2.

Run: ~/.conda/envs/gromacs/bin/python3 -u 01_e2_geometry_analysis.py
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
POS_DIR    = SCRIPT_DIR.parent / 'geometric_gnn' / 'positions'
REFDIR     = Path.home() / 'protac-md' / 'reference_structures'
E2_ENS     = REFDIR / 'e2_c85_ensemble.json'
OUT_DIR    = SCRIPT_DIR / 'results'
OUT_DIR.mkdir(parents=True, exist_ok=True)

LABELS = {
    'DDR2':1,'ABL1':0,'AXL':0,'EPHA2':1,
    'EPHB2':1,'MAPK14_full':1,'MET':1,'RIPK2_full':1,
}
POI_ORDER = ['DDR2','ABL1','AXL','EPHA2','EPHB2','MAPK14_full','MET','RIPK2_full']

E2_CUTOFFS = [30, 50, 80]


def load_e2_ensemble():
    with open(E2_ENS) as f:
        data = json.load(f)
    if isinstance(data, list):
        return np.array(data, dtype=np.float32)
    return np.array([p for v in data.values() for p in v], dtype=np.float32)


def angle_deg(v1, v2):
    """Angle between two vectors in degrees."""
    c = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9)
    return float(np.degrees(np.arccos(np.clip(c, -1, 1))))


def min_dist_to_cloud(points, cloud):
    """Minimum distance from any point in `points` to any point in `cloud`."""
    # Compute pairwise distances in chunks to avoid huge memory
    min_d = np.inf
    chunk = 500
    for i in range(0, len(points), chunk):
        d = np.sqrt(((points[i:i+chunk, None] - cloud[None]) ** 2).sum(2))
        min_d = min(min_d, float(d.min()))
    return min_d


def fraction_within(points, cloud, cutoff):
    """Fraction of `points` with any E2 point within `cutoff` Å."""
    count = 0
    chunk = 200
    for i in range(0, len(points), chunk):
        d = np.sqrt(((points[i:i+chunk, None] - cloud[None]) ** 2).sum(2))
        count += int((d.min(1) < cutoff).sum())
    return count / max(len(points), 1)


def radius_of_gyration(points):
    cen = points.mean(0)
    return float(np.sqrt(((points - cen) ** 2).sum(1).mean()))


def main():
    e2 = load_e2_ensemble()
    e2_centroid = e2.mean(0)
    print(f"E2 ensemble: {e2.shape[0]} points, centroid={e2_centroid.round(1)}")

    rows = []

    for poi in POI_ORDER:
        p = POS_DIR / f'{poi}.npz'
        if not p.exists():
            print(f"  [WARN] missing {p}")
            continue

        d = np.load(p)
        vhl_pos = d['vhl_pos'].astype(np.float32)
        kin_pos = d['kin_pos'].astype(np.float32)
        kin_res = d['kin_resids']
        label   = int(d['label'][0])

        vhl_cen = vhl_pos.mean(0)
        kin_cen = kin_pos.mean(0)

        # ── H1: Orientation ───────────────────────────────────────────────────
        v_vhl_e2  = e2_centroid - vhl_cen          # VHL → E2
        v_vhl_kin = kin_cen - vhl_cen              # VHL → kinase
        ang_kin_e2 = angle_deg(v_vhl_kin, v_vhl_e2)

        # Kinase-to-E2 distance (centroid to centroid)
        dist_kin_e2_cen = float(np.linalg.norm(kin_cen - e2_centroid))

        # ── H2: Bulk / shadow ─────────────────────────────────────────────────
        rg_kin = radius_of_gyration(kin_pos)
        n_kin  = len(kin_pos)

        # Min distance from any kinase CA to any E2 point
        min_d_kin_e2 = min_dist_to_cloud(kin_pos, e2)

        # Fraction of kinase CAs reachable from E2 at different cutoffs
        frac = {c: fraction_within(kin_pos, e2, c) for c in E2_CUTOFFS}

        # Does the kinase CLOUD physically overlap E2 cloud?
        # Use convex proxy: fraction of E2 points within Rg+10 of kinase centroid
        e2_near_kin = float((np.linalg.norm(e2 - kin_cen, axis=1) < rg_kin + 10).mean())

        # ── H3: Which face (N-lobe vs C-lobe) ────────────────────────────────
        # N-lobe ≈ first 40% of kinase residues, C-lobe ≈ last 60%
        n_total = len(kin_res)
        split   = int(n_total * 0.4)
        nlobe_pos = kin_pos[:split]
        clobe_pos = kin_pos[split:]

        nlobe_min_d  = min_dist_to_cloud(nlobe_pos, e2)
        clobe_min_d  = min_dist_to_cloud(clobe_pos, e2)
        nlobe_frac50 = fraction_within(nlobe_pos, e2, 50)
        clobe_frac50 = fraction_within(clobe_pos, e2, 50)

        # Which lobe centroid is closer to E2?
        nlobe_cen_d = float(np.linalg.norm(nlobe_pos.mean(0) - e2_centroid))
        clobe_cen_d = float(np.linalg.norm(clobe_pos.mean(0) - e2_centroid))
        closer_lobe = 'N' if nlobe_cen_d < clobe_cen_d else 'C'

        row = {
            'poi':         poi,
            'label':       label,
            # H1
            'ang_kin_e2':      round(ang_kin_e2, 1),
            'dist_kin_e2_cen': round(dist_kin_e2_cen, 1),
            # H2
            'n_kin_res':   n_kin,
            'rg_kin':      round(rg_kin, 1),
            'min_d_kin_e2': round(min_d_kin_e2, 1),
            'frac_e2_30A': round(frac[30], 3),
            'frac_e2_50A': round(frac[50], 3),
            'frac_e2_80A': round(frac[80], 3),
            'e2_near_kin': round(e2_near_kin, 3),
            # H3
            'nlobe_min_d': round(nlobe_min_d, 1),
            'clobe_min_d': round(clobe_min_d, 1),
            'nlobe_frac50': round(nlobe_frac50, 3),
            'clobe_frac50': round(clobe_frac50, 3),
            'nlobe_cen_d': round(nlobe_cen_d, 1),
            'clobe_cen_d': round(clobe_cen_d, 1),
            'closer_lobe': closer_lobe,
        }
        rows.append(row)

        status = "DEGRADE" if label == 1 else "NO-DEG "
        print(f"\n[{poi}] {status}  n_kin={n_kin}  Rg={rg_kin:.1f}")
        print(f"  H1: ang(VHL→kin, VHL→E2) = {ang_kin_e2:.1f}°  |  dist(kin_cen, E2_cen) = {dist_kin_e2_cen:.1f} Å")
        print(f"  H2: min_d(kin→E2)={min_d_kin_e2:.1f} Å  |  frac@50A={frac[50]:.3f}  frac@80A={frac[80]:.3f}")
        print(f"  H3: closer_lobe={closer_lobe}  N_cen={nlobe_cen_d:.1f}  C_cen={clobe_cen_d:.1f}  "
              f"N_frac50={nlobe_frac50:.3f}  C_frac50={clobe_frac50:.3f}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / 'e2_geometry.tsv', sep='\t', index=False)
    print(f"\n{'='*70}")
    print("SUMMARY — sorted by ang_kin_e2")
    print('='*70)
    print(df[['poi','label','ang_kin_e2','dist_kin_e2_cen','n_kin_res','rg_kin',
              'min_d_kin_e2','frac_e2_50A','closer_lobe']].sort_values('ang_kin_e2').to_string(index=False))

    # Print degrader vs non-degrader means
    print(f"\n{'':20s}  {'degrader':>12}  {'non-deg':>12}")
    deg  = df[df['label']==1]
    ndeg = df[df['label']==0]
    for col in ['ang_kin_e2','dist_kin_e2_cen','min_d_kin_e2','frac_e2_50A','rg_kin']:
        print(f"  {col:<20s}  {deg[col].mean():>12.2f}  {ndeg[col].mean():>12.2f}")

    print(f"\nSaved: {OUT_DIR / 'e2_geometry.tsv'}")


if __name__ == '__main__':
    main()
