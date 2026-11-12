#!/usr/bin/env python3
"""
Build two new graph types from MD trajectories:

A) Per-pose contact graphs
   One contact frequency graph per individual 20ns MD pose.
   Gives 16-19 graphs per POI instead of 1.
   Saves: geometric_gnn/graphs/per_pose/<POI>_<pose_id>.npz

B) Bridging-conditional contact graphs
   For each POI/pose, splits trajectory into:
     G_bridge   — frames where PROTAC contacts BOTH VHL (<5Å) AND kinase (<5Å)
     G_nonbr    — all other frames
   Then per-POI aggregation: mean over poses
   Saves: geometric_gnn/graphs/bridge/<POI>_bridge.npz
           geometric_gnn/graphs/bridge/<POI>_nonbr.npz

Run with: ~/.conda/envs/gromacs/bin/python3 04_build_per_pose_and_bridge_graphs.py

Requires: MDAnalysis (available in gromacs env)
"""
import sys, argparse
import numpy as np
import pandas as pd
from pathlib import Path

try:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import capped_distance
except ImportError:
    print("[ERROR] MDAnalysis not found. Run: ~/.conda/envs/gromacs/bin/python3 ...")
    sys.exit(1)

RANDOM_SCREEN = Path.home() / 'protac-md' / 'random_screen'
MANIFEST      = RANDOM_SCREEN / 'manifest.tsv'
OUT_DIR       = Path(__file__).resolve().parent / 'graphs'

LABELS = {
    'DDR2': 1, 'ABL1': 0, 'AXL': 0, 'EPHA2': 1,
    'EPHB2': 1, 'MAPK14_full': 1, 'MET': 1, 'RIPK2_full': 1
}
POI_KIN_MIN = {
    'DDR2': 553, 'ABL1': 235, 'AXL': 521, 'EPHA2': 600,
    'EPHB2': 615, 'MAPK14_full': 200, 'MET': 1053, 'RIPK2_full': 200
}

CONTACT_THRESH  = 5.0    # Å
BRIDGE_THRESH   = 5.0    # Å — same threshold for bridging definition
EDGE_THRESH     = 0.01   # include edge if contacted in ≥ 1% of frames
STRIDE          = 10     # subsampling


# ── Topology helpers (copied from 04c) ───────────────────────────────────────

def _find_chain_boundaries(gro_path):
    lines = open(gro_path).readlines()
    n_atoms = int(lines[1])
    atoms = lines[2: 2 + n_atoms]
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
    for j in range(len(boundaries) - 1):
        s, e = boundaries[j], boundaries[j+1]
        chains.append({'start': s+1, 'end': e,
                       'first_res': int(atoms[s][0:5]),
                       'last_res':  int(atoms[e-1][0:5])})
    return chains


def _detect_atom_ranges(gro_path):
    chains = _find_chain_boundaries(gro_path)
    try:
        vhl = next(c for c in chains if c['first_res'] == 1)
    except StopIteration:
        return None
    kin = next((c for c in chains if c['first_res'] > 103), None)
    if kin is None:
        return None
    return (vhl['start'], vhl['end']), (kin['start'], kin['end'])


# ── Core contact computation (single trajectory) ──────────────────────────────

def compute_contacts(tpr, xtc, gro, stride=STRIDE, bridge_split=False):
    """
    Compute contact frequency between PROTAC and each VHL/kinase CA.

    If bridge_split=True, return separate arrays for bridging and non-bridging frames.

    Returns dict or None on failure.
    """
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

    if len(protac) == 0 or len(vhl_ca) == 0 or len(kin_ca) == 0:
        return None

    vhl_resids = np.array([a.resid for a in vhl_ca])
    kin_resids = np.array([a.resid for a in kin_ca])
    n_vhl, n_kin = len(vhl_ca), len(kin_ca)

    # Arrays for overall, bridging, non-bridging
    def _zero(): return (np.zeros(n_vhl), np.zeros(n_kin),
                         np.zeros(n_vhl), np.zeros(n_kin))
    vc_all, kc_all, vd_all, kd_all = _zero()
    vc_br,  kc_br,  vd_br,  kd_br  = _zero()
    vc_nb,  kc_nb,  vd_nb,  kd_nb  = _zero()
    n_all = n_br = n_nb = 0

    for ts in u.trajectory[::stride]:
        pp = protac.positions
        vp = vhl_ca.positions
        kp = kin_ca.positions

        dv = np.sqrt(((pp[:, None] - vp[None]) ** 2).sum(2)).min(0)
        dk = np.sqrt(((pp[:, None] - kp[None]) ** 2).sum(2)).min(0)

        vc_all += dv < CONTACT_THRESH
        kc_all += dk < CONTACT_THRESH
        vd_all += dv
        kd_all += dk
        n_all  += 1

        if bridge_split:
            is_bridging = (dv.min() < BRIDGE_THRESH) and (dk.min() < BRIDGE_THRESH)
            if is_bridging:
                vc_br += dv < CONTACT_THRESH
                kc_br += dk < CONTACT_THRESH
                vd_br += dv
                kd_br += dk
                n_br  += 1
            else:
                vc_nb += dv < CONTACT_THRESH
                kc_nb += dk < CONTACT_THRESH
                vd_nb += dv
                kd_nb += dk
                n_nb  += 1

    if n_all == 0:
        return None

    def _norm(vc, kc, vd, kd, n):
        return (vc / n, kc / n, vd / (n * 10), kd / (n * 10))

    result = {
        'vhl_resids': vhl_resids, 'kin_resids': kin_resids,
        'all': _norm(vc_all, kc_all, vd_all, kd_all, n_all),
        'n_all': n_all,
    }
    if bridge_split:
        result['bridge'] = _norm(vc_br, kc_br, vd_br, kd_br, max(n_br,1)) if n_br else None
        result['nonbr']  = _norm(vc_nb, kc_nb, vd_nb, kd_nb, max(n_nb,1)) if n_nb else None
        result['n_br'] = n_br
        result['n_nb'] = n_nb

    return result


# ── Graph builder ─────────────────────────────────────────────────────────────

def build_graph_npz(vhl_resids, kin_resids, vhl_freq, kin_freq,
                    vhl_d, kin_d, label):
    """Build edge list and return dict ready for np.savez."""
    n_vhl = len(vhl_resids)
    n_kin = len(kin_resids)
    n_total = n_vhl + n_kin

    # Node features [is_vhl, is_kin, seq_pos_norm]
    node_feats = np.zeros((n_total, 3), dtype=np.float32)
    node_feats[:n_vhl, 0] = 1.0
    node_feats[n_vhl:, 1] = 1.0
    vp = (vhl_resids - vhl_resids.min()) / max(vhl_resids.max() - vhl_resids.min(), 1)
    kp = (kin_resids - kin_resids.min()) / max(kin_resids.max() - kin_resids.min(), 1)
    node_feats[:n_vhl, 2] = vp.astype(np.float32)
    node_feats[n_vhl:, 2] = kp.astype(np.float32)

    src, dst, ef, edv, edk = [], [], [], [], []
    for i in range(n_vhl):
        for j in range(n_kin):
            co = vhl_freq[i] * kin_freq[j]
            if co > EDGE_THRESH ** 2:
                vn, kn = i, n_vhl + j
                src += [vn, kn]; dst += [kn, vn]
                ef  += [co, co]
                edv += [vhl_d[i], vhl_d[i]]
                edk += [kin_d[j], kin_d[j]]

    if not src:
        # Sparse graph — force at least one edge using max-contact pair
        bi = int(np.argmax(vhl_freq))
        bj = int(np.argmax(kin_freq))
        vn, kn = bi, n_vhl + bj
        co = max(vhl_freq[bi] * kin_freq[bj], 1e-6)
        src += [vn, kn]; dst += [kn, vn]
        ef += [co, co]
        edv += [vhl_d[bi], vhl_d[bi]]
        edk += [kin_d[bj], kin_d[bj]]

    ea = np.column_stack([ef, edv, edk]).astype(np.float32)
    return dict(
        x=node_feats,
        edge_index=np.array([src, dst], dtype=np.int64),
        edge_attr=ea,
        label=np.array([label], dtype=np.int64),
    )


# ── Per-POI aggregation ───────────────────────────────────────────────────────

def process_poi(poi, done_poses, build_bridge=True):
    kin_min  = POI_KIN_MIN[poi]
    label    = LABELS[poi]
    per_pose_all    = []
    per_pose_bridge = []
    per_pose_nonbr  = []
    pose_meta = []

    for pose in done_poses:
        pid     = pose['pose_id']
        pdir    = RANDOM_SCREEN / poi / pid
        tpr     = pdir / 'prod/prod.tpr'
        xtc     = pdir / 'prod/prod.xtc'
        gro     = pdir / 'topol/complex_protein.gro'

        print(f"    {pid} ...", end=' ', flush=True)
        result = compute_contacts(tpr, xtc, gro, bridge_split=build_bridge)
        if result is None:
            print("FAILED")
            continue

        vhl_r = result['vhl_resids']
        kin_r = result['kin_resids']
        vf, kf, vd, kd = result['all']

        per_pose_all.append((vf, kf, vd, kd))
        pose_meta.append({'poi': poi, 'pose_id': pid,
                          'n_frames': result['n_all'],
                          'vhl_max': float(vf.max()), 'kin_max': float(kf.max())})

        if build_bridge:
            if result['bridge']:
                per_pose_bridge.append(result['bridge'])
            else:
                per_pose_bridge.append(None)

            if result['nonbr']:
                per_pose_nonbr.append(result['nonbr'])
            else:
                per_pose_nonbr.append(None)

        print(f"n={result['n_all']} vhl={vf.max():.3f} kin={kf.max():.3f}"
              + (f" br_frames={result.get('n_br',0)}" if build_bridge else ""),
              flush=True)

    if not per_pose_all:
        return None

    return {
        'poi': poi, 'label': label,
        'vhl_resids': vhl_r, 'kin_resids': kin_r,
        'per_pose_all': per_pose_all,
        'per_pose_bridge': per_pose_bridge if build_bridge else None,
        'per_pose_nonbr':  per_pose_nonbr  if build_bridge else None,
        'pose_meta': pose_meta,
    }


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--poi', nargs='+', default=list(POI_KIN_MIN.keys()))
    parser.add_argument('--no-bridge', action='store_true')
    args = parser.parse_args()

    manifest = pd.read_csv(MANIFEST, sep='\t')
    done     = manifest[manifest['stage'] == 'ANALYSIS_DONE']

    per_pose_dir = OUT_DIR / 'per_pose'
    bridge_dir   = OUT_DIR / 'bridge'
    per_pose_dir.mkdir(parents=True, exist_ok=True)
    bridge_dir.mkdir(parents=True, exist_ok=True)

    all_meta = []

    for poi in args.poi:
        if poi not in LABELS:
            print(f"[{poi}] not in label map — skipping")
            continue

        poi_done = done[done['poi'] == poi].to_dict('records')
        print(f"\n[{poi}] label={LABELS[poi]}  poses={len(poi_done)}", flush=True)

        result = process_poi(poi, poi_done, build_bridge=not args.no_bridge)
        if result is None:
            print(f"  [{poi}] no valid poses")
            continue

        vhl_r = result['vhl_resids']
        kin_r = result['kin_resids']
        label = result['label']
        all_meta.extend(result['pose_meta'])

        # ── A: Per-pose individual graphs ──
        for i, (vf, kf, vd, kd) in enumerate(result['per_pose_all']):
            pid = result['pose_meta'][i]['pose_id']
            g = build_graph_npz(vhl_r, kin_r, vf, kf, vd, kd, label)
            np.savez(per_pose_dir / f'{poi}_{pid}.npz', **g)

        print(f"  Saved {len(result['per_pose_all'])} per-pose graphs")

        # ── B: Bridging-conditional aggregated graphs ──
        if not args.no_bridge:
            bridge_poses = [b for b in result['per_pose_bridge'] if b is not None]
            nonbr_poses  = [b for b in result['per_pose_nonbr']  if b is not None]

            if bridge_poses:
                vf = np.mean([b[0] for b in bridge_poses], 0)
                kf = np.mean([b[1] for b in bridge_poses], 0)
                vd = np.mean([b[2] for b in bridge_poses], 0)
                kd = np.mean([b[3] for b in bridge_poses], 0)
                g_br = build_graph_npz(vhl_r, kin_r, vf, kf, vd, kd, label)
                np.savez(bridge_dir / f'{poi}_bridge.npz', **g_br)
                print(f"  Bridge graph: n_edges={g_br['edge_index'].shape[1]//2}  "
                      f"vhl_max={vf.max():.3f}  kin_max={kf.max():.3f}")
            else:
                print(f"  [WARN] No bridging frames for {poi}")

            if nonbr_poses:
                vf = np.mean([b[0] for b in nonbr_poses], 0)
                kf = np.mean([b[1] for b in nonbr_poses], 0)
                vd = np.mean([b[2] for b in nonbr_poses], 0)
                kd = np.mean([b[3] for b in nonbr_poses], 0)
                g_nb = build_graph_npz(vhl_r, kin_r, vf, kf, vd, kd, label)
                np.savez(bridge_dir / f'{poi}_nonbr.npz', **g_nb)
                print(f"  NonBr graph:  n_edges={g_nb['edge_index'].shape[1]//2}  "
                      f"vhl_max={vf.max():.3f}  kin_max={kf.max():.3f}")

    pd.DataFrame(all_meta).to_csv(OUT_DIR / 'per_pose_meta.tsv', sep='\t', index=False)
    print(f"\nDone. Meta: {OUT_DIR / 'per_pose_meta.tsv'}")


if __name__ == '__main__':
    main()
