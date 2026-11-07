#!/usr/bin/env python3
"""
Tier2b MD branch: Build MD-derived contact frequency graphs.

For each POI, aggregates all ANALYSIS_DONE random-screen production trajectories
(20 ns each, ~2000 frames) into a per-POI contact frequency graph:

  - Nodes  = VHL residues (Cα) + kinase residues (Cα)  [no separate PROTAC nodes]
  - Edges  = PROTAC-to-residue contact frequency (min distance < CONTACT_THRESH)
  - Edge weight = fraction of trajectory frames where any PROTAC heavy atom
                  is within CONTACT_THRESH Å of the residue's heavy atoms
  - Graph label  = 1 (degrader) / 0 (non-degrader)

This is the DYNAMIC complement to the pRosettaC ensemble graph (04a):
  pRosettaC: can the geometry be achieved? (static)
  MD:        is it thermodynamically stable when achieved? (dynamic)

Outputs:
  analysis/gnn/md_contact_graphs.pt       — list of PyTorch Geometric Data objects
  analysis/gnn/md_contact_graphs_np.npz   — numpy fallback
  analysis/gnn/md_graph_metadata.tsv      — per-POI graph stats
  analysis/gnn/md_contact_per_pose.tsv    — per-pose VHL/kinase contact frequencies

Node features (6 dims): [is_vhl, is_kin, seq_pos_norm, aa_onehot_C, aa_onehot_N, mean_b_factor]
Edge features (3 dims): [contact_freq, mean_dist_nm, std_dist_nm]

Note: PROTAC resname in GROMACS is 'UNL' (acpype default).
VHL/kinase ranges detected from complex_protein.gro residue number jumps.

Usage:
  python3 scripts/02_analysis/04c_build_md_contact_graphs.py
  python3 scripts/02_analysis/04c_build_md_contact_graphs.py --poi DDR2 ABL1 --stride 5
"""
import argparse
import sys
import numpy as np
import pandas as pd
from pathlib import Path

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances as mdadist
    HAS_MDA = True
except ImportError:
    HAS_MDA = False

try:
    import torch
    from torch_geometric.data import Data
    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False

RANDOM_SCREEN = Path.home() / 'protac-md' / 'random_screen'
ANALYSIS_DIR  = Path.home() / 'protac-md' / 'analysis'
GNN_DIR       = ANALYSIS_DIR / 'gnn'
MANIFEST      = RANDOM_SCREEN / 'manifest.tsv'
OUT_PT        = GNN_DIR / 'md_contact_graphs.pt'
OUT_NPZ       = GNN_DIR / 'md_contact_graphs_np.npz'
OUT_META      = GNN_DIR / 'md_graph_metadata.tsv'
OUT_POSE      = GNN_DIR / 'md_contact_per_pose.tsv'

VHL_MAX       = 103     # VHL always residues 1-103
CONTACT_THRESH = 5.0    # Å — residue is "contacted" if any PROTAC atom within this
EDGE_FREQ_THRESH = 0.01 # include edge if contacted in ≥ 1% of frames
STRIDE        = 10      # subsample every N frames (20ns/2ps step = 10000 frames; stride 10 = 1000 frames)

KINOME_SCREEN = Path.home() / 'protac-md' / 'kinome_md_screen'
KINOME_MANIFEST = KINOME_SCREEN / 'manifest.tsv'

POI_CONFIG = {
    'DDR2':        {'kin_min': 553,  'degrader': True},
    'RIPK2_full':  {'kin_min': 200,  'degrader': True},
    'MAPK14_full': {'kin_min': 200,  'degrader': True},
    'MET':         {'kin_min': 1053, 'degrader': True},
    'EPHB2':       {'kin_min': 615,  'degrader': True},
    'ABL1':        {'kin_min': 235,  'degrader': False},
    'AXL':         {'kin_min': 521,  'degrader': False},
    'EPHA2':       {'kin_min': 600,  'degrader': False},
    # Kinome extension — prediction targets (label unknown → degrader=None)
    # graph_name must match the pRosettaC ensemble graph poi name used in 04a/04d
    'SIK3':    {'kin_min': 200, 'degrader': None, 'graph_name': 'SIK3_full'},
    'RPS6KB2': {'kin_min': 200, 'degrader': None, 'graph_name': 'RPS6KB2_full'},
}


# ── Topology detection (reused from 02_run_analysis.py) ───────────────────────

def _find_chain_boundaries(gro_path: str) -> list[dict]:
    lines = open(gro_path).readlines()
    n_atoms = int(lines[1])
    atoms   = lines[2: 2 + n_atoms]
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
        s, e = boundaries[j], boundaries[j + 1]
        chains.append({
            'start':     s + 1, 'end': e,
            'first_res': int(atoms[s][0:5]),
            'last_res':  int(atoms[e - 1][0:5]),
        })
    return chains, atoms


def _detect_vhl_kin_resid(gro_path: str, kin_min: int) -> tuple[tuple[int, int], tuple[int, int]]:
    """Return (vhl_first, vhl_last), (kin_first, kin_last) residue numbers from GRO."""
    chains, atoms = _find_chain_boundaries(gro_path)
    vhl = next(c for c in chains if c['first_res'] == 1)

    try:
        kin = next(c for c in chains if c['first_res'] > 103)
    except StopIteration:
        # kin_min=104 case: no residue jump; scan for first res=104
        kin = next((c for c in chains if c['first_res'] <= 104 <= c['last_res']), None)
        if kin is None:
            kin = {'first_res': kin_min, 'last_res': max(c['last_res'] for c in chains)}
        else:
            kin = {**kin, 'first_res': 104}
        vhl = {**vhl, 'last_res': 103}

    return (vhl['first_res'], vhl['last_res']), (kin['first_res'], kin['last_res'])


def _detect_vhl_kin_atom_ranges(gro_path: str) -> tuple[tuple[int, int], tuple[int, int]] | None:
    """Return (vhl_start, vhl_end), (kin_start, kin_end) as 1-indexed atom index ranges.

    Uses GRO atom positions directly, bypassing MDAnalysis resid renumbering.
    Works for all POIs including those with large kin_min (MET=1053, EPHB2=615)
    where MDAnalysis sequential resids don't match GRO residue numbers.
    Returns None if chains cannot be identified.
    """
    chains, _ = _find_chain_boundaries(gro_path)
    try:
        vhl_chain = next(c for c in chains if c['first_res'] == 1)
    except StopIteration:
        return None

    # Kinase: first chain with first_res > 103 (has a residue-number jump from VHL)
    kin_chain = next((c for c in chains if c['first_res'] > 103), None)
    if kin_chain is None:
        # kin_min=104 case: no jump; kinase continues after VHL in single chain
        # GRO resids are sequential so MDAnalysis resid selection still works — caller falls back
        return None

    return (vhl_chain['start'], vhl_chain['end']), (kin_chain['start'], kin_chain['end'])


# ── Per-pose contact frequency computation ────────────────────────────────────

def compute_pose_contacts(
    tpr: Path, xtc: Path, gro: Path, kin_min: int, stride: int = STRIDE
) -> dict | None:
    """
    Compute per-frame contact frequency between PROTAC and each VHL/kinase residue.
    Returns dict: {
        'vhl_resids':  np.array,         # VHL Cα residue numbers
        'kin_resids':  np.array,         # kinase Cα residue numbers
        'vhl_freq':    np.array,         # contact frequency per VHL residue
        'kin_freq':    np.array,         # contact frequency per kinase residue
        'vhl_mean_d':  np.array,         # mean distance per VHL residue (Å)
        'kin_mean_d':  np.array,
        'n_frames':    int,
    }
    """
    if not tpr.exists() or not xtc.exists() or not gro.exists():
        return None

    try:
        atom_ranges = _detect_vhl_kin_atom_ranges(str(gro))
    except Exception as e:
        print(f"      topology detection failed: {e}")
        return None

    if atom_ranges is None:
        # kin_min=104 fallback: GRO resids are sequential, use resid selection
        try:
            (vhl_r0, vhl_r1), (kin_r0, kin_r1) = _detect_vhl_kin_resid(str(gro), kin_min)
        except Exception as e:
            print(f"      topology detection (resid fallback) failed: {e}")
            return None
        resid_mode = True
    else:
        (vhl_a0, vhl_a1), (kin_a0, kin_a1) = atom_ranges
        resid_mode = False

    try:
        u = mda.Universe(str(tpr), str(xtc))
    except Exception as e:
        print(f"      MDAnalysis load failed: {e}")
        return None

    # Select atom groups — bynum uses 1-indexed GRO atom positions, immune to MDAnalysis resid renumbering
    protac = u.select_atoms("resname UNL and not name H*")
    if resid_mode:
        vhl_ca = u.select_atoms(f"protein and resid {vhl_r0}:{vhl_r1} and name CA")
        kin_ca = u.select_atoms(f"protein and resid {kin_r0}:{kin_r1} and name CA")
    else:
        vhl_ca = u.select_atoms(f"bynum {vhl_a0}:{vhl_a1} and name CA")
        kin_ca = u.select_atoms(f"bynum {kin_a0}:{kin_a1} and name CA")

    if len(protac) == 0:
        print(f"      PROTAC (resname UNL) not found")
        return None
    if len(vhl_ca) == 0 or len(kin_ca) == 0:
        print(f"      VHL ({len(vhl_ca)}) or kinase ({len(kin_ca)}) CA not found")
        return None

    vhl_resids = np.array([a.resid for a in vhl_ca])
    kin_resids = np.array([a.resid for a in kin_ca])
    n_vhl = len(vhl_ca)
    n_kin = len(kin_ca)

    # For contact: use all heavy atoms of each residue vs all PROTAC heavy atoms
    # Faster: use Cα only for speed, acceptable for residue-level contact frequency
    vhl_contact_count  = np.zeros(n_vhl)
    kin_contact_count  = np.zeros(n_kin)
    vhl_dist_sum       = np.zeros(n_vhl)
    kin_dist_sum       = np.zeros(n_kin)
    n_frames = 0

    for ts in u.trajectory[::stride]:
        ptc_pos = protac.positions   # (n_ptc, 3)
        vhl_pos = vhl_ca.positions   # (n_vhl, 3)
        kin_pos = kin_ca.positions   # (n_kin, 3)

        # Min PROTAC–VHL distance per VHL residue
        dv = np.sqrt(((ptc_pos[:, None, :] - vhl_pos[None, :, :]) ** 2).sum(axis=2))  # (n_ptc, n_vhl)
        min_dv = dv.min(axis=0)  # (n_vhl,)

        dk = np.sqrt(((ptc_pos[:, None, :] - kin_pos[None, :, :]) ** 2).sum(axis=2))  # (n_ptc, n_kin)
        min_dk = dk.min(axis=0)  # (n_kin,)

        vhl_contact_count  += (min_dv < CONTACT_THRESH)
        kin_contact_count  += (min_dk < CONTACT_THRESH)
        vhl_dist_sum       += min_dv
        kin_dist_sum       += min_dk
        n_frames += 1

    if n_frames == 0:
        return None

    return {
        'vhl_resids':  vhl_resids,
        'kin_resids':  kin_resids,
        'vhl_freq':    vhl_contact_count / n_frames,
        'kin_freq':    kin_contact_count / n_frames,
        'vhl_mean_d':  vhl_dist_sum / n_frames,
        'kin_mean_d':  kin_dist_sum / n_frames,
        'n_frames':    n_frames,
        'n_ptc_atoms': len(protac),
    }


# ── POI-level aggregation ─────────────────────────────────────────────────────

def build_poi_md_graph(poi: str, done_poses: list[dict], stride: int, cfg: dict) -> dict | None:
    kin_min = cfg['kin_min']

    pose_results = []
    pose_meta    = []

    screen_dir = cfg.get('screen_dir', RANDOM_SCREEN)
    for pose in done_poses:
        pose_id = pose['pose_id']
        pose_dir = screen_dir / poi / pose_id

        tpr = pose_dir / 'prod' / 'prod.tpr'
        xtc = pose_dir / 'prod' / 'prod.xtc'
        gro = pose_dir / 'topol' / 'complex_protein.gro'

        print(f"    {pose_id} ...", end=' ', flush=True)
        result = compute_pose_contacts(tpr, xtc, gro, kin_min, stride)

        if result is None:
            print("FAILED")
            continue

        print(f"n_frames={result['n_frames']}  "
              f"vhl_max_freq={result['vhl_freq'].max():.3f}  "
              f"kin_max_freq={result['kin_freq'].max():.3f}", flush=True)

        pose_results.append(result)
        pose_meta.append({
            'poi': poi, 'pose_id': pose_id,
            'n_frames': result['n_frames'],
            'vhl_max_freq': float(result['vhl_freq'].max()),
            'kin_max_freq': float(result['kin_freq'].max()),
            'vhl_mean_freq': float(result['vhl_freq'].mean()),
            'kin_mean_freq': float(result['kin_freq'].mean()),
        })

    if not pose_results:
        return None

    # Use canonical residue lists from first valid result
    vhl_resids = pose_results[0]['vhl_resids']
    kin_resids = pose_results[0]['kin_resids']
    n_vhl = len(vhl_resids)
    n_kin = len(kin_resids)

    # Average across poses
    vhl_freq_mean = np.mean([r['vhl_freq'] for r in pose_results], axis=0)
    kin_freq_mean = np.mean([r['kin_freq'] for r in pose_results], axis=0)
    vhl_d_mean    = np.mean([r['vhl_mean_d'] for r in pose_results], axis=0)
    kin_d_mean    = np.mean([r['kin_mean_d'] for r in pose_results], axis=0)
    vhl_d_std     = np.std([r['vhl_freq'] for r in pose_results], axis=0)
    kin_d_std     = np.std([r['kin_freq'] for r in pose_results], axis=0)

    # ── Build node features ──
    # [is_vhl, is_kin, seq_pos_norm]
    n_total = n_vhl + n_kin
    node_feats = np.zeros((n_total, 3), dtype=np.float32)
    node_feats[:n_vhl, 0] = 1.0   # is_vhl
    node_feats[n_vhl:, 1] = 1.0   # is_kin

    vhl_pos_norm = (vhl_resids - vhl_resids.min()) / max(vhl_resids.max() - vhl_resids.min(), 1)
    kin_pos_norm = (kin_resids - kin_resids.min()) / max(kin_resids.max() - kin_resids.min(), 1)
    node_feats[:n_vhl, 2] = vhl_pos_norm.astype(np.float32)
    node_feats[n_vhl:, 2] = kin_pos_norm.astype(np.float32)

    # ── Build edges: VHL-node and kinase-node connected to a virtual PROTAC node ──
    # OR: connect VHL nodes to kinase nodes weighted by co-contact frequency
    # (both contacted in same frame)
    #
    # Simpler: per-node PROTAC contact frequency as a node attribute, then
    # VHL→kinase edges weighted by product of their contact frequencies
    # (high PROTAC-VHL × high PROTAC-kinase = good bridging pair)
    #
    # Graph structure: VHL nodes connected to kinase nodes with edge weight
    # = product of their PROTAC contact frequencies (bridging co-occurrence proxy)

    edge_src, edge_dst, edge_freq, edge_dvhl, edge_dkin = [], [], [], [], []

    for i in range(n_vhl):
        for j in range(n_kin):
            fv = vhl_freq_mean[i]
            fk = kin_freq_mean[j]
            co_freq = fv * fk   # proxy for simultaneous PROTAC contact
            if co_freq > EDGE_FREQ_THRESH ** 2:   # threshold both
                vhl_node = i
                kin_node = n_vhl + j
                edge_src.extend([vhl_node, kin_node])
                edge_dst.extend([kin_node, vhl_node])
                edge_freq.extend([co_freq, co_freq])
                edge_dvhl.extend([vhl_d_mean[i] / 10.0, vhl_d_mean[i] / 10.0])   # to nm
                edge_dkin.extend([kin_d_mean[j] / 10.0, kin_d_mean[j] / 10.0])

    n_edges = len(edge_src) // 2

    print(f"  [{poi}] {len(pose_results)} poses | "
          f"n_nodes={n_total}  n_edges={n_edges}  "
          f"vhl_max={vhl_freq_mean.max():.3f}  kin_max={kin_freq_mean.max():.3f}", flush=True)

    graph_dict = {
        'poi':          cfg.get('graph_name', poi),
        'degrader':     int(cfg['degrader']) if cfg['degrader'] is not None else -1,
        'n_vhl':        n_vhl,
        'n_kin':        n_kin,
        'n_poses':      len(pose_results),
        'n_edges':      n_edges,
        'node_feats':   node_feats,
        'vhl_resids':   vhl_resids,
        'kin_resids':   kin_resids,
        'vhl_freq_mean': vhl_freq_mean.astype(np.float32),
        'kin_freq_mean': kin_freq_mean.astype(np.float32),
        'vhl_d_mean':   (vhl_d_mean / 10.0).astype(np.float32),   # nm
        'kin_d_mean':   (kin_d_mean / 10.0).astype(np.float32),
        'edge_src':     np.array(edge_src,   dtype=np.int64),
        'edge_dst':     np.array(edge_dst,   dtype=np.int64),
        'edge_freq':    np.array(edge_freq,  dtype=np.float32),
        'pose_meta':    pose_meta,
    }

    if HAS_TORCH and edge_src:
        ea = np.column_stack([edge_freq, edge_dvhl, edge_dkin])
        graph_dict['pyg'] = Data(
            x          = torch.tensor(node_feats, dtype=torch.float32),
            edge_index = torch.tensor([edge_src, edge_dst], dtype=torch.long),
            edge_attr  = torch.tensor(ea, dtype=torch.float32),
            y          = torch.tensor([int(cfg['degrader']) if cfg['degrader'] is not None else -1], dtype=torch.long),
            poi        = cfg.get('graph_name', poi),
            n_vhl      = n_vhl,
            n_kin      = n_kin,
            n_poses    = len(pose_results),
        )

    return graph_dict, pose_meta


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    if not HAS_MDA:
        print("[ERROR] MDAnalysis not available. Install: conda install mdanalysis")
        sys.exit(1)

    parser = argparse.ArgumentParser()
    parser.add_argument('--poi',    nargs='+', default=None)
    parser.add_argument('--stride', type=int,  default=STRIDE,
                        help=f'Frame subsampling stride (default: {STRIDE})')
    parser.add_argument('--merge', action='store_true',
                        help='Load existing output files and append new POIs instead of overwriting')
    args = parser.parse_args()

    if not MANIFEST.exists():
        print(f"[ERROR] manifest.tsv not found at {MANIFEST}")
        sys.exit(1)

    manifest = pd.read_csv(MANIFEST, sep='\t')
    done = manifest[manifest['stage'] == 'ANALYSIS_DONE']

    pois = args.poi or list(POI_CONFIG.keys())
    GNN_DIR.mkdir(parents=True, exist_ok=True)

    print(f"\nMD Contact Graph Builder — stride={args.stride}  "
          f"contact_thresh={CONTACT_THRESH}Å  edge_freq_thresh={EDGE_FREQ_THRESH}")
    print(f"POIs: {pois}\n")

    all_graphs  = []
    pyg_graphs  = []
    meta_rows   = []
    all_pose_meta = []

    # --merge: load existing outputs so we append rather than overwrite
    if args.merge:
        existing_pois = set()
        if HAS_TORCH and OUT_PT.exists():
            existing = torch.load(OUT_PT, map_location='cpu', weights_only=False)
            for eg in existing:
                if eg.poi not in existing_pois:
                    pyg_graphs.append(eg)
                    existing_pois.add(eg.poi)
            print(f"[merge] Loaded {len(existing)} existing PyG graphs: {sorted(existing_pois)}")
        if OUT_META.exists():
            meta_rows = [r for r in pd.read_csv(OUT_META, sep='\t').to_dict('records')
                         if r['poi'] in existing_pois]
        pois = [p for p in pois if p not in existing_pois]
        print(f"[merge] New POIs to build: {pois}\n")

    for poi in pois:
        if poi not in POI_CONFIG:
            print(f"[{poi}] not in POI_CONFIG — skipping")
            continue

        poi_done = done[done['poi'] == poi].to_dict('records')
        print(f"\n[{poi}] {len(poi_done)} ANALYSIS_DONE poses:", flush=True)

        if not poi_done:
            print(f"  no completed poses — skipping")
            continue

        result = build_poi_md_graph(poi, poi_done, args.stride, POI_CONFIG[poi])
        if result is None:
            continue

        graph_dict, pose_meta = result
        all_graphs.append(graph_dict)
        all_pose_meta.extend(pose_meta)

        if 'pyg' in graph_dict:
            pyg_graphs.append(graph_dict['pyg'])

        meta_rows.append({
            'poi':          graph_dict['poi'],
            'degrader':     graph_dict['degrader'],
            'n_vhl':        graph_dict['n_vhl'],
            'n_kin':        graph_dict['n_kin'],
            'n_nodes':      graph_dict['n_vhl'] + graph_dict['n_kin'],
            'n_edges':      graph_dict['n_edges'],
            'n_poses':      graph_dict['n_poses'],
            'vhl_max_freq': float(graph_dict['vhl_freq_mean'].max()),
            'kin_max_freq': float(graph_dict['kin_freq_mean'].max()),
        })

    if not all_graphs:
        print("No graphs built.")
        return

    # Save outputs
    pd.DataFrame(meta_rows).to_csv(OUT_META, sep='\t', index=False)
    pd.DataFrame(all_pose_meta).to_csv(OUT_POSE, sep='\t', index=False)
    print(f"\nMetadata saved: {OUT_META}")
    print(f"Per-pose data saved: {OUT_POSE}")

    if HAS_TORCH and pyg_graphs:
        torch.save(pyg_graphs, OUT_PT)
        print(f"PyG graphs saved: {OUT_PT}  ({len(pyg_graphs)} graphs)")

    np_dict = {}
    for g in all_graphs:
        p = g['poi']
        for k in ('node_feats', 'vhl_resids', 'kin_resids', 'vhl_freq_mean',
                  'kin_freq_mean', 'vhl_d_mean', 'kin_d_mean',
                  'edge_src', 'edge_dst', 'edge_freq'):
            if k in g:
                np_dict[f"{p}/{k}"] = g[k]
        np_dict[f"{p}/degrader"] = np.array([g['degrader']])
    np.savez_compressed(OUT_NPZ, **np_dict)
    print(f"NumPy graphs saved: {OUT_NPZ}")

    # Summary table
    print(f"\n{'='*70}")
    print(f"MD CONTACT GRAPH SUMMARY")
    print(f"{'='*70}")
    print(f"\n{'POI':<16} {'label':<7} {'n_poses':>8} {'n_nodes':>8} "
          f"{'n_edges':>8} {'vhl_max':>8} {'kin_max':>8}")
    print('-' * 60)
    for r in meta_rows:
        label = 'DEG' if r['degrader'] else 'NON'
        print(f"{r['poi']:<16} {label:<7} {r['n_poses']:>8} {r['n_nodes']:>8} "
              f"{r['n_edges']:>8} {r['vhl_max_freq']:>8.3f} {r['kin_max_freq']:>8.3f}")


if __name__ == '__main__':
    main()
