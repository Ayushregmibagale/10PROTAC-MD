#!/usr/bin/env python3
"""
Tier2b: Build pRosettaC ensemble contact frequency graphs.

For each POI, aggregates all ~20,000 combined PDBs into a single graph where:
  - Nodes  = VHL residues (Cα) + kinase residues (Cα) + PROTAC heavy atoms
  - Edges  = pairs in contact in ≥ EDGE_FREQ_THRESH fraction of poses
  - Edge features = [contact_frequency, mean_distance, min_distance]
  - Node features = [is_vhl, is_kin, is_ptc, seq_pos_norm, elem_C, elem_N, elem_O, elem_F, elem_S]
  - Graph label   = 1 (degrader) / 0 (non-degrader)

Key insight: this is the static structural analogue of an MD contact-frequency graph.
For CDK2/4/6/DAPK1/IRAK1/TLK2 (zero bridging), edges near the bridging region will
all be near-zero — the graph encodes the geometric failure mode.

Outputs:
  analysis/gnn/ensemble_graphs.pt  — list of PyTorch Geometric Data objects
  analysis/gnn/graph_metadata.tsv  — node counts, edge counts, label per POI

Usage:
  python3 scripts/02_analysis/04a_build_ensemble_graphs.py
  python3 scripts/02_analysis/04a_build_ensemble_graphs.py --poi DDR2 CDK2 --workers 12
"""
import argparse
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from multiprocessing import Pool

try:
    import torch
    from torch_geometric.data import Data
    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False
    print("[WARNING] torch / torch_geometric not available — will save as numpy dicts (.npz)")

ANALYSIS_DIR = Path.home() / 'protac-md' / 'analysis'
ROSETTA_BASE = Path.home() / 'PROTAC degradation Prediction' / 'pRosettaC' / 'runs'
GNN_DIR      = ANALYSIS_DIR / 'gnn'
OUT_PT       = GNN_DIR / 'ensemble_graphs.pt'
OUT_NPZ      = GNN_DIR / 'ensemble_graphs_np.npz'
OUT_META     = GNN_DIR / 'graph_metadata.tsv'

# ── Contact frequency parameters ─────────────────────────────────────────────
CONTACT_THRESH  = 8.0   # Å — include in contact count if within this distance
EDGE_FREQ_THRESH = 0.005 # include edge if in contact in ≥ 0.5% of poses (~100/20k)
VHL_MAX         = 103   # VHL residues: chain A, 1-103

# ── POI configuration ─────────────────────────────────────────────────────────
POIS_CONFIG = {
    'DDR2':        {'base': 'high_affinity_degraders',      'kin_chain': 'A', 'kin_min': 553,  'degrader': True,  'artifact': False},
    'RIPK2_full':  {'base': 'high_affinity_degraders',      'kin_chain': 'A', 'kin_min': 200,  'degrader': True,  'artifact': False},
    'MAPK14_full': {'base': 'high_affinity_degraders',      'kin_chain': 'A', 'kin_min': 200,  'degrader': True,  'artifact': False},
    'MET':         {'base': 'high_affinity_degraders',      'kin_chain': 'A', 'kin_min': 1053, 'degrader': True,  'artifact': False},
    'EPHB2':       {'base': 'high_affinity_degraders',      'kin_chain': 'A', 'kin_min': 615,  'degrader': True,  'artifact': False},
    'ABL1':        {'base': 'high_affinity_no_degradation', 'kin_chain': 'A', 'kin_min': 235,  'degrader': False, 'artifact': False},
    'AXL':         {'base': 'high_affinity_no_degradation', 'kin_chain': 'A', 'kin_min': 521,  'degrader': False, 'artifact': False},
    'EPHA2':       {'base': 'high_affinity_no_degradation', 'kin_chain': 'A', 'kin_min': 600,  'degrader': False, 'artifact': False},
    # Chain-D reruns (Aug 12): 0% bridging confirmed, geometry-blocked non-degraders
    'MAP4K5':      {'base': 'high_affinity_no_degradation', 'kin_chain': 'D', 'kin_min': 0,    'degrader': False, 'artifact': False, 'src_subdir': 'Results'},
    'SLK':         {'base': 'high_affinity_no_degradation', 'kin_chain': 'D', 'kin_min': 0,    'degrader': False, 'artifact': False, 'src_subdir': 'Results'},
    # CDK17: confirmed NON-DEGRADER Aug 12 (was incorrectly labelled degrader)
    'CDK17':       {'base': 'low_affinity_degraders',       'kin_chain': 'D', 'kin_min': 0,    'degrader': False, 'artifact': False},
    'MAPK14':      {'base': 'low_affinity_degraders',       'kin_chain': 'A', 'kin_min': 104,  'degrader': True,  'artifact': True},
    'TNIK':        {'base': 'low_affinity_degraders',       'kin_chain': 'D', 'kin_min': 0,    'degrader': True,  'artifact': True,  'src_subdir': 'Results'},
    'CDK2':        {'base': 'low_affinity_no_degradation',  'kin_chain': 'D', 'kin_min': 0,    'degrader': False, 'artifact': False},
    'CDK4':        {'base': 'low_affinity_no_degradation',  'kin_chain': 'D', 'kin_min': 0,    'degrader': False, 'artifact': False},
    'CDK6':        {'base': 'low_affinity_no_degradation',  'kin_chain': 'D', 'kin_min': 0,    'degrader': False, 'artifact': False},
    'DAPK1':       {'base': 'low_affinity_no_degradation',  'kin_chain': 'D', 'kin_min': 0,    'degrader': False, 'artifact': False},
    'IRAK1':       {'base': 'low_affinity_no_degradation',  'kin_chain': 'D', 'kin_min': 0,    'degrader': False, 'artifact': False},
    'TLK2':        {'base': 'low_affinity_no_degradation',  'kin_chain': 'D', 'kin_min': 0,    'degrader': False, 'artifact': False},
    # Kinome screen Phase 2 — prediction targets (label unknown → degrader=None)
    # Renumbered to res 200+ (Group 1); data on D drive after rsync from Ernie
    'SIK3_full':   {'kin_chain': 'A', 'kin_min': 200, 'degrader': None, 'artifact': False,
                    'path': Path('/mnt/d/kinome_results/sik3_full/Patchdock_Results')},
    'RPS6KB2_full':{'kin_chain': 'A', 'kin_min': 200, 'degrader': None, 'artifact': False,
                    'path': Path('/mnt/d/kinome_results/rps6kb2_full/Patchdock_Results')},
}

# One-hot element index for PROTAC heavy atoms
ELEM_ORDER = ['C', 'N', 'O', 'F', 'S']
ELEM_IDX   = {e: i for i, e in enumerate(ELEM_ORDER)}


# ── Node discovery ────────────────────────────────────────────────────────────

def discover_nodes(pdb_path: str, kin_chain: str, kin_min: int) -> dict:
    """Parse one PDB to discover canonical node sets for this POI."""
    vhl_res    = {}   # resnum → residue_name
    kin_res    = {}   # resnum → residue_name (canonical order)
    ptc_atoms  = {}   # atom_idx → (atom_name, element)

    with open(pdb_path) as f:
        for line in f:
            rec = line[:6].strip()
            if rec not in ('ATOM', 'HETATM'):
                continue
            aname  = line[12:16].strip()
            resname = line[17:20].strip()
            chain  = line[21]
            resnum = int(line[22:26])
            elem   = line[76:78].strip() if len(line) > 77 else aname[0]
            serial = int(line[6:11])

            # VHL residues
            if chain == 'A' and 1 <= resnum <= VHL_MAX and aname == 'CA':
                vhl_res[resnum] = resname

            # Kinase residues
            if kin_chain == 'D' and chain == 'D' and aname == 'CA':
                kin_res[resnum] = resname
            elif kin_chain == 'A' and chain == 'A' and resnum >= kin_min and aname == 'CA':
                kin_res[resnum] = resname

            # PROTAC heavy atoms (chain X, no H)
            if chain == 'X' and not aname.startswith('H'):
                ptc_atoms[serial - 1] = (aname, elem if elem in ELEM_IDX else aname[0])

    return {
        'vhl_res': dict(sorted(vhl_res.items())),
        'kin_res': dict(sorted(kin_res.items())),
        'ptc_atoms': dict(sorted(ptc_atoms.items())),
    }


def build_node_features(nodes: dict) -> np.ndarray:
    """
    Build node feature matrix from discovered nodes.
    Features: [is_vhl, is_kin, is_ptc, seq_pos_norm, C, N, O, F, S]  — 9 dims
    """
    rows = []
    n_vhl = len(nodes['vhl_res'])
    n_kin = len(nodes['kin_res'])
    n_ptc = len(nodes['ptc_atoms'])

    # VHL residues
    for i, (resnum, resname) in enumerate(nodes['vhl_res'].items()):
        feat = np.zeros(9, dtype=np.float32)
        feat[0] = 1.0                            # is_vhl
        feat[3] = resnum / VHL_MAX               # seq_pos_norm
        feat[4] = 1.0                            # Cα = C
        rows.append(feat)

    # Kinase residues
    kin_resnums = list(nodes['kin_res'].keys())
    rmin, rmax = (min(kin_resnums), max(kin_resnums)) if kin_resnums else (1, 1)
    for i, (resnum, resname) in enumerate(nodes['kin_res'].items()):
        feat = np.zeros(9, dtype=np.float32)
        feat[1] = 1.0                                              # is_kin
        feat[3] = (resnum - rmin) / max(rmax - rmin, 1)           # seq_pos_norm
        feat[4] = 1.0                                              # Cα = C
        rows.append(feat)

    # PROTAC heavy atoms
    n_ptc_total = max(nodes['ptc_atoms'].keys()) + 1 if nodes['ptc_atoms'] else 1
    for i, (atom_idx, (aname, elem)) in enumerate(nodes['ptc_atoms'].items()):
        feat = np.zeros(9, dtype=np.float32)
        feat[2] = 1.0                                       # is_ptc
        feat[3] = atom_idx / n_ptc_total                    # seq_pos_norm
        eidx = ELEM_IDX.get(elem, ELEM_IDX.get(elem[0], 0))
        feat[4 + eidx] = 1.0                                # element one-hot
        rows.append(feat)

    return np.array(rows, dtype=np.float32)


# ── Contact accumulation (parallel workers) ───────────────────────────────────

def _parse_pdb_contacts(args):
    """
    Worker: parse one PDB → return (vhl_ptc_dists, kin_ptc_dists, vhl_kin_dists).
    Each is a 2D float32 array of min-atom distances within CONTACT_THRESH.
    Returns None if parsing fails.
    """
    pdb_path, kin_chain, kin_min, n_vhl, n_kin, n_ptc, vhl_resnums, kin_resnums, ptc_order = args

    vhl_coords  = np.full((n_vhl, 3), np.nan, dtype=np.float32)
    kin_coords  = np.full((n_kin, 3), np.nan, dtype=np.float32)
    ptc_coords  = np.full((n_ptc, 3), np.nan, dtype=np.float32)

    vhl_idx = {r: i for i, r in enumerate(vhl_resnums)}
    kin_idx = {r: i for i, r in enumerate(kin_resnums)}
    ptc_map = {s: i for i, s in enumerate(ptc_order)}

    try:
        with open(pdb_path) as f:
            for line in f:
                rec = line[:6].strip()
                if rec not in ('ATOM', 'HETATM'):
                    continue
                aname  = line[12:16].strip()
                chain  = line[21]
                resnum = int(line[22:26])
                serial = int(line[6:11])
                coord  = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])],
                                  dtype=np.float32)

                # VHL Cα
                if chain == 'A' and resnum in vhl_idx and aname == 'CA':
                    vhl_coords[vhl_idx[resnum]] = coord

                # Kinase Cα
                if kin_chain == 'D' and chain == 'D' and resnum in kin_idx and aname == 'CA':
                    kin_coords[kin_idx[resnum]] = coord
                elif kin_chain == 'A' and chain == 'A' and resnum >= kin_min and resnum in kin_idx and aname == 'CA':
                    kin_coords[kin_idx[resnum]] = coord

                # PROTAC heavy atoms (by serial order)
                if chain == 'X' and not aname.startswith('H'):
                    pidx = serial - 1  # 0-indexed
                    if pidx in ptc_map:
                        ptc_coords[ptc_map[pidx]] = coord

    except Exception:
        return None

    # Compute pairwise distances; return contact indicator matrices
    # PROTAC-VHL: (n_ptc, n_vhl)
    ptc_vhl = np.full((n_ptc, n_vhl), np.inf, dtype=np.float32)
    ptc_kin = np.full((n_ptc, n_kin), np.inf, dtype=np.float32)
    vhl_kin = np.full((n_vhl, n_kin), np.inf, dtype=np.float32)

    valid_ptc = ~np.any(np.isnan(ptc_coords), axis=1)
    valid_vhl = ~np.any(np.isnan(vhl_coords), axis=1)
    valid_kin = ~np.any(np.isnan(kin_coords), axis=1)

    if valid_ptc.any() and valid_vhl.any():
        p = ptc_coords[valid_ptc]
        v = vhl_coords[valid_vhl]
        diffs = p[:, None, :] - v[None, :, :]          # (n_ptc_ok, n_vhl_ok, 3)
        d = np.sqrt((diffs ** 2).sum(axis=2))           # (n_ptc_ok, n_vhl_ok)
        pidx_ok = np.where(valid_ptc)[0]
        vidx_ok = np.where(valid_vhl)[0]
        ptc_vhl[np.ix_(pidx_ok, vidx_ok)] = d

    if valid_ptc.any() and valid_kin.any():
        p = ptc_coords[valid_ptc]
        k = kin_coords[valid_kin]
        diffs = p[:, None, :] - k[None, :, :]
        d = np.sqrt((diffs ** 2).sum(axis=2))
        pidx_ok = np.where(valid_ptc)[0]
        kidx_ok = np.where(valid_kin)[0]
        ptc_kin[np.ix_(pidx_ok, kidx_ok)] = d

    if valid_vhl.any() and valid_kin.any():
        v = vhl_coords[valid_vhl]
        k = kin_coords[valid_kin]
        diffs = v[:, None, :] - k[None, :, :]
        d = np.sqrt((diffs ** 2).sum(axis=2))
        vidx_ok = np.where(valid_vhl)[0]
        kidx_ok = np.where(valid_kin)[0]
        vhl_kin[np.ix_(vidx_ok, kidx_ok)] = d

    # Only return pairs within threshold (for memory efficiency: return binary masks + distances)
    ptc_vhl[ptc_vhl > CONTACT_THRESH] = 0.0
    ptc_kin[ptc_kin > CONTACT_THRESH] = 0.0
    vhl_kin[vhl_kin > CONTACT_THRESH] = 0.0

    return (ptc_vhl, ptc_kin, vhl_kin)


# ── Graph builder ─────────────────────────────────────────────────────────────

def build_poi_graph(poi: str, n_workers: int = 12) -> dict | None:
    """Build ensemble contact frequency graph for one POI."""
    cfg = POIS_CONFIG[poi]
    src_subdir = cfg.get('src_subdir', 'Patchdock_Results')

    # Kinome POIs use a direct path override instead of ROSETTA_BASE hierarchy
    if 'path' in cfg:
        pd_dir = cfg['path']
        if not pd_dir.exists():
            print(f"  [{poi}] path not found: {pd_dir} — skipping")
            return None
        pdb_files = sorted(pd_dir.glob('combined_*.pdb'))
    elif src_subdir == 'Results':
        poi_dir = ROSETTA_BASE / cfg['base'] / poi
        results_dir = poi_dir / 'Results'
        if not results_dir.exists():
            print(f"  [{poi}] Results/ not found — skipping")
            return None
        pdb_files = sorted(results_dir.glob('cluster*/combined_*.pdb'))
    else:
        poi_dir = ROSETTA_BASE / cfg['base'] / poi
        pd_dir = poi_dir / 'Patchdock_Results'
        if not pd_dir.exists():
            print(f"  [{poi}] Patchdock_Results not found — skipping")
            return None
        pdb_files = sorted(pd_dir.glob('combined_*.pdb'))

    if not pdb_files:
        print(f"  [{poi}] no combined_*.pdb files in {src_subdir} — skipping")
        return None

    print(f"  [{poi}] {len(pdb_files)} PDBs — discovering node set...", flush=True)

    # ── Step 1: discover canonical nodes from first PDB ──
    nodes = discover_nodes(str(pdb_files[0]), cfg['kin_chain'], cfg['kin_min'])
    n_vhl = len(nodes['vhl_res'])
    n_kin = len(nodes['kin_res'])
    n_ptc = len(nodes['ptc_atoms'])
    n_total = n_vhl + n_kin + n_ptc

    vhl_resnums = list(nodes['vhl_res'].keys())
    kin_resnums = list(nodes['kin_res'].keys())
    ptc_order   = list(nodes['ptc_atoms'].keys())   # serial-1 indices

    print(f"  [{poi}] nodes: VHL={n_vhl}  KIN={n_kin}  PTC={n_ptc}  total={n_total}", flush=True)

    # ── Step 2: accumulate contacts across all PDBs ──
    sum_ptc_vhl  = np.zeros((n_ptc, n_vhl), dtype=np.float64)
    sum_ptc_kin  = np.zeros((n_ptc, n_kin), dtype=np.float64)
    sum_vhl_kin  = np.zeros((n_vhl, n_kin), dtype=np.float64)
    # For mean distance: need sum of distances where in contact
    dsum_ptc_vhl = np.zeros((n_ptc, n_vhl), dtype=np.float64)
    dsum_ptc_kin = np.zeros((n_ptc, n_kin), dtype=np.float64)
    dsum_vhl_kin = np.zeros((n_vhl, n_kin), dtype=np.float64)
    dmin_ptc_vhl = np.full((n_ptc, n_vhl), np.inf, dtype=np.float32)
    dmin_ptc_kin = np.full((n_ptc, n_kin), np.inf, dtype=np.float32)
    dmin_vhl_kin = np.full((n_vhl, n_kin), np.inf, dtype=np.float32)

    args_iter = [
        (str(p), cfg['kin_chain'], cfg['kin_min'],
         n_vhl, n_kin, n_ptc, vhl_resnums, kin_resnums, ptc_order)
        for p in pdb_files
    ]

    n_parsed = 0
    batch_size = 2000
    for batch_start in range(0, len(args_iter), batch_size):
        batch = args_iter[batch_start:batch_start + batch_size]
        with Pool(n_workers) as pool:
            for result in pool.imap_unordered(_parse_pdb_contacts, batch, chunksize=50):
                if result is None:
                    continue
                pv, pk, vk = result
                n_parsed += 1

                # Accumulate: pv/pk/vk have distance where in contact, 0 where not
                contact_pv = pv > 0
                contact_pk = pk > 0
                contact_vk = vk > 0

                sum_ptc_vhl[contact_pv] += 1
                sum_ptc_kin[contact_pk] += 1
                sum_vhl_kin[contact_vk] += 1

                dsum_ptc_vhl[contact_pv] += pv[contact_pv]
                dsum_ptc_kin[contact_pk] += pk[contact_pk]
                dsum_vhl_kin[contact_vk] += vk[contact_vk]

                dmin_ptc_vhl = np.where((pv > 0) & (pv < dmin_ptc_vhl), pv.astype(np.float32), dmin_ptc_vhl)
                dmin_ptc_kin = np.where((pk > 0) & (pk < dmin_ptc_kin), pk.astype(np.float32), dmin_ptc_kin)
                dmin_vhl_kin = np.where((vk > 0) & (vk < dmin_vhl_kin), vk.astype(np.float32), dmin_vhl_kin)

        print(f"  [{poi}] processed {min(batch_start + batch_size, len(pdb_files))}/{len(pdb_files)} "
              f"({n_parsed} valid)", flush=True)

    # ── Step 3: compute frequencies ──
    freq_ptc_vhl  = sum_ptc_vhl / n_parsed
    freq_ptc_kin  = sum_ptc_kin / n_parsed
    freq_vhl_kin  = sum_vhl_kin / n_parsed
    mean_ptc_vhl  = np.where(sum_ptc_vhl > 0, dsum_ptc_vhl / sum_ptc_vhl, CONTACT_THRESH)
    mean_ptc_kin  = np.where(sum_ptc_kin > 0, dsum_ptc_kin / sum_ptc_kin, CONTACT_THRESH)
    mean_vhl_kin  = np.where(sum_vhl_kin > 0, dsum_vhl_kin / sum_vhl_kin, CONTACT_THRESH)
    dmin_ptc_vhl[dmin_ptc_vhl == np.inf] = CONTACT_THRESH
    dmin_ptc_kin[dmin_ptc_kin == np.inf] = CONTACT_THRESH
    dmin_vhl_kin[dmin_vhl_kin == np.inf] = CONTACT_THRESH

    # ── Step 4: build edge index and features ──
    # Node index layout: [VHL_0..VHL_n_vhl-1 | KIN_0..KIN_n_kin-1 | PTC_0..PTC_n_ptc-1]
    vhl_off = 0
    kin_off = n_vhl
    ptc_off = n_vhl + n_kin

    edge_src, edge_dst, edge_freq, edge_mean_d, edge_min_d = [], [], [], [], []

    def add_edges(freq_mat, mean_mat, min_mat, src_off, dst_off):
        rows, cols = np.where(freq_mat >= EDGE_FREQ_THRESH)
        for r, c in zip(rows, cols):
            i, j = int(r) + src_off, int(c) + dst_off
            f = float(freq_mat[r, c])
            md = float(mean_mat[r, c])
            mn = float(min_mat[r, c])
            # Undirected: add both directions
            edge_src.extend([i, j])
            edge_dst.extend([j, i])
            edge_freq.extend([f, f])
            edge_mean_d.extend([md, md])
            edge_min_d.extend([mn, mn])

    add_edges(freq_ptc_vhl, mean_ptc_vhl, dmin_ptc_vhl, ptc_off, vhl_off)
    add_edges(freq_ptc_kin, mean_ptc_kin, dmin_ptc_kin, ptc_off, kin_off)
    add_edges(freq_vhl_kin, mean_vhl_kin, dmin_vhl_kin, vhl_off, kin_off)

    n_edges = len(edge_src) // 2

    print(f"  [{poi}] {n_parsed} poses processed | "
          f"n_edges={n_edges}  "
          f"ptc_vhl_edges={(freq_ptc_vhl >= EDGE_FREQ_THRESH).sum()}  "
          f"ptc_kin_edges={(freq_ptc_kin >= EDGE_FREQ_THRESH).sum()}  "
          f"max_ptc_vhl_freq={freq_ptc_vhl.max():.3f}  "
          f"max_ptc_kin_freq={freq_ptc_kin.max():.3f}", flush=True)

    # ── Step 5: node features ──
    node_feats = build_node_features(nodes)  # (n_total, 9)

    # ── Step 6: package as dict (PyG Data if available) ──
    graph_dict = {
        'poi':           poi,
        'degrader':      int(cfg['degrader']) if cfg['degrader'] is not None else -1,
        'artifact':      int(cfg['artifact']),
        'n_vhl':         n_vhl,
        'n_kin':         n_kin,
        'n_ptc':         n_ptc,
        'n_poses':       n_parsed,
        'n_edges':       n_edges,
        'node_feats':    node_feats.astype(np.float32),       # (n_nodes, 9)
        'edge_src':      np.array(edge_src,    dtype=np.int64),
        'edge_dst':      np.array(edge_dst,    dtype=np.int64),
        'edge_freq':     np.array(edge_freq,   dtype=np.float32),
        'edge_mean_d':   np.array(edge_mean_d, dtype=np.float32),
        'edge_min_d':    np.array(edge_min_d,  dtype=np.float32),
        # Raw frequency matrices (for analysis)
        'freq_ptc_vhl':  freq_ptc_vhl.astype(np.float32),
        'freq_ptc_kin':  freq_ptc_kin.astype(np.float32),
        'freq_vhl_kin':  freq_vhl_kin.astype(np.float32),
        'vhl_resnums':   np.array(vhl_resnums, dtype=np.int32),
        'kin_resnums':   np.array(kin_resnums, dtype=np.int32),
    }

    if HAS_TORCH:
        ei = torch.tensor([edge_src, edge_dst], dtype=torch.long)
        ea = torch.tensor(
            np.column_stack([edge_freq, edge_mean_d, edge_min_d]), dtype=torch.float32
        )
        graph_dict['pyg'] = Data(
            x           = torch.tensor(node_feats, dtype=torch.float32),
            edge_index  = ei,
            edge_attr   = ea,
            y           = torch.tensor([int(cfg['degrader']) if cfg['degrader'] is not None else -1], dtype=torch.long),
            poi         = poi,
            artifact    = int(cfg['artifact']),
            n_vhl       = n_vhl,
            n_kin       = n_kin,
            n_ptc       = n_ptc,
        )

    return graph_dict


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--poi',     nargs='+', default=None)
    parser.add_argument('--workers', type=int,  default=12)
    parser.add_argument('--merge',   action='store_true',
                        help='Load existing .pt and merge new POIs (skip POIs already present)')
    args = parser.parse_args()

    GNN_DIR.mkdir(parents=True, exist_ok=True)

    pois = args.poi or list(POIS_CONFIG.keys())

    print(f"\nTier2b: pRosettaC Ensemble Graph Builder")
    print(f"Contact threshold: {CONTACT_THRESH} Å | Min edge freq: {EDGE_FREQ_THRESH}")
    print(f"Workers: {args.workers} | POIs: {pois}\n")

    # Load existing graphs if --merge requested
    existing_pyg   = []
    existing_pois  = set()
    if args.merge and HAS_TORCH and OUT_PT.exists():
        existing_pyg  = torch.load(OUT_PT, weights_only=False)
        existing_pois = {g.poi for g in existing_pyg}
        print(f"Loaded {len(existing_pyg)} existing graphs: {sorted(existing_pois)}")
        pois = [p for p in pois if p not in existing_pois]
        print(f"Building {len(pois)} new POIs: {pois}\n")

    all_graphs  = []
    meta_rows   = []
    pyg_graphs  = list(existing_pyg)   # start with existing

    for poi in pois:
        g = build_poi_graph(poi, args.workers)
        if g is None:
            continue
        all_graphs.append(g)
        if 'pyg' in g:
            pyg_graphs.append(g['pyg'])
        meta_rows.append({
            'poi':      g['poi'],
            'degrader': g['degrader'],
            'artifact': g['artifact'],
            'n_vhl':    g['n_vhl'],
            'n_kin':    g['n_kin'],
            'n_ptc':    g['n_ptc'],
            'n_nodes':  g['n_vhl'] + g['n_kin'] + g['n_ptc'],
            'n_edges':  g['n_edges'],
            'n_poses':  g['n_poses'],
        })

    if not all_graphs and not existing_pyg:
        print("No graphs built.")
        return

    # Add existing POIs to meta_rows (from graph attrs) if merging
    if args.merge and existing_pyg and OUT_META.exists():
        ex_meta = pd.read_csv(OUT_META, sep='\t')
        ex_meta = ex_meta[ex_meta['poi'].isin(existing_pois)]
        meta_rows = ex_meta.to_dict('records') + meta_rows

    # Save metadata
    meta_df = pd.DataFrame(meta_rows)
    meta_df.to_csv(OUT_META, sep='\t', index=False)
    print(f"\nMetadata saved: {OUT_META}")

    # Save PyG graphs
    if HAS_TORCH and pyg_graphs:
        torch.save(pyg_graphs, OUT_PT)
        print(f"PyG graphs saved: {OUT_PT}  ({len(pyg_graphs)} graphs)")

    # Save numpy fallback
    np_dict = {}
    for g in all_graphs:
        p = g['poi']
        for k in ('node_feats', 'edge_src', 'edge_dst', 'edge_freq', 'edge_mean_d', 'edge_min_d',
                  'freq_ptc_vhl', 'freq_ptc_kin', 'freq_vhl_kin', 'vhl_resnums', 'kin_resnums'):
            if k in g:
                np_dict[f"{p}/{k}"] = g[k]
        np_dict[f"{p}/degrader"] = np.array([g['degrader']])
        np_dict[f"{p}/artifact"] = np.array([g['artifact']])
    np.savez_compressed(OUT_NPZ, **np_dict)
    print(f"NumPy graphs saved: {OUT_NPZ}")

    # Summary
    print(f"\n{'='*70}")
    print(f"ENSEMBLE GRAPH BUILD COMPLETE")
    print(f"{'='*70}")
    print(f"\n{'POI':<16} {'label':<7} {'art':>4} {'n_nodes':>8} {'n_edges':>8} {'n_poses':>8}")
    print('-' * 55)
    for r in meta_rows:
        label = 'DEG' if r['degrader'] else 'NON'
        art = 'Y' if r['artifact'] else 'N'
        print(f"{r['poi']:<16} {label:<7} {art:>4} {r['n_nodes']:>8} {r['n_edges']:>8} {r['n_poses']:>8}")


if __name__ == '__main__':
    main()
