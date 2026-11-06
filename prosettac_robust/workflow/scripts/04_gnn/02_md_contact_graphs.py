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
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from gnn_contract import periodic_minima,require_verified_ligand_graph
from protac_qc import ValidationError,sha256
from md_mapping import load_mapped_production
from gnn_contract import configure_native_pyg
DISABLED_OPTIONAL_GNN_BACKENDS=configure_native_pyg()
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

VERIFIED_SMILES=json.loads((Path(__file__).resolve().parents[2]/'config/molecules.json').read_text())['compounds']['SJF8240']['smiles']

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
    work=tpr.parent
    tracked=[tpr,xtc,work/'md_preparation.json',work/'prod.receipt.json']
    input_hashes={str(p):sha256(p) for p in tracked}
    u,groups=load_mapped_production(work/'md_preparation.json',tpr,xtc,work/'prod.receipt.json')
    protac=require_verified_ligand_graph(u,VERIFIED_SMILES)
    if set(map(int,protac.indices))!=set(map(int,groups['protac'].indices)):
        raise ValidationError('Ligand residue differs from executed atom mapping')
    vhl_ca=groups['vhl'].select_atoms('name CA');kin_ca=groups['poi'].select_atoms('name CA')
    if not len(vhl_ca) or not len(kin_ca):raise ValidationError('Mapped protein role has no Cα atoms')
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

    joint_counts=np.zeros((n_vhl,n_kin),dtype=np.uint64)
    last_time=None
    for ts in u.trajectory[::stride]:
        if not np.isfinite(ts.time) or (last_time is not None and ts.time<=last_time):raise ValidationError('Invalid or repeated trajectory time')
        last_time=ts.time
        ptc_pos = protac.positions   # (n_ptc, 3)
        vhl_pos = vhl_ca.positions   # (n_vhl, 3)
        kin_pos = kin_ca.positions   # (n_kin, 3)

        # This schema uses ligand-heavy to protein Cα distances, explicitly.
        min_dv=periodic_minima(ptc_pos,vhl_pos,ts.dimensions)
        min_dk=periodic_minima(ptc_pos,kin_pos,ts.dimensions)
        joint_counts+=np.outer(min_dv<CONTACT_THRESH,min_dk<CONTACT_THRESH)
        vhl_contact_count  += (min_dv < CONTACT_THRESH)
        kin_contact_count  += (min_dk < CONTACT_THRESH)
        vhl_dist_sum       += min_dv
        kin_dist_sum       += min_dk
        n_frames += 1

    if any(sha256(p)!=input_hashes[str(p)] for p in tracked):raise ValidationError('MD graph inputs changed during extraction')
    if n_frames == 0:
        raise ValidationError('No usable production frames')

    return {
        'vhl_resids':  vhl_resids,
        'kin_resids':  kin_resids,
        'vhl_freq':    vhl_contact_count / n_frames,
        'kin_freq':    kin_contact_count / n_frames,
        'vhl_mean_d':  vhl_dist_sum / n_frames,
        'kin_mean_d':  kin_dist_sum / n_frames,
        'n_frames':    n_frames,
        'n_ptc_atoms': len(protac),
        'joint_freq':joint_counts/n_frames,
        'input_sha256':input_hashes,
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

        tpr=pose_dir/'solute.tpr'
        xtc=pose_dir/'prod.xtc'
        gro=pose_dir/'system.gro' # Mapping comes from the executed preparation, not residue-number heuristics.

        print(f"    {pose_id} ...", end=' ', flush=True)
        result = compute_pose_contacts(tpr, xtc, gro, kin_min, stride)

        if result is None:
            print("FAILED")
            raise ValidationError('A declared completed MD pose failed graph extraction')

        print(f"n_frames={result['n_frames']}  "
              f"vhl_max_freq={result['vhl_freq'].max():.3f}  "
              f"kin_max_freq={result['kin_freq'].max():.3f}", flush=True)

        pose_results.append(result)
        pose_meta.append({
            'poi': poi, 'pose_id': pose_id,
            'n_frames': result['n_frames'],
            'input_sha256':json.dumps(result['input_sha256'],sort_keys=True),
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
    for result in pose_results[1:]:
        if not np.array_equal(result['vhl_resids'],vhl_resids) or not np.array_equal(result['kin_resids'],kin_resids):
            raise ValidationError('Unaligned residue identities across MD graph poses')
    joint_freq_mean=np.mean([r['joint_freq'] for r in pose_results],axis=0)
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
    # Historical graph used product-of-marginals. This derivative uses observed joint contact.
    # Graph structure: VHL nodes connected to kinase nodes with edge weight
    # = product of their PROTAC contact frequencies (bridging co-occurrence proxy)

    edge_src, edge_dst, edge_freq, edge_dvhl, edge_dkin = [], [], [], [], []

    for i in range(n_vhl):
        for j in range(n_kin):
            fv = vhl_freq_mean[i]
            fk = kin_freq_mean[j]
            co_freq = joint_freq_mean[i,j]  # Observed simultaneous contact, averaged equally over poses
            if co_freq > EDGE_FREQ_THRESH:
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
        'edge_attr':np.column_stack([edge_freq,edge_dvhl,edge_dkin]).astype(np.float32),
    }

    if HAS_TORCH:
        ea=graph_dict['edge_attr']
        graph_dict['pyg'] = Data(
            x          = torch.tensor(node_feats, dtype=torch.float32),
            edge_index = torch.tensor([edge_src, edge_dst], dtype=torch.long),
            edge_attr  = torch.tensor(ea, dtype=torch.float32),
            y          = torch.tensor([int(cfg['degrader']) if cfg['degrader'] is not None else -1], dtype=torch.long),
            poi        = cfg.get('graph_name', poi),
            feature_schema='SJF8240-periodic-contacts-v1',
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
    parser.add_argument('--root',type=Path,required=True,help='Explicit regenerated screen directory')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--input-provenance',type=Path,required=True,help='Reviewed chemistry/manifest receipt for regenerated MD')
    args = parser.parse_args()
    if args.merge:raise ValidationError('Implicit merge is disabled; regenerate an explicit versioned graph packet')
    if args.stride<1:raise ValidationError('Positive frame stride required')
    for original in ['/home/ayushregmibagale/protac-md','/mnt/d/PROTAC degradation Prediction']:
        if args.output.resolve().is_relative_to(Path(original)):raise ValidationError('Output points into original work')
        if args.root.resolve().is_relative_to(Path(original)):raise ValidationError('Feature production requires separately rebuilt systems; use read-only audit tools for originals')
    if args.output.exists():raise ValidationError('Graph output already exists')
    global RANDOM_SCREEN,MANIFEST,GNN_DIR,OUT_PT,OUT_NPZ,OUT_META,OUT_POSE
    RANDOM_SCREEN=args.root.resolve();MANIFEST=RANDOM_SCREEN/'manifest.tsv';GNN_DIR=args.output.resolve()
    OUT_PT=GNN_DIR/'md_contact_graphs.pt';OUT_NPZ=GNN_DIR/'md_contact_graphs_np.npz'
    OUT_META=GNN_DIR/'md_graph_metadata.tsv';OUT_POSE=GNN_DIR/'md_contact_per_pose.tsv'
    provenance=json.loads(args.input_provenance.read_text())
    if provenance.get('compound')!='SJF8240' or provenance.get('ligand_inchikey')!='SCUZPTNDZHEHLB-ZWUDBOSPSA-N' or provenance.get('manifest_sha256')!=sha256(MANIFEST):
        raise ValidationError('Reviewed verified-SJF8240 screen provenance required')
    source_paths=[Path(__file__),Path(__file__).resolve().parents[1]/'lib/gnn_contract.py',Path(__file__).resolve().parents[1]/'lib/md_mapping.py',args.input_provenance,MANIFEST]
    source_hashes={str(p):sha256(p) for p in source_paths}
    labels=provenance.get('experimental_labels',{})
    if not labels or any(type(v) is not int or v not in (-1,0,1) for v in labels.values()):
        raise ValidationError('Explicit assay-reviewed labels or -1 unknown required; no hardcoded labels')
    POI_CONFIG.clear()
    POI_CONFIG.update({p:{'kin_min':0,'degrader':None if y==-1 else bool(y)} for p,y in labels.items()})

    if not MANIFEST.exists():
        print(f"[ERROR] manifest.tsv not found at {MANIFEST}")
        sys.exit(1)

    manifest = pd.read_csv(MANIFEST, sep='\t')
    done=manifest[manifest['stage'].isin(['PRODUCTION_COMPLETE','PROD_COMPLETE','ANALYSIS_DONE'])]

    pois=args.poi or sorted(POI_CONFIG)
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
            raise ValidationError('Target is outside explicit reviewed cohort: '+poi)

        poi_done = done[done['poi'] == poi].to_dict('records')
        print(f"\n[{poi}] {len(poi_done)} ANALYSIS_DONE poses:", flush=True)

        if not poi_done:
            raise ValidationError('Declared target has no completed production poses: '+poi)

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
        raise ValidationError('No valid graphs built')

    if any(sha256(p)!=source_hashes[str(p)] for p in source_paths):raise ValidationError('Graph extraction source/provenance/manifest changed during execution')
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
                  'edge_src', 'edge_dst', 'edge_freq','edge_attr'):
            if k in g:
                np_dict[f"{p}/{k}"] = g[k]
        np_dict[f"{p}/degrader"] = np.array([g['degrader']])
    np.savez_compressed(OUT_NPZ, **np_dict)
    print(f"NumPy graphs saved: {OUT_NPZ}")
    receipt={'status':'GRAPH_PACKET_REQUIRES_REVIEW','feature_schema':'SJF8240-periodic-contacts-v1',
        'distance_definition':'Ligand heavy atoms to protein Cα, triclinic minimum image, Angstrom',
        'joint_edges':'Observed same-frame co-contact, then equal weighting over poses; no independence assumption',
        'edge_frequency_threshold':EDGE_FREQ_THRESH,'stride':args.stride,
        'source_and_plan_sha256':source_hashes,'source_sha256':sha256(Path(__file__)),'manifest_sha256':sha256(MANIFEST),
        'input_provenance_sha256':sha256(args.input_provenance),
        'output_sha256':{p.name:sha256(p) for p in GNN_DIR.iterdir() if p.is_file()},
        'automatic_training_release':False}
    (GNN_DIR/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

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
