#!/usr/bin/env python3
"""
Static feature extraction for MET–VHL–SJF8240 ternary complex clusters.

Chain map (from pRosettaC combined PDB):
  Chain A res   1-103  = VHL (4W9H chain A)
  Chain A res 1053-1345 = MET kinase domain
  Chain B               = ElonginB
  Chain C               = ElonginC
  Chain X               = SJF8240 (PTC residue)
"""

import os, sys, glob, math, warnings, csv, json
import numpy as np
from pathlib import Path

warnings.filterwarnings('ignore')

RESULTS = Path('/home/s1638308/protac/pRosettaC/runs/high_affinity_degraders/MET/Results')
SCORE_SC = Path('/home/s1638308/protac/pRosettaC/runs/high_affinity_degraders/MET/Patchdock_Results/score.sc')
SCORE_SHIFT = 9576.0   # added during score-shifting step; subtract to recover original

# ── Residue ranges ────────────────────────────────────────────────────────────
VHL_CHAIN, VHL_MIN, VHL_MAX   = 'A',    1,  103
MET_CHAIN, MET_MIN, MET_MAX   = 'A', 1053, 1345
ELONG_B_CHAIN                 = 'B'
PROTAC_CHAIN                  = 'X'

# VHL E2-facing proxy: C-terminal region of VHL (connects to ElonginB scaffold;
# CRL2 positions E2 on the ElonginB side).  Use ElonginB Cα centroid as anchor.
# Also: literature E2–substrate Lys optimal distance ~45-60 Å (CRL2).

# ── Score lookup ─────────────────────────────────────────────────────────────
def load_scores():
    scores = {}
    if not SCORE_SC.exists():
        return scores
    with open(SCORE_SC) as f:
        for line in f:
            if line.startswith('SCORE:') and 'total_score' not in line:
                parts = line.split()
                # desc like 'combined_48_10_0001'
                desc = parts[-1]
                raw = float(parts[1])
                scores[desc] = raw + SCORE_SHIFT   # restore original
    return scores

# ── PDB parsing ──────────────────────────────────────────────────────────────
def parse_pdb(path):
    """Return list of atom dicts from ATOM/HETATM records."""
    atoms = []
    with open(path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            try:
                rec    = line[:6].strip()
                serial = int(line[6:11])
                name   = line[12:16].strip()
                resn   = line[17:20].strip()
                chain  = line[21]
                resi   = int(line[22:26])
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except (ValueError, IndexError):
                continue
            atoms.append(dict(rec=rec, serial=serial, name=name, resn=resn,
                               chain=chain, resi=resi, xyz=np.array([x, y, z])))
    return atoms


def select(atoms, chain=None, resi_range=None, resn=None, name=None, exclude_chain=None):
    out = atoms
    if chain:
        out = [a for a in out if a['chain'] == chain]
    if resi_range:
        lo, hi = resi_range
        out = [a for a in out if lo <= a['resi'] <= hi]
    if resn:
        if isinstance(resn, str):
            out = [a for a in out if a['resn'] == resn]
        else:
            out = [a for a in out if a['resn'] in resn]
    if name:
        if isinstance(name, str):
            out = [a for a in out if a['name'] == name]
        else:
            out = [a for a in out if a['name'] in name]
    if exclude_chain:
        out = [a for a in out if a['chain'] != exclude_chain]
    return out


def coords(atoms):
    return np.array([a['xyz'] for a in atoms])


# ── Feature functions ─────────────────────────────────────────────────────────

def interface_contacts(atoms, cutoff=5.0):
    """
    Count heavy-atom pairs < cutoff Å between MET and VHL (chain A res 1-103 + B + C).
    """
    met = select(atoms, chain=MET_CHAIN, resi_range=(MET_MIN, MET_MAX))
    vhl_a = select(atoms, chain=VHL_CHAIN, resi_range=(VHL_MIN, VHL_MAX))
    vhl_bc = [a for a in atoms if a['chain'] in (ELONG_B_CHAIN, 'C')]
    vhl = vhl_a + vhl_bc

    met_co = coords(met)
    vhl_co = coords(vhl)
    if len(met_co) == 0 or len(vhl_co) == 0:
        return 0, 0.0

    # pairwise distances (vectorized)
    diff = met_co[:, None, :] - vhl_co[None, :, :]  # (Nmet, Nvhl, 3)
    dists = np.sqrt((diff**2).sum(axis=2))           # (Nmet, Nvhl)
    n_contacts = int((dists < cutoff).sum())

    # Approximate interface area: each contact pair contributes ~2 Å² (rough)
    # A better estimate: count unique residues on each side with at least 1 contact
    met_mask  = (dists < cutoff).any(axis=1)
    vhl_mask  = (dists < cutoff).any(axis=0)
    n_met_iface  = int(met_mask.sum())
    n_vhl_iface  = int(vhl_mask.sum())
    # Rough surface area: ~15 Å² per interface heavy atom
    iface_area_est = (n_met_iface + n_vhl_iface) * 7.5

    return n_contacts, iface_area_est, n_met_iface, n_vhl_iface


POLAR_NAMES = ['N', 'O', 'ND1', 'ND2', 'NE', 'NE2', 'NH1', 'NH2',
               'NZ', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH', 'SG']

def hbond_count(atoms, d_cutoff=3.5):
    """
    Count N/O donor–acceptor pairs < d_cutoff between MET and the full VHL/ElonginB/C.
    No H-atom angle check (upper bound).
    """
    met_don = select(atoms, chain=MET_CHAIN, resi_range=(MET_MIN, MET_MAX),
                     name=POLAR_NAMES)
    # VHL = chain A (1-103) + chain B (ElonginB) + chain C (ElonginC)
    vhl_a  = select(atoms, chain=VHL_CHAIN, resi_range=(VHL_MIN, VHL_MAX),
                    name=POLAR_NAMES)
    vhl_bc = [a for a in atoms
              if a['chain'] in (ELONG_B_CHAIN, 'C') and a['name'] in POLAR_NAMES]
    vhl_acc = vhl_a + vhl_bc
    if not met_don or not vhl_acc:
        return 0
    mc = coords(met_don)
    vc = coords(vhl_acc)
    diff = mc[:, None, :] - vc[None, :, :]
    dists = np.sqrt((diff**2).sum(axis=2))
    return int((dists < d_cutoff).sum())


def protac_geometry(atoms):
    """
    Check PROTAC (chain X) for geometry issues.
    - Fragmentation: any atom whose nearest PROTAC neighbor is > 3 Å = disconnected piece
    - Steric overlap: any PROTAC atom pair < 1.3 Å (catastrophic overlap)
    Returns: (n_atoms, min_nearest_neighbor, n_disconnected, n_overlaps, centroid)
    """
    ptc = select(atoms, chain=PROTAC_CHAIN)
    if not ptc:
        return 0, 0.0, 0, 0, None
    co = coords(ptc)
    n = len(co)

    # Pairwise distances (excluding self)
    diff = co[:, None, :] - co[None, :, :]
    pdists = np.sqrt((diff**2).sum(axis=2))
    np.fill_diagonal(pdists, 999.0)

    # Nearest neighbor for each atom
    nearest = pdists.min(axis=1)
    max_nearest = float(nearest.max())         # >3 Å = likely disconnected fragment
    n_disconnected = int((nearest > 3.0).sum())

    # Overlaps (atoms < 1.3 Å from another PROTAC atom, can't be bonded)
    np.fill_diagonal(pdists, 999.0)
    n_overlaps = int((pdists < 1.3).sum()) // 2

    centroid = co.mean(axis=0)
    return n, max_nearest, n_disconnected, n_overlaps, centroid


def exposed_lysines_near_vhl(atoms, sasa_proxy_cutoff=10.0, near_vhl_cutoff=30.0):
    """
    Find MET lysines whose NZ is:
    1. 'Exposed' (proxy: no non-PROTAC atom within sasa_proxy_cutoff Å within the same chain)
       Actually we use: NZ with few neighbors (< 5 heavy atoms within 5 Å on MET side)
    2. Within near_vhl_cutoff Å of any VHL chain A atom

    Returns list of (resi, NZ_xyz, dist_to_vhl, n_buried_neighbors)
    """
    met_lys_nz = [a for a in atoms
                  if a['chain'] == MET_CHAIN and a['resi'] >= MET_MIN
                  and a['resn'] == 'LYS' and a['name'] == 'NZ']
    vhl_atoms = select(atoms, chain=VHL_CHAIN, resi_range=(VHL_MIN, VHL_MAX))
    met_heavy = select(atoms, chain=MET_CHAIN, resi_range=(MET_MIN, MET_MAX))
    elonb = [a for a in atoms if a['chain'] == ELONG_B_CHAIN]
    all_vhl = vhl_atoms + elonb

    if not met_lys_nz or not all_vhl:
        return []

    vhl_co = coords(all_vhl)
    met_co = coords(met_heavy)

    results = []
    for lys in met_lys_nz:
        nz = lys['xyz']

        # Burial proxy: count MET heavy atoms within 5 Å of NZ
        d_met = np.sqrt(((met_co - nz)**2).sum(axis=1))
        n_buried = int((d_met < 5.0).sum()) - 1  # -1 for self

        # Distance to nearest VHL atom
        d_vhl = np.sqrt(((vhl_co - nz)**2).sum(axis=1)).min()

        results.append((lys['resi'], nz, float(d_vhl), n_buried))

    # Sort by distance to VHL (closest first)
    results.sort(key=lambda x: x[2])
    return results


def e2_facing_distances(atoms, e2_lys_optimal=50.0):
    """
    Proxy for E2-facing geometry:
    Use ElonginB Cα centroid as anchor for the CRL2 scaffold direction.
    The E2 active-site Cys is ~50-60 Å from VHL in the CRL2 complex.
    Return distances from MET Lys NZ to ElonginB centroid.
    """
    elonb_ca = select(atoms, chain=ELONG_B_CHAIN, name='CA')
    met_lys_nz = [a for a in atoms
                  if a['chain'] == MET_CHAIN and a['resi'] >= MET_MIN
                  and a['resn'] == 'LYS' and a['name'] == 'NZ']
    if not elonb_ca or not met_lys_nz:
        return []

    eb_centroid = coords(elonb_ca).mean(axis=0)
    results = []
    for lys in met_lys_nz:
        d = float(np.linalg.norm(lys['xyz'] - eb_centroid))
        results.append((lys['resi'], d))
    results.sort(key=lambda x: x[1])
    return results


def clash_check(atoms, clash_cutoff=2.0):
    """
    Count heavy-atom pairs < clash_cutoff Å between MET and VHL, and PROTAC vs protein.
    H atoms are excluded — PyRosetta relax adds H to the combined PDB, and H-H
    contacts at <2 Å are not real steric clashes.
    """
    def _ha(lst):
        return [a for a in lst if a.get('elem', '') not in ('H', '')]

    met = _ha(select(atoms, chain=MET_CHAIN, resi_range=(MET_MIN, MET_MAX)))
    vhl_a = _ha(select(atoms, chain=VHL_CHAIN, resi_range=(VHL_MIN, VHL_MAX)))
    vhl_bc = _ha([a for a in atoms if a['chain'] in (ELONG_B_CHAIN, 'C')])
    vhl = vhl_a + vhl_bc
    ptc = _ha(select(atoms, chain=PROTAC_CHAIN))

    n_clash_interface = 0
    if met and vhl:
        mc = coords(met)
        vc = coords(vhl)
        diff = mc[:, None, :] - vc[None, :, :]
        dists = np.sqrt((diff**2).sum(axis=2))
        n_clash_interface = int((dists < clash_cutoff).sum())

    n_clash_protac = 0
    prot_atoms = _ha([a for a in atoms if a['chain'] != PROTAC_CHAIN])
    if ptc and prot_atoms:
        pc = coords(ptc)
        prco = coords(prot_atoms)
        diff2 = pc[:, None, :] - prco[None, :, :]
        dists2 = np.sqrt((diff2**2).sum(axis=2))
        n_clash_protac = int((dists2 < clash_cutoff).sum())

    return n_clash_interface, n_clash_protac


# ── Main analysis ─────────────────────────────────────────────────────────────

def pick_representative(cluster_dir, scores):
    """Pick best-scored PDB in a cluster directory."""
    pdbs = sorted(glob.glob(str(cluster_dir / 'combined_*_0001.pdb')))
    if not pdbs:
        return None
    best_pdb, best_score = pdbs[0], float('inf')
    for p in pdbs:
        desc = Path(p).stem          # combined_48_10_0001
        sc = scores.get(desc, float('inf'))
        if sc < best_score:
            best_score = sc
            best_pdb = p
    return best_pdb, best_score


def analyse_cluster(cluster_id, cluster_dir, scores):
    rep, rep_score = pick_representative(cluster_dir, scores)
    if rep is None:
        return None

    atoms = parse_pdb(rep)
    n_poses = len(list(cluster_dir.glob('combined_*_0001.pdb')))

    n_contacts, iface_area, n_met_i, n_vhl_i = interface_contacts(atoms)
    n_hbonds  = hbond_count(atoms)
    n_ptc, max_nn, n_disconn, n_overlaps, ptc_cen = protac_geometry(atoms)
    lys_list  = exposed_lysines_near_vhl(atoms)
    e2_dists  = e2_facing_distances(atoms)
    n_clash_if, n_clash_ptc = clash_check(atoms)

    # Lysines within 25 Å of VHL (accessible)
    close_lys = [(r, d) for r, _, d, nb in lys_list if d < 25.0 and nb < 8]
    # Optimal E2 candidates: 40-70 Å from ElonginB centroid
    e2_candidates = [(r, d) for r, d in e2_dists if 35.0 < d < 75.0]

    return {
        'cluster': cluster_id,
        'n_poses': n_poses,
        'representative': Path(rep).name,
        'score_orig': rep_score,
        'interface_contacts_5A': n_contacts,
        'interface_area_est_A2': round(iface_area),
        'met_iface_atoms': n_met_i,
        'vhl_iface_atoms': n_vhl_i,
        'polar_hbond_pairs': n_hbonds,
        'protac_n_atoms': n_ptc,
        'protac_max_nn_dist': round(max_nn, 2),
        'protac_disconnected': n_disconn,
        'protac_overlaps': n_overlaps,
        'clash_interface': n_clash_if,
        'clash_protac_protein': n_clash_ptc,
        'n_lys_within_25A_vhl': len(close_lys),
        'closest_lys_to_vhl': close_lys[:3],          # top-3 closest
        'n_lys_e2_optimal_range': len(e2_candidates),
        'best_e2_lys': e2_candidates[:3],              # top-3 optimal
    }


def print_report(features):
    print('\n' + '='*80)
    print(f"CLUSTER {features['cluster']}  ({features['n_poses']} poses)")
    print(f"  Representative : {features['representative']}")
    print(f"  Rosetta score  : {features['score_orig']:.1f} REU (pre-shift)")
    print()
    print(f"  Interface:")
    print(f"    Contacts (≤5 Å)    : {features['interface_contacts_5A']}")
    print(f"    Est. area          : {features['interface_area_est_A2']} Å²")
    print(f"    MET iface atoms    : {features['met_iface_atoms']}")
    print(f"    VHL iface atoms    : {features['vhl_iface_atoms']}")
    print(f"    Polar H-bond pairs : {features['polar_hbond_pairs']}")
    print()
    print(f"  PROTAC geometry:")
    print(f"    Atoms             : {features['protac_n_atoms']}")
    nn = features['protac_max_nn_dist']
    print(f"    Max nearest-neigh : {nn:.2f} Å "
          f"{'⚠ FRAGMENT?' if nn > 3.0 else '✓'}")
    print(f"    Disconnected atoms: {features['protac_disconnected']} "
          f"{'⚠' if features['protac_disconnected'] > 0 else '✓'}")
    print(f"    Atom overlaps     : {features['protac_overlaps']}")
    print()
    print(f"  Clashes:")
    print(f"    MET–VHL (<2Å)     : {features['clash_interface']} "
          f"{'⚠ BAD' if features['clash_interface'] > 5 else '✓'}")
    print(f"    PROTAC–protein    : {features['clash_protac_protein']} "
          f"{'⚠ BAD' if features['clash_protac_protein'] > 10 else '✓'}")
    print()
    print(f"  Ubiquitination geometry:")
    print(f"    MET Lys within 25Å of VHL : {features['n_lys_within_25A_vhl']}")
    if features['closest_lys_to_vhl']:
        lys_str = ', '.join(f"K{r}({d:.0f}Å)" for r, d in features['closest_lys_to_vhl'])
        print(f"    Closest MET Lys to VHL   : {lys_str}")
    print(f"    MET Lys in E2 range 35-75Å: {features['n_lys_e2_optimal_range']}")
    if features['best_e2_lys']:
        e2_str = ', '.join(f"K{r}({d:.0f}Å)" for r, d in features['best_e2_lys'])
        print(f"    Best E2-accessible Lys   : {e2_str}")


def main():
    scores = load_scores()
    all_features = []

    for cluster_dir in sorted(RESULTS.glob('cluster*')):
        cid = cluster_dir.name
        feat = analyse_cluster(cid, cluster_dir, scores)
        if feat:
            all_features.append(feat)
            print_report(feat)

    # ── MD candidate selection ─────────────────────────────────────────────
    print('\n' + '='*80)
    print('MD CANDIDATE SELECTION')
    print('='*80)

    # Criterion 1: largest cluster (most populated = most sampled, robust)
    cand1 = max(all_features, key=lambda x: x['n_poses'])

    # Filter to structurally valid (not massive clashes, not broken PROTAC)
    valid = [f for f in all_features
             if f['clash_interface'] <= 20
             and f['clash_protac_protein'] <= 50
             and f['protac_disconnected'] == 0]

    # Criterion 2: best Rosetta score among valid, excluding cand1
    others = [f for f in valid if f['cluster'] != cand1['cluster']]
    cand2 = min(others, key=lambda x: x['score_orig']) if others else None

    # Criterion 3: best Lys accessibility (most E2-optimal lysines, structurally distinct)
    already = {cand1['cluster']}
    if cand2:
        already.add(cand2['cluster'])
    lys_candidates = [f for f in valid if f['cluster'] not in already]
    cand3 = max(lys_candidates, key=lambda x: x['n_lys_e2_optimal_range']) if lys_candidates else None

    candidates = [c for c in [cand1, cand2, cand3] if c]

    for i, c in enumerate(candidates):
        reasons = []
        if c == cand1:
            reasons.append(f'largest cluster ({c["n_poses"]} poses)')
        if c == cand2:
            reasons.append(f'best Rosetta score ({c["score_orig"]:.0f} REU)')
        if c == cand3:
            reasons.append(f'most E2-optimal Lys ({c["n_lys_e2_optimal_range"]} in range)')
        print(f'\n  MD candidate {i+1}: {c["cluster"].upper()}')
        print(f'    Reason     : {"; ".join(reasons)}')
        print(f'    Rep struct : {c["representative"]}')
        print(f'    Score      : {c["score_orig"]:.0f} REU')
        print(f'    Interface  : {c["interface_contacts_5A"]} contacts, '
              f'{c["interface_area_est_A2"]} Å², {c["polar_hbond_pairs"]} polar pairs')
        print(f'    Clashes    : interface={c["clash_interface"]}, '
              f'PROTAC-prot={c["clash_protac_protein"]}')
        print(f'    E2 Lys     : {c["n_lys_e2_optimal_range"]} in 35-75Å range')
        if c['best_e2_lys']:
            print(f'    Best Lys   : ' +
                  ', '.join(f"K{r}({d:.0f}Å)" for r, d in c['best_e2_lys']))

    # ── Save outputs ───────────────────────────────────────────────────────────
    out_dir = RESULTS.parent   # MET run directory
    _save_outputs(all_features, candidates, out_dir)

    print('\n')
    return all_features, candidates


def _save_outputs(all_features, candidates, out_dir):
    """Write CSV feature table and JSON summary to the MET run directory."""
    out_dir = Path(out_dir)

    # Flat CSV (one row per cluster)
    csv_cols = [
        'cluster', 'n_poses', 'representative', 'score_orig_REU',
        'interface_contacts_5A', 'interface_area_est_A2',
        'met_iface_atoms', 'vhl_iface_atoms', 'polar_hbond_pairs',
        'protac_disconnected', 'protac_overlaps',
        'clash_interface', 'clash_protac_protein',
        'n_lys_within_25A_vhl', 'n_lys_e2_optimal_range',
        'best_e2_lys_top3', 'md_candidate',
    ]
    cand_clusters = {c['cluster']: i + 1 for i, c in enumerate(candidates)}

    csv_path = out_dir / 'cluster_features.csv'
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=csv_cols)
        w.writeheader()
        for feat in all_features:
            e2_str = ';'.join(f"K{r}({d:.0f}A)" for r, d in feat['best_e2_lys'])
            w.writerow({
                'cluster':                  feat['cluster'],
                'n_poses':                  feat['n_poses'],
                'representative':           feat['representative'],
                'score_orig_REU':           round(feat['score_orig'], 1),
                'interface_contacts_5A':    feat['interface_contacts_5A'],
                'interface_area_est_A2':    feat['interface_area_est_A2'],
                'met_iface_atoms':          feat['met_iface_atoms'],
                'vhl_iface_atoms':          feat['vhl_iface_atoms'],
                'polar_hbond_pairs':        feat['polar_hbond_pairs'],
                'protac_disconnected':      feat['protac_disconnected'],
                'protac_overlaps':          feat['protac_overlaps'],
                'clash_interface':          feat['clash_interface'],
                'clash_protac_protein':     feat['clash_protac_protein'],
                'n_lys_within_25A_vhl':     feat['n_lys_within_25A_vhl'],
                'n_lys_e2_optimal_range':   feat['n_lys_e2_optimal_range'],
                'best_e2_lys_top3':         e2_str,
                'md_candidate':             cand_clusters.get(feat['cluster'], ''),
            })
    print(f'[saved] {csv_path}')

    # Full JSON (includes all nested lists)
    json_path = out_dir / 'cluster_features.json'
    with open(json_path, 'w') as f:
        json.dump({
            'poi': 'MET',
            'clusters': all_features,
            'md_candidates': [
                {'rank': i + 1, 'cluster': c['cluster'],
                 'representative': c['representative']}
                for i, c in enumerate(candidates)
            ],
        }, f, indent=2)
    print(f'[saved] {json_path}')


if __name__ == '__main__':
    all_feat, candidates = main()
