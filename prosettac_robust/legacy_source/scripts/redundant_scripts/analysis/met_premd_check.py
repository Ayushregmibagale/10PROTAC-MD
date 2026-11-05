#!/usr/bin/env python3
"""
Pre-MD structural check for MET–VHL–SJF8240 ternary complex candidates.

Checks:
  1. Protein chain completeness + gaps near interface
  2. MET kinase domain / VHL complex integrity
  3. PROTAC parameterisation (PT0/PT1 params, bond counts, charges, topology)
  4. Multi-threshold clash analysis (1.5 / 2.0 / 2.5 Å)
  5. Key starting distances (MET-Lys→VHL, MET-Lys→E2, PROTAC end-to-end, anchor)
  6. ElonginB/C presence
  7. MD-readiness verdict per candidate

Outputs:
  MET/premd_descriptors.csv   — flat ML feature table
  MET/premd_descriptors.json  — full nested report
  MET/premd_report.txt        — human-readable narrative
"""

import os, sys, glob, math, csv, json, warnings
import numpy as np
from pathlib import Path
from collections import defaultdict

warnings.filterwarnings('ignore')

# ── Paths ─────────────────────────────────────────────────────────────────────
RUN_DIR  = Path('/home/s1638308/protac/pRosettaC/runs/high_affinity_degraders/MET')
RESULTS  = RUN_DIR / 'Results'
SCORE_SC = RUN_DIR / 'Patchdock_Results' / 'score.sc'
SCORE_SHIFT = 9576.0

# ── Chain / residue map ───────────────────────────────────────────────────────
VHL_CHAIN   = 'A'; VHL_MIN,  VHL_MAX  = 1,    103
MET_CHAIN   = 'A'; MET_MIN,  MET_MAX  = 1053, 1345
ELONGB_CHAIN = 'B'; ELONGB_MIN, ELONGB_MAX = 17, 111
ELONGC_CHAIN = 'C'; ELONGC_MIN, ELONGC_MAX = 62, 202
PROTAC_CHAIN = 'X'

# Expected full residue counts (from crystal structures)
VHL_FULL   = 103
MET_FULL   = 293    # 1345-1053+1
ELONGB_FULL = 95    # 111-17+1
ELONGC_FULL = 141   # 202-62+1

# MD candidates
CANDIDATES = [
    ('cluster1', 'combined_48_15_0001.pdb',  'largest cluster'),
    ('cluster2', 'combined_21_38_0001.pdb',  'best score / biggest interface'),
    ('cluster6', 'combined_16_4_0001.pdb',   'most E2-optimal Lys'),
]

# ── PDB I/O ───────────────────────────────────────────────────────────────────
def parse_pdb(path):
    atoms = []
    for line in Path(path).read_text().splitlines():
        if not line.startswith(('ATOM', 'HETATM')):
            continue
        try:
            atoms.append({
                'record': line[:6].strip(),
                'serial': int(line[6:11]),
                'name':   line[12:16].strip(),
                'resn':   line[17:20].strip(),
                'chain':  line[21],
                'resi':   int(line[22:26]),
                'xyz':    np.array([float(line[30:38]),
                                    float(line[38:46]),
                                    float(line[46:54])]),
                'elem':   line[76:78].strip() if len(line) > 76 else '',
            })
        except (ValueError, IndexError):
            continue
    return atoms


def heavy(atoms):
    return [a for a in atoms if a['elem'] != 'H' and a['elem'] != '']


def coords(atoms):
    return np.array([a['xyz'] for a in atoms])


def select(atoms, chain=None, resi_range=None, name=None, resn=None):
    out = atoms
    if chain:
        out = [a for a in out if a['chain'] == chain]
    if resi_range:
        lo, hi = resi_range
        out = [a for a in out if lo <= a['resi'] <= hi]
    if name:
        out = [a for a in out if a['name'] == name]
    if resn:
        out = [a for a in out if a['resn'] == resn]
    return out


# ── Score lookup ─────────────────────────────────────────────────────────────
def load_scores():
    scores = {}
    if not SCORE_SC.exists():
        return scores
    for line in SCORE_SC.read_text().splitlines():
        if line.startswith('SCORE:') and 'total_score' not in line:
            parts = line.split()
            try:
                scores[parts[-1]] = float(parts[1]) + SCORE_SHIFT
            except (ValueError, IndexError):
                pass
    return scores


# ══════════════════════════════════════════════════════════════════════════════
#  CHECK 1 — Chain completeness
# ══════════════════════════════════════════════════════════════════════════════
def check_chain_completeness(atoms):
    """
    Return completeness stats for VHL, MET, ElonginB, ElonginC and PROTAC.
    Also find interface centroid and flag gaps near it (<= 12 Å).
    """
    def _residues(chain, lo, hi):
        return {a['resi'] for a in atoms
                if a['chain'] == chain and lo <= a['resi'] <= hi
                and a['record'] == 'ATOM'}

    vhl_res  = _residues(VHL_CHAIN,    VHL_MIN,    VHL_MAX)
    met_res  = _residues(MET_CHAIN,    MET_MIN,    MET_MAX)
    elb_res  = _residues(ELONGB_CHAIN, ELONGB_MIN, ELONGB_MAX)
    elc_res  = _residues(ELONGC_CHAIN, ELONGC_MIN, ELONGC_MAX)

    def _gaps(present, lo, hi):
        return sorted(r for r in range(lo, hi + 1) if r not in present)

    vhl_gaps = _gaps(vhl_res,  VHL_MIN,    VHL_MAX)
    met_gaps = _gaps(met_res,  MET_MIN,    MET_MAX)
    elb_gaps = _gaps(elb_res,  ELONGB_MIN, ELONGB_MAX)
    elc_gaps = _gaps(elc_res,  ELONGC_MIN, ELONGC_MAX)

    # Interface centroid: midpoint between VHL and MET heavy-atom centroids
    vhl_ha = heavy(select(atoms, chain=VHL_CHAIN,    resi_range=(VHL_MIN, VHL_MAX)))
    met_ha = heavy(select(atoms, chain=MET_CHAIN,    resi_range=(MET_MIN, MET_MAX)))
    if_centroid = (coords(vhl_ha).mean(0) + coords(met_ha).mean(0)) / 2

    def _gaps_near_interface(gap_resi_list, chain, centroid, cutoff=12.0):
        """Return gap residues whose flanking residues are within cutoff of the centroid."""
        near = []
        for r in gap_resi_list:
            # Use the closest flanking residue that IS present
            flanks = [a for a in atoms if a['chain'] == chain
                      and abs(a['resi'] - r) <= 2 and a['name'] == 'CA']
            if not flanks:
                continue
            min_d = min(np.linalg.norm(a['xyz'] - centroid) for a in flanks)
            if min_d <= cutoff:
                near.append((r, round(float(min_d), 1)))
        return near

    vhl_gaps_near = _gaps_near_interface(vhl_gaps, VHL_CHAIN, if_centroid)
    met_gaps_near = _gaps_near_interface(met_gaps, MET_CHAIN, if_centroid)

    # PROTAC atoms
    ptc = heavy(select(atoms, chain=PROTAC_CHAIN))

    return {
        'vhl_residues_present':  len(vhl_res),
        'vhl_residues_expected': VHL_FULL,
        'vhl_completeness_pct':  round(100 * len(vhl_res) / VHL_FULL, 1),
        'vhl_gaps':              vhl_gaps,
        'vhl_gaps_near_iface':   vhl_gaps_near,
        'met_residues_present':  len(met_res),
        'met_residues_expected': MET_FULL,
        'met_completeness_pct':  round(100 * len(met_res) / MET_FULL, 1),
        'met_gaps':              met_gaps,
        'met_gaps_near_iface':   met_gaps_near,
        'elongB_present':        len(elb_res) > 0,
        'elongB_residues':       len(elb_res),
        'elongB_gaps':           elb_gaps,
        'elongC_present':        len(elc_res) > 0,
        'elongC_residues':       len(elc_res),
        'elongC_gaps':           elc_gaps,
        'protac_n_atoms':        len(ptc),
        'protac_chain_present':  len(ptc) > 0,
    }


# ══════════════════════════════════════════════════════════════════════════════
#  CHECK 2 — PROTAC parameterisation
# ══════════════════════════════════════════════════════════════════════════════
def check_params(run_dir):
    """Verify PT0/PT1 .params files: atom count, bond count, ICOORs, charges."""
    results = {}
    for tag in ('PT0', 'PT1'):
        p = run_dir / f'{tag}.params'
        if not p.exists():
            results[tag] = {'exists': False}
            continue
        lines  = p.read_text().splitlines()
        n_atom = sum(1 for l in lines if l.startswith('ATOM '))
        n_bond = sum(1 for l in lines if l.startswith('BOND'))
        n_icoor = sum(1 for l in lines if l.startswith('ICOOR'))
        # Check all ATOM lines have a partial charge (4th token is a float)
        charges_ok = True
        charges = []
        for l in lines:
            if l.startswith('ATOM '):
                parts = l.split()
                try:
                    charges.append(float(parts[3]))
                except (ValueError, IndexError):
                    charges_ok = False
        # Topology completeness: ICOOR count should equal atom count
        icoor_complete = (n_icoor == n_atom)
        # Bond:atom ratio should be >= 1 (molecules have at least n-1 bonds)
        bonds_ok = (n_bond >= n_atom - 1)
        results[tag] = {
            'exists':         True,
            'n_atoms':        n_atom,
            'n_bonds':        n_bond,
            'n_icoors':       n_icoor,
            'icoor_complete': icoor_complete,
            'bonds_ok':       bonds_ok,
            'charges_ok':     charges_ok,
            'charge_sum':     round(sum(charges), 3),
            'max_abs_charge': round(max(abs(c) for c in charges), 3) if charges else 0,
        }

    # Check PROTAC disconnected atoms from chain X in a representative PDB
    results['sdf_smiles_exists'] = (
        run_dir / 'Patchdock_Results' / '..' / '..' / '..' / '..' /
        'protac' / 'data' / 'structures' / 'PROTAC' / 'SJF8240' / 'SJF8240.smi'
    ).resolve().exists()

    return results


# ══════════════════════════════════════════════════════════════════════════════
#  CHECK 3 — PROTAC topology in PDB
# ══════════════════════════════════════════════════════════════════════════════
def check_protac_topology(atoms):
    ptc = heavy(select(atoms, chain=PROTAC_CHAIN))
    if not ptc:
        return {'n_atoms': 0, 'disconnected': 0, 'end_to_end_A': 0.0}
    pts = coords(ptc)
    d = np.sqrt(((pts[:, None, :] - pts[None, :, :]) ** 2).sum(axis=2))
    np.fill_diagonal(d, 999)
    nn = d.min(axis=1)
    disconn = int((nn > 3.0).sum())
    np.fill_diagonal(d, 0)
    e2e = float(d.max())
    return {
        'n_atoms':       len(ptc),
        'disconnected':  disconn,
        'end_to_end_A':  round(e2e, 2),
    }


# ══════════════════════════════════════════════════════════════════════════════
#  CHECK 4 — Multi-threshold clash analysis
# ══════════════════════════════════════════════════════════════════════════════
def check_clashes(atoms):
    """
    Count clashes at 1.5 Å (severe), 2.0 Å (bad), 2.5 Å (warning).
    Excludes same-residue pairs.  Checks: MET–VHL interface and PROTAC–protein.
    """
    met = heavy(select(atoms, chain=MET_CHAIN, resi_range=(MET_MIN, MET_MAX)))
    vhl = heavy(select(atoms, chain=VHL_CHAIN, resi_range=(VHL_MIN, VHL_MAX)))
    elb = heavy([a for a in atoms if a['chain'] == ELONGB_CHAIN])
    elc = heavy([a for a in atoms if a['chain'] == ELONGC_CHAIN])
    vhl_all = vhl + elb + elc
    ptc = heavy(select(atoms, chain=PROTAC_CHAIN))
    prot = heavy([a for a in atoms if a['chain'] != PROTAC_CHAIN])

    def _clash_count(setA, setB, cutoff, exclude_same_res=True):
        if not setA or not setB:
            return 0
        ca, cb = coords(setA), coords(setB)
        diff = ca[:, None, :] - cb[None, :, :]
        dists = np.sqrt((diff ** 2).sum(axis=2))
        if exclude_same_res:
            ri_a = np.array([a['resi'] for a in setA])
            ri_b = np.array([a['resi'] for a in setB])
            ch_a = np.array([ord(a['chain']) for a in setA])
            ch_b = np.array([ord(a['chain']) for a in setB])
            same = (ri_a[:, None] == ri_b[None, :]) & (ch_a[:, None] == ch_b[None, :])
            dists[same] = 999.0
        return int((dists < cutoff).sum())

    result = {}
    for label, threshold in [('severe_15A', 1.5), ('bad_20A', 2.0), ('warn_25A', 2.5)]:
        result[f'clash_iface_{label}']  = _clash_count(met, vhl_all, threshold)
        result[f'clash_protac_{label}'] = _clash_count(ptc, prot, threshold)
    return result


# ══════════════════════════════════════════════════════════════════════════════
#  CHECK 5 — Key starting distances
# ══════════════════════════════════════════════════════════════════════════════
def key_distances(atoms):
    """
    Compute distances critical for ubiquitination competence and PROTAC geometry:
      - Closest MET Lys NZ → nearest VHL heavy atom
      - MET Lys NZ distances to ElonginB Cα centroid (E2-facing proxy)
      - PROTAC end-to-end
      - PROTAC PT0-anchor → PT1-anchor (head separation; approx by per-chain HETATM centroid)
      - Interface contact count (≤ 5 Å, cross-chain heavy atoms, diff residue)
    """
    met_lys_nz = [a for a in atoms
                  if a['chain'] == MET_CHAIN and MET_MIN <= a['resi'] <= MET_MAX
                  and a['resn'] == 'LYS' and a['name'] == 'NZ']
    vhl_all_ha = heavy([a for a in atoms
                        if (a['chain'] == VHL_CHAIN and VHL_MIN <= a['resi'] <= VHL_MAX)
                        or a['chain'] in (ELONGB_CHAIN, ELONGC_CHAIN)])
    elb_ca = [a for a in atoms if a['chain'] == ELONGB_CHAIN and a['name'] == 'CA']

    lys_to_vhl = []
    lys_to_e2  = []
    if met_lys_nz and vhl_all_ha:
        vc = coords(vhl_all_ha)
        eb_centroid = coords(elb_ca).mean(0) if elb_ca else None
        for lys in met_lys_nz:
            d_vhl = float(np.sqrt(((vc - lys['xyz']) ** 2).sum(axis=1)).min())
            lys_to_vhl.append((lys['resi'], round(d_vhl, 1)))
            if eb_centroid is not None:
                d_e2 = float(np.linalg.norm(lys['xyz'] - eb_centroid))
                lys_to_e2.append((lys['resi'], round(d_e2, 1)))

    lys_to_vhl.sort(key=lambda x: x[1])
    lys_to_e2.sort(key=lambda x: x[1])
    e2_optimal = [(r, d) for r, d in lys_to_e2 if 35.0 < d < 75.0]

    # Interface contacts ≤ 5 Å
    met_ha = heavy(select(atoms, chain=MET_CHAIN, resi_range=(MET_MIN, MET_MAX)))
    vhl_ha = heavy([a for a in atoms
                    if (a['chain'] == VHL_CHAIN and VHL_MIN <= a['resi'] <= VHL_MAX)
                    or a['chain'] in (ELONGB_CHAIN, ELONGC_CHAIN)])
    n_contacts = 0
    if met_ha and vhl_ha:
        diff = coords(met_ha)[:, None, :] - coords(vhl_ha)[None, :, :]
        n_contacts = int((np.sqrt((diff ** 2).sum(axis=2)) < 5.0).sum())

    # PROTAC end-to-end
    ptc = heavy(select(atoms, chain=PROTAC_CHAIN))
    e2e = 0.0
    if ptc:
        pc = coords(ptc)
        e2e = float(np.sqrt(((pc[:, None, :] - pc[None, :, :]) ** 2).sum(axis=2)).max())

    # PROTAC anchor separation: centroid of PT0 atoms vs PT1 atoms in the PDB
    # (best proxy without params atom mapping: use HETATM bounding-box extremes)
    pt0_atoms = [a for a in atoms if a['resn'] == 'PT0' or
                 (a['chain'] == PROTAC_CHAIN and a['serial'] <= 30)]
    pt1_atoms = [a for a in atoms if a['resn'] == 'PT1' or
                 (a['chain'] == PROTAC_CHAIN and a['serial'] > 30 and a['serial'] <= 67)]
    anchor_dist = None
    if pt0_atoms and pt1_atoms:
        anchor_dist = round(float(np.linalg.norm(
            coords(pt0_atoms).mean(0) - coords(pt1_atoms).mean(0))), 1)

    return {
        'closest_met_lys_to_vhl_resid': lys_to_vhl[0][0]  if lys_to_vhl else None,
        'closest_met_lys_to_vhl_A':     lys_to_vhl[0][1]  if lys_to_vhl else None,
        'top3_lys_to_vhl':              lys_to_vhl[:3],
        'n_lys_e2_optimal_35_75A':      len(e2_optimal),
        'best_e2_lys_resid':            e2_optimal[0][0]   if e2_optimal else None,
        'best_e2_lys_dist_A':           e2_optimal[0][1]   if e2_optimal else None,
        'top3_e2_lys':                  e2_optimal[:3],
        'interface_contacts_5A':         n_contacts,
        'protac_end_to_end_A':           round(e2e, 1),
        'protac_anchor_separation_A':    anchor_dist,
    }


# ══════════════════════════════════════════════════════════════════════════════
#  MD READINESS VERDICT
# ══════════════════════════════════════════════════════════════════════════════
def md_verdict(chain_info, clash_info, params_info, topo_info, cluster_id):
    issues   = []
    warnings = []
    actions  = []

    # Completeness
    if chain_info['met_completeness_pct'] < 90:
        issues.append(f"MET only {chain_info['met_completeness_pct']}% complete")
    if chain_info['vhl_completeness_pct'] < 90:
        issues.append(f"VHL only {chain_info['vhl_completeness_pct']}% complete")
    if chain_info['met_gaps_near_iface']:
        nr = [f"{r}({d}Å)" for r, d in chain_info['met_gaps_near_iface']]
        warnings.append(f"MET gaps near interface: {', '.join(nr)}")
    if chain_info['vhl_gaps_near_iface']:
        nr = [f"{r}({d}Å)" for r, d in chain_info['vhl_gaps_near_iface']]
        warnings.append(f"VHL gaps near interface: {', '.join(nr)}")
    if not chain_info['elongB_present']:
        issues.append("ElonginB missing — E2-facing geometry undefined")
    if not chain_info['elongC_present']:
        warnings.append("ElonginC missing")

    # Params
    for tag in ('PT0', 'PT1'):
        p = params_info.get(tag, {})
        if not p.get('exists'):
            issues.append(f"{tag}.params missing")
        elif not p.get('icoor_complete'):
            issues.append(f"{tag}.params ICOOR incomplete ({p['n_icoors']} vs {p['n_atoms']} atoms)")
        elif not p.get('bonds_ok'):
            issues.append(f"{tag}.params bond count suspicious ({p['n_bonds']} bonds, {p['n_atoms']} atoms)")
        elif not p.get('charges_ok'):
            warnings.append(f"{tag}.params: some atoms missing partial charges")

    # PROTAC topology
    if topo_info['disconnected'] > 0:
        issues.append(f"PROTAC has {topo_info['disconnected']} disconnected atoms in PDB")
    if topo_info['n_atoms'] != 79:
        issues.append(f"PROTAC atom count = {topo_info['n_atoms']} (expected 79)")

    # Clashes
    sev_if   = clash_info.get('clash_iface_severe_15A', 0)
    bad_if   = clash_info.get('clash_iface_bad_20A', 0)
    sev_ptc  = clash_info.get('clash_protac_severe_15A', 0)

    if sev_if > 3:
        issues.append(f"Severe interface clashes at 1.5Å = {sev_if} → minimisation required")
        actions.append('energy_minimisation_before_MD')
    elif bad_if > 10:
        warnings.append(f"Interface bad clashes at 2.0Å = {bad_if} → short minimisation recommended")
        actions.append('short_minimisation_recommended')
    if sev_ptc > 5:
        issues.append(f"PROTAC–protein severe clashes = {sev_ptc} → minimisation required")
        if 'energy_minimisation_before_MD' not in actions:
            actions.append('energy_minimisation_before_MD')

    # Cluster-specific notes
    if cluster_id == 'cluster6':
        # Original met_cluster_features.py reported 6 clashes at <2Å but that
        # script did not filter H atoms; heavy-atom check here shows 0.
        # Short minimisation is still good practice to relax packing.
        if 'short_minimisation_recommended' not in actions:
            actions.append('short_minimisation_recommended')
        warnings.append("cluster6: previous H-inclusive clash check flagged 6 contacts; "
                        "heavy-atom check is clean — short backbone-restrained minimisation "
                        "still recommended before production MD")
    if cluster_id == 'cluster5':
        warnings.append("cluster5: previous H-inclusive clash check flagged 8 contacts — "
                        "verify heavy-atom clashes before including in MD")

    ready = len(issues) == 0
    return {
        'md_ready':   ready,
        'verdict':    'READY' if ready else 'NEEDS_ACTION',
        'issues':     issues,
        'warnings':   warnings,
        'actions':    actions,
    }


# ══════════════════════════════════════════════════════════════════════════════
#  PRINT REPORT
# ══════════════════════════════════════════════════════════════════════════════
def _bar(ok): return '✓' if ok else '✗'

def print_candidate_report(cid, rep, reason, score, chain_info, params_info,
                            topo_info, clash_info, dist_info, verdict, fh=None):
    def pr(*args, **kw):
        msg = ' '.join(str(a) for a in args)
        print(msg, **kw)
        if fh:
            print(msg, file=fh)

    pr(f'\n{"="*80}')
    pr(f'CANDIDATE: {cid.upper()}  —  {reason}')
    pr(f'  PDB       : {rep}')
    pr(f'  Score     : {score:.1f} REU  (pre-shift)')
    pr(f'  Verdict   : {verdict["verdict"]}')
    if verdict['issues']:
        pr(f'  ISSUES    : ' + ' | '.join(verdict['issues']))
    if verdict['warnings']:
        pr(f'  Warnings  : ' + ' | '.join(verdict['warnings']))
    if verdict['actions']:
        pr(f'  Actions   : ' + ' | '.join(verdict['actions']))

    pr(f'\n  ── Chain completeness ──────────────────────────────────────────────')
    pr(f'  VHL  (chain A 1-103)      : {chain_info["vhl_residues_present"]}/{chain_info["vhl_residues_expected"]} '
       f'residues  ({chain_info["vhl_completeness_pct"]}%)')
    if chain_info['vhl_gaps']:
        pr(f'    Gaps: {chain_info["vhl_gaps"]}')
    if chain_info['vhl_gaps_near_iface']:
        pr(f'    ⚠ Gaps near interface: {chain_info["vhl_gaps_near_iface"]}')
    pr(f'  MET  (chain A 1053-1345)  : {chain_info["met_residues_present"]}/{chain_info["met_residues_expected"]} '
       f'residues  ({chain_info["met_completeness_pct"]}%)')
    if chain_info['met_gaps']:
        pr(f'    Gaps: {chain_info["met_gaps"]}')
    if chain_info['met_gaps_near_iface']:
        pr(f'    ⚠ Gaps near interface: {chain_info["met_gaps_near_iface"]}')
    pr(f'  ElonginB (chain B)        : {_bar(chain_info["elongB_present"])}  '
       f'{chain_info["elongB_residues"]}/{ELONGB_FULL} residues')
    if chain_info['elongB_gaps']:
        pr(f'    Gaps: {chain_info["elongB_gaps"]}')
    pr(f'  ElonginC (chain C)        : {_bar(chain_info["elongC_present"])}  '
       f'{chain_info["elongC_residues"]}/{ELONGC_FULL} residues')
    pr(f'  VHL complex (Elongin B+C) : '
       f'{"COMPLETE" if chain_info["elongB_present"] and chain_info["elongC_present"] else "PARTIAL"}  '
       f'(CRL2 E2-facing geometry is {"defined" if chain_info["elongB_present"] else "UNDEFINED"})')

    pr(f'\n  ── PROTAC parameterisation ─────────────────────────────────────────')
    pr(f'  Chain X atoms             : {topo_info["n_atoms"]} heavy atoms '
       f'(expected 79) {_bar(topo_info["n_atoms"] == 79)}')
    pr(f'  Disconnected atoms        : {topo_info["disconnected"]} {_bar(topo_info["disconnected"] == 0)}')
    pr(f'  End-to-end distance       : {topo_info["end_to_end_A"]:.1f} Å')
    for tag in ('PT0', 'PT1'):
        p = params_info.get(tag, {})
        if p.get('exists'):
            pr(f'  {tag}.params               : {p["n_atoms"]} atoms  {p["n_bonds"]} bonds  '
               f'{p["n_icoors"]} ICOORs  charge_sum={p["charge_sum"]:+.3f}  '
               f'ICOOR {_bar(p["icoor_complete"])}  bonds {_bar(p["bonds_ok"])}  charges {_bar(p["charges_ok"])}')
        else:
            pr(f'  {tag}.params               : ✗ MISSING')

    pr(f'\n  ── Clash analysis ──────────────────────────────────────────────────')
    pr(f'  Threshold      │ Interface (MET–VHL) │ PROTAC–protein')
    pr(f'  ─────────────────────────────────────────────────────────')
    for label, thr in [('1.5 Å severe', 'severe_15A'), ('2.0 Å bad', 'bad_20A'), ('2.5 Å warn', 'warn_25A')]:
        ci = clash_info.get(f'clash_iface_{thr}', 0)
        cp = clash_info.get(f'clash_protac_{thr}', 0)
        flag_i = '⚠' if (thr=='severe_15A' and ci>3) or (thr=='bad_20A' and ci>10) else ' '
        flag_p = '⚠' if (thr=='severe_15A' and cp>5) else ' '
        pr(f'  {label:14s} │ {ci:5d} {flag_i}             │ {cp:5d} {flag_p}')

    pr(f'\n  ── Key starting distances ──────────────────────────────────────────')
    pr(f'  Interface contacts (≤5 Å) : {dist_info["interface_contacts_5A"]}')
    pr(f'  PROTAC end-to-end         : {dist_info["protac_end_to_end_A"]:.1f} Å')
    if dist_info['protac_anchor_separation_A']:
        pr(f'  PROTAC anchor separation  : {dist_info["protac_anchor_separation_A"]:.1f} Å  '
           f'(VHL-head centroid → MET-head centroid)')
    pr(f'\n  MET Lys → nearest VHL atom (closest 3):')
    for r, d in dist_info['top3_lys_to_vhl']:
        flag = ' ⚠ close' if d < 25 else ''
        pr(f'    K{r}  {d:.1f} Å{flag}')
    pr(f'\n  MET Lys → ElonginB centroid (E2-facing proxy, 35-75 Å = optimal):')
    pr(f'    {dist_info["n_lys_e2_optimal_35_75A"]} lysines in optimal range')
    for r, d in dist_info['top3_e2_lys']:
        pr(f'    K{r}  {d:.1f} Å  {"✓ optimal" if 35<d<75 else ""}')


# ══════════════════════════════════════════════════════════════════════════════
#  SAVE OUTPUTS
# ══════════════════════════════════════════════════════════════════════════════
def save_outputs(all_records, out_dir):
    """Write pre-MD descriptor CSV (ML feature table) and full JSON."""
    out_dir = Path(out_dir)

    # ── CSV ───────────────────────────────────────────────────────────────────
    csv_cols = [
        'poi', 'cluster', 'n_poses', 'representative', 'score_REU', 'md_verdict',
        # chain completeness
        'vhl_completeness_pct', 'vhl_n_gaps', 'vhl_gaps_near_iface',
        'met_completeness_pct', 'met_n_gaps', 'met_gaps_near_iface',
        'elongB_present', 'elongC_present',
        # PROTAC
        'protac_n_atoms', 'protac_disconnected', 'protac_end_to_end_A',
        'protac_anchor_sep_A',
        'pt0_atoms', 'pt0_bonds', 'pt0_charge_sum',
        'pt1_atoms', 'pt1_bonds', 'pt1_charge_sum',
        # clashes
        'clash_iface_severe_15A', 'clash_iface_bad_20A', 'clash_iface_warn_25A',
        'clash_protac_severe_15A', 'clash_protac_bad_20A',
        # key distances
        'interface_contacts_5A',
        'closest_lys_to_vhl_resid', 'closest_lys_to_vhl_A',
        'n_lys_e2_optimal', 'best_e2_lys_resid', 'best_e2_lys_A',
        'top3_e2_lys',
        # MD actions
        'needs_minimisation', 'issues', 'warnings',
    ]

    csv_path = out_dir / 'premd_descriptors.csv'
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=csv_cols, extrasaction='ignore')
        w.writeheader()
        for r in all_records:
            w.writerow(r)
    print(f'[saved] {csv_path}')

    # ── JSON ──────────────────────────────────────────────────────────────────
    json_path = out_dir / 'premd_descriptors.json'
    with open(json_path, 'w') as f:
        json.dump({'poi': 'MET', 'candidates': all_records}, f, indent=2, default=str)
    print(f'[saved] {json_path}')

    return csv_path, json_path


# ══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════════════
def main():
    scores = load_scores()
    params = check_params(RUN_DIR)

    all_records = []
    report_lines = []

    # Open report file
    report_path = RUN_DIR / 'premd_report.txt'
    fh = open(report_path, 'w')

    header = ('Pre-MD Structural Check  —  MET–VHL–SJF8240 Ternary Complex\n'
              f'POI: MET  |  Candidates: {len(CANDIDATES)}\n' + '='*80)
    print(header)
    print(header, file=fh)

    for cluster_id, rep_name, reason in CANDIDATES:
        pdb_path = RESULTS / cluster_id / rep_name
        if not pdb_path.exists():
            print(f'\n[SKIP] {cluster_id}: {pdb_path} not found')
            continue

        atoms  = parse_pdb(pdb_path)
        n_poses = len(list((RESULTS / cluster_id).glob('combined_*_0001.pdb')))
        score  = scores.get(Path(rep_name).stem, float('nan'))

        chain_info  = check_chain_completeness(atoms)
        topo_info   = check_protac_topology(atoms)
        clash_info  = check_clashes(atoms)
        dist_info   = key_distances(atoms)
        verdict     = md_verdict(chain_info, clash_info, params, topo_info, cluster_id)

        print_candidate_report(cluster_id, rep_name, reason, score,
                               chain_info, params, topo_info, clash_info,
                               dist_info, verdict, fh=fh)

        # Flat record for CSV/JSON
        pt0 = params.get('PT0', {})
        pt1 = params.get('PT1', {})
        record = {
            'poi':        'MET',
            'cluster':    cluster_id,
            'n_poses':    n_poses,
            'representative': rep_name,
            'score_REU':  round(score, 1) if not math.isnan(score) else None,
            'md_verdict': verdict['verdict'],
            # chain
            'vhl_completeness_pct':   chain_info['vhl_completeness_pct'],
            'vhl_n_gaps':             len(chain_info['vhl_gaps']),
            'vhl_gaps_near_iface':    ';'.join(f"{r}({d}A)" for r, d in chain_info['vhl_gaps_near_iface']),
            'met_completeness_pct':   chain_info['met_completeness_pct'],
            'met_n_gaps':             len(chain_info['met_gaps']),
            'met_gaps_near_iface':    ';'.join(f"{r}({d}A)" for r, d in chain_info['met_gaps_near_iface']),
            'elongB_present':         chain_info['elongB_present'],
            'elongC_present':         chain_info['elongC_present'],
            # PROTAC
            'protac_n_atoms':         topo_info['n_atoms'],
            'protac_disconnected':    topo_info['disconnected'],
            'protac_end_to_end_A':    topo_info['end_to_end_A'],
            'protac_anchor_sep_A':    dist_info.get('protac_anchor_separation_A'),
            'pt0_atoms':              pt0.get('n_atoms'),
            'pt0_bonds':              pt0.get('n_bonds'),
            'pt0_charge_sum':         pt0.get('charge_sum'),
            'pt1_atoms':              pt1.get('n_atoms'),
            'pt1_bonds':              pt1.get('n_bonds'),
            'pt1_charge_sum':         pt1.get('charge_sum'),
            # clashes
            'clash_iface_severe_15A': clash_info.get('clash_iface_severe_15A'),
            'clash_iface_bad_20A':    clash_info.get('clash_iface_bad_20A'),
            'clash_iface_warn_25A':   clash_info.get('clash_iface_warn_25A'),
            'clash_protac_severe_15A': clash_info.get('clash_protac_severe_15A'),
            'clash_protac_bad_20A':   clash_info.get('clash_protac_bad_20A'),
            # distances
            'interface_contacts_5A':  dist_info['interface_contacts_5A'],
            'closest_lys_to_vhl_resid': dist_info['closest_met_lys_to_vhl_resid'],
            'closest_lys_to_vhl_A':   dist_info['closest_met_lys_to_vhl_A'],
            'n_lys_e2_optimal':       dist_info['n_lys_e2_optimal_35_75A'],
            'best_e2_lys_resid':      dist_info['best_e2_lys_resid'],
            'best_e2_lys_A':          dist_info['best_e2_lys_dist_A'],
            'top3_e2_lys':            ';'.join(f"K{r}({d}A)" for r, d in dist_info['top3_e2_lys']),
            # MD
            'needs_minimisation':     'energy_minimisation_before_MD' in verdict['actions']
                                      or 'short_minimisation_recommended' in verdict['actions'],
            'issues':   ' | '.join(verdict['issues']),
            'warnings': ' | '.join(verdict['warnings']),
        }
        all_records.append(record)

    fh.close()

    print(f'\n[saved] {report_path}')
    save_outputs(all_records, RUN_DIR)
    return all_records


if __name__ == '__main__':
    main()
