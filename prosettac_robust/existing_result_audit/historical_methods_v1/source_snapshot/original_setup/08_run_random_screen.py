#!/usr/bin/env python3
"""
run_random_screen.py — N=20 random bridging pose MD campaign
=============================================================
Samples 20 random bridging poses per POI from pRosettaC output and runs
20 ns MD on each, sequentially on one GPU. Fully resumable — can be stopped
and restarted at any time without losing work.

Pipeline per pose:
  STRUCTURE_DONE  → split pRosettaC combined PDB into complex_protein + PROTAC
  TOPOLOGY_DONE   → pdb2gmx + PROTAC merge + box + solvate + ions
  EM_DONE         → energy minimisation (CPU)
  NVT_DONE        → 200 ps NVT thermalisation (GPU, restrained)
  NPT1_DONE       → 500 ps NPT pressure equilibration (GPU, restrained)
  NPT2_DONE       → 500 ps NPT PROTAC release (GPU, protein restrained)
  PROD_DONE       → 20 ns production MD (GPU)
  ANALYSIS_DONE   → GROMACS analysis + metrics.csv + dashboard.png

Usage:
  python run_random_screen.py --plan              # sample poses, write manifest (dry run)
  python run_random_screen.py --run               # execute / resume campaign
  python run_random_screen.py --status            # show progress table
  python run_random_screen.py --run --poi DDR2    # run only one POI
  python run_random_screen.py --run --pose DDR2 pose003  # run one specific pose

Outputs:
  ~/protac-md/random_screen/manifest.tsv          — master state file
  ~/protac-md/random_screen/{POI}/pose{N:03d}/    — simulation files
  ~/protac-md/analysis/random_screen/{POI}/pose{N:03d}/ — analysis outputs
"""

import sys
import os
import json
import random
import shutil
import argparse
import textwrap
import subprocess
import tempfile
import time
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT         = Path.home() / 'protac-md'
RS_DIR       = ROOT / 'random_screen'          # simulation root
ANA_RS_DIR   = ROOT / 'analysis' / 'random_screen'
MANIFEST     = RS_DIR / 'manifest.tsv'
CAMPAIGN_LOG = RS_DIR / 'campaign.log'
MDP_DIR      = ROOT / 'scripts' / 'mdp'
LIGPAR       = ROOT / 'ligand_params' / 'SJF8240_3d.acpype'
SCRIPTS_DIR  = ROOT / 'scripts' / '01_setup'
ROSETTA_BASE = Path.home() / 'PROTAC degradation Prediction' / 'pRosettaC' / 'runs'
ANA_SCRIPTS  = ROOT / 'scripts' / '02_analysis'

# ── POI configuration ──────────────────────────────────────────────────────────
# vhl_max:  last VHL residue in chain A of combined PDB
# kin_min:  first kinase residue in chain A (becomes chain D after split)
# n_sample: number of random bridging poses to sample
#
# Excluded POIs:
#   CDK17, MAPK14-low — 0 real bridging poses (truly geometry-blocked)
#   TNIK, MAP4K5, SLK — kin_min=104: VHL (res 103) and kinase (res 104) are
#     covalently bonded in the combined PDB (1.33 Å gap = same chain), so pRosettaC
#     never performed genuine ternary docking. Bridging fraction 67–74% is a
#     structural artifact, not real geometry. Existing cherry-picked MDs (TNIK rank5
#     99.5%) remain valid but P_stable from random screen would be meaningless.
POI_CONFIG = {
    # ── Pulse-SILAC degraders ─────────────────────────────────────────────────
    'DDR2':        {'base': 'high_affinity_degraders',      'vhl_max': 103, 'kin_min': 553,  'n_sample': 20},
    'RIPK2_full':  {'base': 'high_affinity_degraders',      'vhl_max': 103, 'kin_min': 200,  'n_sample': 20},
    'MAPK14_full': {'base': 'high_affinity_degraders',      'vhl_max': 103, 'kin_min': 200,  'n_sample': 20},
    'MET':         {'base': 'high_affinity_degraders',      'vhl_max': 103, 'kin_min': 1053, 'n_sample': 20},
    'EPHB2':       {'base': 'high_affinity_degraders',      'vhl_max': 103, 'kin_min': 615,  'n_sample': 20},
    # ── Pulse-SILAC non-degraders ─────────────────────────────────────────────
    'ABL1':        {'base': 'high_affinity_no_degradation', 'vhl_max': 103, 'kin_min': 235,  'n_sample': 20},
    'AXL':         {'base': 'high_affinity_no_degradation', 'vhl_max': 103, 'kin_min': 521,  'n_sample': 20},
    'EPHA2':       {'base': 'high_affinity_no_degradation', 'vhl_max': 103, 'kin_min': 600,  'n_sample': 20},
    # ── Low-affinity non-degraders (chain D kinase — kin_min ignored by split_chains) ──
    'CDK2':        {'base': 'low_affinity_no_degradation',  'vhl_max': 103, 'kin_min': 1,    'n_sample': 20},
    'CDK4':        {'base': 'low_affinity_no_degradation',  'vhl_max': 103, 'kin_min': 1,    'n_sample': 20},
    'CDK6':        {'base': 'low_affinity_no_degradation',  'vhl_max': 103, 'kin_min': 1,    'n_sample': 20},
    'DAPK1':       {'base': 'low_affinity_no_degradation',  'vhl_max': 103, 'kin_min': 1,    'n_sample': 20},
    'IRAK1':       {'base': 'low_affinity_no_degradation',  'vhl_max': 103, 'kin_min': 1,    'n_sample': 20},
    'TLK2':        {'base': 'low_affinity_no_degradation',  'vhl_max': 103, 'kin_min': 1,    'n_sample': 20},
}

RANDOM_SEED  = 42
MAX_RETRIES  = 3
GENION_SOL_GROUP = 15   # Group 15 = SOL in our standard topology (verified from existing logs)

STAGES = [
    'PENDING',
    'STRUCTURE_DONE',
    'TOPOLOGY_DONE',
    'EM_DONE',
    'NVT_DONE',
    'NPT1_DONE',
    'NPT2_DONE',
    'PROD_DONE',
    'ANALYSIS_DONE',
    'FAILED',
]
STAGE_ORDER = {s: i for i, s in enumerate(STAGES)}


# ── Logging ───────────────────────────────────────────────────────────────────
def log(msg, also_print=True):
    ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    line = f'[{ts}] {msg}'
    if also_print:
        print(line, flush=True)
    CAMPAIGN_LOG.parent.mkdir(parents=True, exist_ok=True)
    with open(CAMPAIGN_LOG, 'a') as f:
        f.write(line + '\n')


# ── GROMACS runner ────────────────────────────────────────────────────────────
def gmx(args_str, stdin='', cwd=None, check=True, log_file=None):
    """
    Run a GROMACS command. Returns CompletedProcess.
    Raises RuntimeError on non-zero exit if check=True.
    Optionally tees stdout+stderr to log_file.
    """
    GMX_BIN = Path.home() / '.conda/envs/gromacs/bin/gmx'
    gmx_exe = str(GMX_BIN) if GMX_BIN.exists() else 'gmx'
    argv = [gmx_exe] + args_str.split()
    r = subprocess.run(argv, input=stdin, text=True,
                       capture_output=True, cwd=cwd)
    combined = r.stdout + r.stderr
    if log_file:
        Path(log_file).write_text(combined)
    if check and r.returncode != 0:
        raise RuntimeError(
            f'gmx {args_str[:60]} failed (rc={r.returncode}):\n'
            + combined[-600:]
        )
    return r


# ── Manifest management ───────────────────────────────────────────────────────
def read_manifest():
    if not MANIFEST.exists():
        return pd.DataFrame(columns=['poi','pose_id','source_pdb','stage',
                                      'attempts','last_updated','error'])
    return pd.read_csv(MANIFEST, sep='\t', dtype=str).fillna('')


def write_manifest(df):
    """Atomic write — write to temp file then rename."""
    tmp = MANIFEST.with_suffix('.tsv.tmp')
    df.to_csv(tmp, sep='\t', index=False)
    tmp.replace(MANIFEST)


def update_pose(df, poi, pose_id, stage=None, increment_attempts=False, error=''):
    mask = (df['poi'] == poi) & (df['pose_id'] == pose_id)
    if mask.sum() == 0:
        raise ValueError(f'Pose {poi}/{pose_id} not in manifest')
    idx = df.index[mask][0]
    if stage:
        df.at[idx, 'stage'] = stage
    if increment_attempts:
        df.at[idx, 'attempts'] = str(int(df.at[idx, 'attempts']) + 1)
    df.at[idx, 'last_updated'] = datetime.now().strftime('%Y-%m-%d %H:%M')
    df.at[idx, 'error'] = error
    write_manifest(df)
    return df


def get_pose_row(df, poi, pose_id):
    mask = (df['poi'] == poi) & (df['pose_id'] == pose_id)
    if mask.sum() == 0:
        return None
    return df[mask].iloc[0]


# ── Pose directory helpers ────────────────────────────────────────────────────
def pose_dir(poi, pose_id):
    return RS_DIR / poi / pose_id

def ana_dir(poi, pose_id):
    return ANA_RS_DIR / poi / pose_id

def topol_dir(poi, pose_id):
    return pose_dir(poi, pose_id) / 'topol'

def struct_dir(poi, pose_id):
    return pose_dir(poi, pose_id) / 'structure'

def em_dir(poi, pose_id):
    return pose_dir(poi, pose_id) / 'em'

def nvt_dir(poi, pose_id):
    return pose_dir(poi, pose_id) / 'nvt'

def npt1_dir(poi, pose_id):
    return pose_dir(poi, pose_id) / 'npt_restrained'

def npt2_dir(poi, pose_id):
    return pose_dir(poi, pose_id) / 'npt'

def prod_dir(poi, pose_id):
    return pose_dir(poi, pose_id) / 'prod'


# ── Stage verification (idempotent checks) ────────────────────────────────────
def is_structure_done(poi, pose_id):
    cp = struct_dir(poi, pose_id) / 'complex_protein.pdb'
    sp = struct_dir(poi, pose_id) / 'SJF8240_protac.pdb'
    if not (cp.exists() and sp.exists()):
        return False
    # Verify all four chains present in complex_protein.pdb
    content = cp.read_text()
    chains = {l[21] for l in content.splitlines() if l.startswith('ATOM')}
    return {'A', 'B', 'C', 'D'}.issubset(chains)

def is_topology_done(poi, pose_id):
    td = topol_dir(poi, pose_id)
    return (td / 'complex_ions.gro').exists() and (td / 'topol.top').exists()

def is_em_done(poi, pose_id):
    return (em_dir(poi, pose_id) / 'em.gro').exists()

def is_nvt_done(poi, pose_id):
    nd = nvt_dir(poi, pose_id)
    return (nd / 'nvt.gro').exists() and (nd / 'nvt.cpt').exists()

def is_npt1_done(poi, pose_id):
    return (npt1_dir(poi, pose_id) / 'npt_restrained.gro').exists()

def is_npt2_done(poi, pose_id):
    nd = npt2_dir(poi, pose_id)
    return (nd / 'npt.gro').exists() and (nd / 'npt.cpt').exists()

def is_prod_done(poi, pose_id):
    pd_ = prod_dir(poi, pose_id)
    return (pd_ / 'prod.gro').exists() and (pd_ / 'prod.xtc').exists()

def is_analysis_done(poi, pose_id):
    return (ana_dir(poi, pose_id) / 'metrics.csv').exists()

def prod_crashed_recoverable(poi, pose_id):
    """True if production started but crashed — checkpoint exists but no final GRO."""
    pd_ = prod_dir(poi, pose_id)
    return (pd_ / 'prod.cpt').exists() and not (pd_ / 'prod.gro').exists()

def detect_stage_from_disk(poi, pose_id):
    """Infer current stage by checking output files. Used on resume."""
    if is_analysis_done(poi, pose_id):  return 'ANALYSIS_DONE'
    if is_prod_done(poi, pose_id):      return 'PROD_DONE'
    if is_npt2_done(poi, pose_id):      return 'NPT2_DONE'
    if is_npt1_done(poi, pose_id):      return 'NPT1_DONE'
    if is_nvt_done(poi, pose_id):       return 'NVT_DONE'
    if is_em_done(poi, pose_id):        return 'EM_DONE'
    if is_topology_done(poi, pose_id):  return 'TOPOLOGY_DONE'
    if is_structure_done(poi, pose_id): return 'STRUCTURE_DONE'
    return 'PENDING'


# ── Plan: sample random bridging poses ────────────────────────────────────────
def plan_campaign():
    """
    Read bridging_screen.tsv for each POI, randomly sample n_sample bridging poses,
    and write the manifest. Idempotent — won't overwrite existing manifest entries.
    """
    RS_DIR.mkdir(parents=True, exist_ok=True)
    rng = random.Random(RANDOM_SEED)

    existing = read_manifest()
    new_rows = []

    for poi, cfg in POI_CONFIG.items():
        bscreen = ROOT / 'analysis' / poi / 'bridging_screen.tsv'
        if not bscreen.exists():
            log(f'WARNING: {poi} bridging_screen.tsv not found — skipping')
            continue

        df_screen = pd.read_csv(bscreen, sep='\t')

        mask = df_screen['bridging'] == 1

        bridging_poses = df_screen[mask]['pose'].tolist()
        n_bridging = len(bridging_poses)

        if n_bridging == 0:
            log(f'WARNING: {poi} has 0 real bridging poses — skipping')
            continue

        n_sample = min(cfg['n_sample'], n_bridging)
        sampled = rng.sample(bridging_poses, n_sample)

        log(f'{poi}: {n_bridging} real bridging poses → sampling {n_sample}')

        # Skip poses already in manifest
        existing_for_poi = existing[existing['poi'] == poi]['source_pdb'].tolist()

        for i, src_pdb in enumerate(sampled):
            pose_id = f'pose{i+1:03d}'
            if src_pdb in existing_for_poi:
                log(f'  {poi}/{pose_id} ({src_pdb}) already in manifest — skipping')
                continue
            new_rows.append({
                'poi':          poi,
                'pose_id':      pose_id,
                'source_pdb':   src_pdb,
                'stage':        'PENDING',
                'attempts':     '0',
                'last_updated': '',
                'error':        '',
            })

    if new_rows:
        new_df = pd.DataFrame(new_rows)
        combined = pd.concat([existing, new_df], ignore_index=True)
        write_manifest(combined)
        log(f'Manifest written: {MANIFEST}  ({len(new_rows)} new poses)')
    else:
        log('No new poses added (manifest already complete)')

    # Print summary
    df = read_manifest()
    print_status(df)


# ── Structure preparation ─────────────────────────────────────────────────────
def do_structure(poi, pose_id, source_pdb):
    """
    Copy pRosettaC combined PDB and split into complex_protein.pdb + SJF8240_protac.pdb.
    Replicates split_chains.py logic exactly (same chain reassignment rules).
    """
    cfg = POI_CONFIG[poi]
    vhl_max = cfg['vhl_max']
    kin_min = cfg['kin_min']

    pdb_src = (ROSETTA_BASE / cfg['base'] / poi / 'Patchdock_Results'
               / f'{source_pdb}.pdb')
    if not pdb_src.exists():
        raise FileNotFoundError(f'Source PDB not found: {pdb_src}')

    sd = struct_dir(poi, pose_id)
    sd.mkdir(parents=True, exist_ok=True)

    # Copy source for reference
    shutil.copy2(pdb_src, sd / 'source_combined.pdb')

    # Call existing split_chains.py
    r = subprocess.run(
        [sys.executable,
         str(SCRIPTS_DIR / '01_split_chains.py'),
         str(sd / 'source_combined.pdb'),
         str(vhl_max),
         str(kin_min)],
        capture_output=True, text=True, cwd=str(sd)
    )
    if r.returncode != 0:
        raise RuntimeError(f'split_chains.py failed:\n{r.stdout}\n{r.stderr}')

    # Verify outputs
    if not is_structure_done(poi, pose_id):
        raise RuntimeError(
            f'split_chains produced incomplete output for {poi}/{pose_id}.\n'
            f'stdout: {r.stdout[-300:]}'
        )

    log(f'  {poi}/{pose_id}: structure done ({source_pdb})')


# ── Topology build ────────────────────────────────────────────────────────────
def do_topology(poi, pose_id):
    """
    Replicate build_topology.sh exactly:
      1. pdb2gmx (AMBER99SB-ILDN + TIP3P)
      2. Patch topol.top with SJF8240.itp (atomtypes first — Bug 24)
      3. Reposition PROTAC to ternary complex coordinates + merge GRO
      4. editconf (dodecahedron box, 12 Å buffer)
      5. solvate (TIP3P)
      6. grompp + genion (0.15 M NaCl, neutralise)
    """
    sd = struct_dir(poi, pose_id)
    td = topol_dir(poi, pose_id)
    td.mkdir(parents=True, exist_ok=True)

    complex_pdb  = sd / 'complex_protein.pdb'
    protac_pdb   = sd / 'SJF8240_protac.pdb'
    acpype_itp   = LIGPAR / 'SJF8240_3d_GMX.itp'
    acpype_gro   = LIGPAR / 'SJF8240_3d_GMX.gro'

    for f in [complex_pdb, protac_pdb, acpype_itp, acpype_gro]:
        if not f.exists():
            raise FileNotFoundError(f'Required file missing: {f}')

    # ── 1. pdb2gmx ──────────────────────────────────────────────────────────
    gmx(f'pdb2gmx -f {complex_pdb} -o {td}/complex_protein.gro '
        f'-p {td}/topol.top -i {td}/posre.itp '
        f'-ff amber99sb-ildn -water tip3p -ignh -missing',
        log_file=td / 'pdb2gmx.log')

    # ── 2. Patch topol.top (Bug 24: atomtypes must precede moleculetype) ────
    shutil.copy2(acpype_itp, td / 'SJF8240.itp')
    shutil.copy2(acpype_gro, td / 'SJF8240.gro')

    top_text = (td / 'topol.top').read_text()
    ff_include = '#include "amber99sb-ildn.ff/forcefield.itp"'
    if ff_include not in top_text:
        raise RuntimeError('Expected forcefield include not found in topol.top')

    sjf_block = (
        '\n\n; GAFF2 atomtypes for SJF8240 — must precede all molecule definitions\n'
        '#include "SJF8240.itp"'
    )
    top_text = top_text.replace(ff_include, ff_include + sjf_block, 1)
    top_text = top_text.rstrip() + '\nSJF8240_3d              1\n'
    (td / 'topol.top').write_text(top_text)

    # ── 3. Reposition PROTAC to ternary complex coordinates ─────────────────
    _reposition_protac(td, sd / 'SJF8240_protac.pdb', poi)

    # ── 4. Box (dodecahedron, 12 Å buffer) ──────────────────────────────────
    gmx(f'editconf -f {td}/complex_merged.gro -o {td}/complex_box.gro '
        f'-bt dodecahedron -d 1.2',
        log_file=td / 'editconf.log')

    # ── 5. Solvate ───────────────────────────────────────────────────────────
    gmx(f'solvate -cp {td}/complex_box.gro -cs spc216.gro '
        f'-o {td}/complex_solv.gro -p {td}/topol.top',
        log_file=td / 'solvate.log')

    # ── 6. Ions: grompp then genion ──────────────────────────────────────────
    gmx(f'grompp -f {MDP_DIR}/ions.mdp -c {td}/complex_solv.gro '
        f'-p {td}/topol.top -o {td}/ions.tpr -maxwarn 2',
        log_file=td / 'grompp_ions.log')

    # Detect correct SOL group number (safety check against GENION_SOL_GROUP)
    sol_group = _detect_sol_group(td / 'ions.tpr')

    gmx(f'genion -s {td}/ions.tpr -o {td}/complex_ions.gro '
        f'-p {td}/topol.top -pname NA -nname CL -neutral -conc 0.15',
        stdin=f'{sol_group}\n',
        log_file=td / 'genion.log')

    if not is_topology_done(poi, pose_id):
        raise RuntimeError('Topology build completed but complex_ions.gro missing')

    log(f'  {poi}/{pose_id}: topology done')


def _reposition_protac(td, protac_pdb, poi):
    """
    Exact Python translation of build_topology.sh Step 3:
    Place heavy atoms at pRosettaC coordinates, H atoms parent-relative.
    Writes complex_merged.gro combining protein + repositioned PROTAC.
    """
    def parse_gro(fname):
        lines = Path(fname).read_text().splitlines()
        n = int(lines[1])
        atoms = []
        for line in lines[2:2+n]:
            name = line[10:15].strip()
            x = float(line[20:28]); y = float(line[28:36]); z = float(line[36:44])
            atoms.append({'name': name, 'x': x, 'y': y, 'z': z, 'raw': line})
        return lines[0], n, atoms, lines[2+n]

    _, n_prot, prot_atoms, box_line = parse_gro(td / 'complex_protein.gro')
    _, n_lig,  lig_atoms,  _        = parse_gro(td / 'SJF8240.gro')

    # Parse heavy atom coordinates from SJF8240_protac.pdb (Å → nm)
    pdb_coords = []
    for line in Path(protac_pdb).read_text().splitlines():
        if line[:4] not in ('ATOM', 'HETA'):
            continue
        pdb_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    pdb_coords = np.array(pdb_coords) / 10.0   # Å → nm

    heavy_idx = [i for i, a in enumerate(lig_atoms) if not a['name'].startswith('H')]
    h_idx     = [i for i, a in enumerate(lig_atoms) if     a['name'].startswith('H')]

    if len(heavy_idx) != len(pdb_coords):
        raise ValueError(
            f'Heavy atom count mismatch: GRO has {len(heavy_idx)}, '
            f'PDB has {len(pdb_coords)}'
        )

    # Parse H→heavy parent bonds from SJF8240.itp
    h_parent = {}
    in_bonds = False
    for line in (td / 'SJF8240.itp').read_text().splitlines():
        s = line.strip()
        if '[ bonds ]' in s or '[bonds]' in s:
            in_bonds = True; continue
        if in_bonds and s.startswith('['):
            in_bonds = False
        if in_bonds and s and not s.startswith(';'):
            parts = s.split()
            if len(parts) >= 2:
                a1, a2 = int(parts[0])-1, int(parts[1])-1
                n1 = lig_atoms[a1]['name']; n2 = lig_atoms[a2]['name']
                if n1.startswith('H') and not n2.startswith('H'):
                    h_parent[a1] = a2
                elif n2.startswith('H') and not n1.startswith('H'):
                    h_parent[a2] = a1

    gro_coords = np.array([[a['x'], a['y'], a['z']] for a in lig_atoms])
    new_coords = gro_coords.copy()

    # Heavy atoms: exact pRosettaC positions (0.0 Å error)
    for seq_i, gro_i in enumerate(heavy_idx):
        new_coords[gro_i] = pdb_coords[seq_i]

    # H atoms: translate relative to parent heavy atom (EM optimises geometry)
    missing_parent = 0
    for h_i in h_idx:
        parent = h_parent.get(h_i)
        if parent is not None:
            bond_vec = gro_coords[h_i] - gro_coords[parent]
            new_coords[h_i] = new_coords[parent] + bond_vec
        else:
            missing_parent += 1

    if missing_parent > 0:
        log(f'    WARNING: {missing_parent} H atoms had no parent bond — kept at acpype position')

    # Sanity: heavy atom placement error should be 0.0 Å
    placed = new_coords[heavy_idx]
    max_err = np.max(np.linalg.norm(placed - pdb_coords, axis=1)) * 10
    if max_err > 0.01:
        raise RuntimeError(f'Heavy atom placement error {max_err:.4f} Å > 0.01 Å threshold')

    # Sanity: PROTAC should be close to protein (< 10 Å)
    prot_coords = np.array([[a['x'], a['y'], a['z']] for a in prot_atoms])
    min_d = np.min(np.linalg.norm(
        prot_coords[:, None, :] - new_coords[None, :, :], axis=2)) * 10
    if min_d > 15.0:
        log(f'    WARNING: PROTAC–protein min dist = {min_d:.1f} Å (expect < 5 Å)')

    # Write complex_merged.gro
    total = n_prot + n_lig
    lines_out = [f'{poi} ternary complex with SJF8240\n', f'{total:5d}\n']
    lines_out.extend(a['raw'] + '\n' for a in prot_atoms)
    for i, a in enumerate(lig_atoms):
        x, y, z = new_coords[i]
        resnum_name = a['raw'][:15]
        atomnum     = a['raw'][15:20]
        lines_out.append(f'{resnum_name}{atomnum}{x:8.3f}{y:8.3f}{z:8.3f}'
                         f'{0.0:8.4f}{0.0:8.4f}{0.0:8.4f}\n')
    lines_out.append(box_line + '\n')
    (td / 'complex_merged.gro').write_text(''.join(lines_out))


def _detect_sol_group(tpr_file):
    """
    Run 'gmx genion -s tpr' with input 'q' to list groups, parse for SOL.
    Returns the group number for SOL (usually 15). Falls back to GENION_SOL_GROUP.
    """
    GMX_BIN = Path.home() / '.conda/envs/gromacs/bin/gmx'
    gmx_exe = str(GMX_BIN) if GMX_BIN.exists() else 'gmx'
    r = subprocess.run(
        [gmx_exe, 'genion', '-s', str(tpr_file)],
        input='q\n', text=True, capture_output=True
    )
    output = r.stdout + r.stderr
    for line in output.splitlines():
        # Lines like: "Group    15 (            SOL) has ..."
        if 'SOL' in line and 'Group' in line:
            parts = line.split()
            try:
                return int(parts[1])
            except (IndexError, ValueError):
                pass
    log(f'    WARNING: could not detect SOL group — using default {GENION_SOL_GROUP}')
    return GENION_SOL_GROUP


# ── Equilibration ─────────────────────────────────────────────────────────────
def do_em(poi, pose_id):
    """Energy minimisation (CPU). Idempotent."""
    td = topol_dir(poi, pose_id)
    ed = em_dir(poi, pose_id)
    ed.mkdir(parents=True, exist_ok=True)

    if (ed / 'em.gro').exists():
        log(f'  {poi}/{pose_id}: EM already done — skipping')
        return

    # Generate PROTAC position restraints if missing
    _ensure_posre_sjf8240(td)

    gmx(f'grompp -f {MDP_DIR}/em.mdp -c {td}/complex_ions.gro '
        f'-p {td}/topol.top -o {ed}/em.tpr -maxwarn 2',
        log_file=ed / 'grompp_em.log')

    gmx(f'mdrun -v -deffnm {ed}/em -ntmpi 1 -ntomp 8',
        log_file=ed / 'mdrun_em.log')

    if not (ed / 'em.gro').exists():
        raise RuntimeError(f'EM failed — em.gro not produced for {poi}/{pose_id}')

    log(f'  {poi}/{pose_id}: EM done')


def do_nvt(poi, pose_id):
    """200 ps NVT thermalisation, protein+PROTAC restrained."""
    td = topol_dir(poi, pose_id)
    ed = em_dir(poi, pose_id)
    nd = nvt_dir(poi, pose_id)
    nd.mkdir(parents=True, exist_ok=True)

    if (nd / 'nvt.gro').exists() and (nd / 'nvt.cpt').exists():
        log(f'  {poi}/{pose_id}: NVT already done — skipping')
        return

    gmx(f'grompp -f {MDP_DIR}/nvt.mdp -c {ed}/em.gro -r {ed}/em.gro '
        f'-p {td}/topol.top -o {nd}/nvt.tpr -maxwarn 2',
        log_file=nd / 'grompp_nvt.log')

    gmx(f'mdrun -deffnm {nd}/nvt -ntmpi 1 -ntomp 8 '
        f'-gpu_id 0 -nb gpu -pme gpu -bonded gpu -pin on',
        log_file=nd / 'mdrun_nvt.log')

    if not (nd / 'nvt.gro').exists():
        raise RuntimeError(f'NVT failed for {poi}/{pose_id}')

    log(f'  {poi}/{pose_id}: NVT done')


def do_npt1(poi, pose_id):
    """500 ps NPT, protein+PROTAC restrained."""
    td = topol_dir(poi, pose_id)
    nd = nvt_dir(poi, pose_id)
    n1 = npt1_dir(poi, pose_id)
    n1.mkdir(parents=True, exist_ok=True)

    if (n1 / 'npt_restrained.gro').exists():
        log(f'  {poi}/{pose_id}: NPT1 already done — skipping')
        return

    gmx(f'grompp -f {MDP_DIR}/npt_restrained.mdp '
        f'-c {nd}/nvt.gro -r {nd}/nvt.gro -t {nd}/nvt.cpt '
        f'-p {td}/topol.top -o {n1}/npt_restrained.tpr -maxwarn 2',
        log_file=n1 / 'grompp_npt1.log')

    gmx(f'mdrun -deffnm {n1}/npt_restrained -ntmpi 1 -ntomp 8 '
        f'-gpu_id 0 -nb gpu -pme gpu -bonded gpu -pin on',
        log_file=n1 / 'mdrun_npt1.log')

    if not (n1 / 'npt_restrained.gro').exists():
        raise RuntimeError(f'NPT1 failed for {poi}/{pose_id}')

    log(f'  {poi}/{pose_id}: NPT1 done')


def do_npt2(poi, pose_id):
    """500 ps NPT, protein restrained, PROTAC released."""
    td = topol_dir(poi, pose_id)
    n1 = npt1_dir(poi, pose_id)
    n2 = npt2_dir(poi, pose_id)
    n2.mkdir(parents=True, exist_ok=True)

    if (n2 / 'npt.gro').exists() and (n2 / 'npt.cpt').exists():
        log(f'  {poi}/{pose_id}: NPT2 already done — skipping')
        return

    gmx(f'grompp -f {MDP_DIR}/npt.mdp '
        f'-c {n1}/npt_restrained.gro -r {n1}/npt_restrained.gro '
        f'-t {n1}/npt_restrained.cpt '
        f'-p {td}/topol.top -o {n2}/npt.tpr -maxwarn 2',
        log_file=n2 / 'grompp_npt2.log')

    gmx(f'mdrun -deffnm {n2}/npt -ntmpi 1 -ntomp 8 '
        f'-gpu_id 0 -nb gpu -pme gpu -bonded gpu -pin on',
        log_file=n2 / 'mdrun_npt2.log')

    if not (n2 / 'npt.gro').exists():
        raise RuntimeError(f'NPT2 failed for {poi}/{pose_id}')

    log(f'  {poi}/{pose_id}: NPT2 done')


def _ensure_posre_sjf8240(td):
    """Generate posre_SJF8240.itp and patch SJF8240.itp — idempotent."""
    posre_itp = td / 'posre_SJF8240.itp'
    sjf_itp   = td / 'SJF8240.itp'

    if not posre_itp.exists():
        gmx(f'genrestr -f {td}/SJF8240.gro -o {posre_itp} -fc 1000 1000 1000',
            stdin='0\n', log_file=td / 'genrestr.log')

    if 'POSRES_LIG' not in sjf_itp.read_text():
        with open(sjf_itp, 'a') as f:
            f.write(textwrap.dedent("""
                ; PROTAC position restraints — active when -DPOSRES_LIG is defined
                #ifdef POSRES_LIG
                #include "posre_SJF8240.itp"
                #endif
            """))


# ── Production MD ─────────────────────────────────────────────────────────────
def do_production(poi, pose_id):
    """
    20 ns production MD. Resumes from checkpoint if crashed mid-run.
    Uses prod_20ns.mdp (nsteps=10000000).
    """
    td = topol_dir(poi, pose_id)
    n2 = npt2_dir(poi, pose_id)
    pd_ = prod_dir(poi, pose_id)
    pd_.mkdir(parents=True, exist_ok=True)

    tpr = pd_ / 'prod.tpr'
    cpt = pd_ / 'prod.cpt'
    gro = pd_ / 'prod.gro'
    xtc = pd_ / 'prod.xtc'

    if is_prod_done(poi, pose_id):
        log(f'  {poi}/{pose_id}: production already done — skipping')
        return

    # grompp (only if tpr doesn't exist — safe to rerun if no cpt)
    if not tpr.exists():
        gmx(f'grompp -f {MDP_DIR}/prod_20ns.mdp -c {n2}/npt.gro -t {n2}/npt.cpt '
            f'-p {td}/topol.top -o {tpr} -maxwarn 1',
            log_file=pd_ / 'grompp_prod.log')

    if cpt.exists() and not gro.exists():
        # Crashed — resume from checkpoint
        log(f'  {poi}/{pose_id}: resuming production from checkpoint')
        gmx(f'mdrun -s {tpr} -deffnm {pd_}/prod -cpi {cpt} -append '
            f'-ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu -pme gpu -bonded gpu -update gpu -pin on',
            log_file=pd_ / 'mdrun_prod.log', check=False)
    else:
        log(f'  {poi}/{pose_id}: starting 20 ns production')
        gmx(f'mdrun -s {tpr} -deffnm {pd_}/prod '
            f'-ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu -pme gpu -bonded gpu -update gpu -pin on',
            log_file=pd_ / 'mdrun_prod.log', check=False)

    if not is_prod_done(poi, pose_id):
        raise RuntimeError(
            f'Production did not complete for {poi}/{pose_id} '
            f'(prod.gro missing). Check {pd_}/mdrun_prod.log'
        )

    log(f'  {poi}/{pose_id}: production done')


# ── Analysis ──────────────────────────────────────────────────────────────────
def do_analysis(poi, pose_id):
    """
    GROMACS analysis + metrics.csv + dashboard.png.
    Imports helper functions from run_analysis.py — uses their exact logic
    with explicit Path arguments (functions internally call .exists() on args).
    """
    if str(ANA_SCRIPTS) not in sys.path:
        sys.path.insert(0, str(ANA_SCRIPTS))
    import importlib.util as _ilu
    _spec = _ilu.spec_from_file_location('run_analysis', ANA_SCRIPTS / '02_run_analysis.py')
    ra = _ilu.module_from_spec(_spec)
    _spec.loader.exec_module(ra)

    td   = topol_dir(poi, pose_id)
    pd_  = prod_dir(poi, pose_id)
    ad   = ana_dir(poi, pose_id)
    ad.mkdir(parents=True, exist_ok=True)

    # All kept as Path objects — run_analysis functions call .exists() on them
    gro         = td / 'complex_ions.gro'
    protein_gro = td / 'complex_protein.gro'
    tpr         = pd_ / 'prod.tpr'
    xtc_in      = pd_ / 'prod.xtc'

    for f in [gro, protein_gro, tpr, xtc_in]:
        if not f.exists():
            raise FileNotFoundError(f'Analysis input missing: {f}')

    # Build index — returns Path to index.ndx
    ndx = ra.build_index(gro, tpr, protein_gro, ad)

    # Preprocess trajectory — returns (xtc_nojump, xtc_stripped) as Path objects
    xtc_nojump, _xtc_stripped = ra.preprocess_traj(tpr, xtc_in, ndx, ad)

    # Compute distances/RMSD — pass Path objects (functions check .exists() to skip)
    xvg_dv = ad / 'dist_protac_vhl.xvg'
    xvg_dk = ad / 'dist_protac_kinase.xvg'
    xvg_rc = ad / 'rmsd_complex.xvg'
    xvg_rp = ad / 'rmsd_protac.xvg'

    ra.run_mindist(tpr, xtc_nojump, ndx, 'Other', 'VHL',          xvg_dv)
    ra.run_mindist(tpr, xtc_nojump, ndx, 'Other', 'Kinase',       xvg_dk)
    ra.run_rms(tpr, xtc_nojump, ndx, 'VHL_Kinase_CA', 'VHL_Kinase_CA', xvg_rc)
    ra.run_rms(tpr, xtc_nojump, ndx, 'VHL_Kinase_CA', 'Other',          xvg_rp)

    # Parse XVGs — distances in nm; convert to Å for all downstream calculations
    t_d, dv = ra.read_xvg(xvg_dv)
    _,   dk = ra.read_xvg(xvg_dk)
    t_r, rc = ra.read_xvg(xvg_rc)
    _,   rp = ra.read_xvg(xvg_rp)

    dv_A = dv[:, 0] * 10
    dk_A = dk[:, 0] * 10
    rc_A = rc[:, 0] * 10
    rp_A = rp[:, 0] * 10

    n = min(len(dv_A), len(dk_A))
    bridging_mask = (dv_A[:n] < 5.0) & (dk_A[:n] < 5.0)
    bridging_pct = float(bridging_mask.mean() * 100)

    label = f'{poi} / {pose_id}'

    # Write metrics.csv
    stats_header = (
        f'# {label} — 20 ns random screen\n'
        f'# PROTAC-VHL:    {dv_A.mean():.2f} ± {dv_A.std():.2f} Å'
        f'  [min {dv_A.min():.2f}, max {dv_A.max():.2f}]\n'
        f'# PROTAC-kinase: {dk_A.mean():.2f} ± {dk_A.std():.2f} Å'
        f'  [min {dk_A.min():.2f}, max {dk_A.max():.2f}]\n'
        f'# Bridging (PROTAC-VHL<5Å AND PROTAC-kinase<5Å): {bridging_pct:.1f}% of frames\n'
        f'time_ns,d_protac_vhl_A,d_protac_kinase_A,rmsd_complex_A,rmsd_protac_A\n'
    )
    rows = np.column_stack([t_d[:n], dv_A[:n], dk_A[:n], rc_A[:n], rp_A[:n]])
    with open(ad / 'metrics.csv', 'w') as f:
        f.write(stats_header)
        np.savetxt(f, rows, delimiter=',', fmt='%.4f')

    # Dashboard plot
    _plot_dashboard_pose(t_d, dv_A, dk_A, t_r, rc_A, rp_A, bridging_pct, label, ad)

    log(f'  {poi}/{pose_id}: analysis done — bridging {bridging_pct:.1f}%')


def _plot_dashboard_pose(t_d, dv_A, dk_A, t_r, rc_A, rp_A,
                          bridging_pct, label, out_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    fig.suptitle(f'{label}\n20 ns random screen | AMBER99SB-ILDN + GAFF2 + TIP3P',
                 fontsize=12, fontweight='bold')

    configs = [
        (axes[0,0], t_d, dv_A,  'PROTAC → VHL (Å)',        'steelblue', True),
        (axes[0,1], t_d, dk_A,  'PROTAC → kinase (Å)',     'coral',     True),
        (axes[1,0], t_r, rc_A,  'VHL+kinase Cα RMSD (Å)',  'seagreen',  False),
        (axes[1,1], t_r, rp_A,  'PROTAC RMSD (Å)',         'mediumpurple', False),
    ]
    for ax, t, y, ylabel, col, thresholds in configs:
        ax.plot(t, y, lw=0.6, color=col, alpha=0.5)
        ax.set(xlabel='Time (ns)', ylabel=ylabel)
        ax.set_ylim(bottom=0)
        if thresholds:
            ax.axhline(5.0,  color='green',  ls='--', lw=1.2, alpha=0.7, label='PASS (5 Å)')
            ax.axhline(10.0, color='orange', ls='--', lw=1.2, alpha=0.7, label='MARGINAL')
            ax.legend(fontsize=8)

    txt = (f'Bridging (<5Å both): {bridging_pct:.1f}%  |  '
           f'PROTAC–VHL: {dv_A.mean():.1f}±{dv_A.std():.1f} Å  |  '
           f'PROTAC–kinase: {dk_A.mean():.1f}±{dk_A.std():.1f} Å')
    fig.text(0.5, -0.01, txt, ha='center', fontsize=9, style='italic',
             bbox=dict(boxstyle='round,pad=0.4', fc='lightyellow', alpha=0.9))

    fig.savefig(str(out_dir / 'dashboard.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)


# ── Main run loop ─────────────────────────────────────────────────────────────
STEP_FNS = {
    'PENDING':       ('STRUCTURE_DONE',  lambda p, r: do_structure(p, r, None)),  # src injected below
    'STRUCTURE_DONE':('TOPOLOGY_DONE',   do_topology),
    'TOPOLOGY_DONE': ('EM_DONE',         do_em),
    'EM_DONE':       ('NVT_DONE',        do_nvt),
    'NVT_DONE':      ('NPT1_DONE',       do_npt1),
    'NPT1_DONE':     ('NPT2_DONE',       do_npt2),
    'NPT2_DONE':     ('PROD_DONE',       do_production),
    'PROD_DONE':     ('ANALYSIS_DONE',   do_analysis),
}


def run_pose(df, poi, pose_id, source_pdb):
    """
    Advance a single pose through all remaining stages.
    Reads current stage from manifest, reconciles with disk state, then steps forward.
    Returns updated df.
    """
    row = get_pose_row(df, poi, pose_id)
    stage = row['stage']
    attempts = int(row['attempts'])

    # Reconcile manifest stage with disk (handles restarts where manifest was not updated)
    disk_stage = detect_stage_from_disk(poi, pose_id)
    if STAGE_ORDER.get(disk_stage, 0) > STAGE_ORDER.get(stage, 0):
        log(f'  {poi}/{pose_id}: disk stage ({disk_stage}) ahead of manifest ({stage}) — syncing')
        stage = disk_stage
        df = update_pose(df, poi, pose_id, stage=stage)

    if stage == 'ANALYSIS_DONE':
        log(f'  {poi}/{pose_id}: already complete')
        return df

    if stage == 'FAILED' and attempts >= MAX_RETRIES:
        log(f'  {poi}/{pose_id}: FAILED (max retries reached)')
        return df

    log(f'\n{"="*60}')
    log(f'  Processing {poi}/{pose_id}  (stage={stage}, attempt {attempts+1})')
    log(f'{"="*60}')

    df = update_pose(df, poi, pose_id, increment_attempts=True)

    # Step through all remaining stages
    current = stage
    while current != 'ANALYSIS_DONE':
        if current not in STEP_FNS:
            log(f'  ERROR: unknown stage {current}')
            break
        next_stage, fn = STEP_FNS[current]

        try:
            if current == 'PENDING':
                do_structure(poi, pose_id, source_pdb)
            else:
                fn(poi, pose_id)
            current = next_stage
            df = update_pose(df, poi, pose_id, stage=current, error='')
        except Exception as e:
            err_msg = str(e)[:300]
            log(f'  FAILED at {current}: {err_msg}')
            df = update_pose(df, poi, pose_id, stage='FAILED', error=err_msg)
            return df

    return df


def run_campaign(filter_poi=None, filter_pose=None):
    """Execute all PENDING/FAILED-retryable poses sequentially."""
    if not MANIFEST.exists():
        log('ERROR: manifest not found. Run --plan first.')
        sys.exit(1)

    df = read_manifest()
    total = len(df)
    done_count = (df['stage'] == 'ANALYSIS_DONE').sum()
    log(f'Campaign: {total} poses total, {done_count} already done')

    # Build ordered list of poses to process
    to_run = []
    for _, row in df.iterrows():
        poi     = row['poi']
        pose_id = row['pose_id']
        stage   = row['stage']
        attempts = int(row['attempts'])

        if filter_poi and poi != filter_poi:
            continue
        if filter_pose and pose_id != filter_pose:
            continue
        if stage == 'ANALYSIS_DONE':
            continue
        if stage == 'FAILED' and attempts >= MAX_RETRIES:
            continue
        to_run.append((poi, pose_id, row['source_pdb']))

    log(f'{len(to_run)} poses to run (or resume)')

    for i, (poi, pose_id, source_pdb) in enumerate(to_run):
        log(f'\n[{i+1}/{len(to_run)}] {poi}/{pose_id} ({source_pdb})')
        df = read_manifest()  # re-read to get latest state
        df = run_pose(df, poi, pose_id, source_pdb)

    # Final summary
    df = read_manifest()
    done  = (df['stage'] == 'ANALYSIS_DONE').sum()
    failed = (df['stage'] == 'FAILED').sum()
    log(f'\nCampaign batch complete: {done}/{len(df)} done, {failed} failed')


# ── Status display ────────────────────────────────────────────────────────────
def print_status(df=None):
    if df is None:
        df = read_manifest()

    if df.empty:
        print('No poses in manifest yet. Run --plan first.')
        return

    print(f'\n{"POI":<14} {"Pose":<10} {"Source PDB":<25} {"Stage":<18} {"Att":>3} {"Updated":<17} {"Error"}')
    print('-' * 110)

    stage_counts = {}
    for _, row in df.sort_values(['poi','pose_id']).iterrows():
        stage = row['stage']
        stage_counts[stage] = stage_counts.get(stage, 0) + 1
        err = row['error'][:40] if row['error'] else ''
        print(f"{row['poi']:<14} {row['pose_id']:<10} {row['source_pdb']:<25} "
              f"{stage:<18} {row['attempts']:>3} {row['last_updated']:<17} {err}")

    print('-' * 110)
    print('Summary:', '  '.join(f'{s}: {n}' for s, n in sorted(stage_counts.items())))
    total = len(df)
    done  = stage_counts.get('ANALYSIS_DONE', 0)
    print(f'Progress: {done}/{total} ({done/total*100:.0f}%)\n')


# ── Entry point ───────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--plan',   action='store_true',
                       help='Sample poses and write manifest (dry run)')
    group.add_argument('--run',    action='store_true',
                       help='Execute / resume campaign')
    group.add_argument('--status', action='store_true',
                       help='Show current progress')

    parser.add_argument('--poi',  default=None,
                        help='Restrict to one POI (use with --run)')
    parser.add_argument('--pose', default=None,
                        help='Restrict to one pose ID (use with --run and --poi)')
    args = parser.parse_args()

    if args.plan:
        plan_campaign()
    elif args.run:
        run_campaign(filter_poi=args.poi, filter_pose=args.pose)
    elif args.status:
        print_status()


if __name__ == '__main__':
    main()
