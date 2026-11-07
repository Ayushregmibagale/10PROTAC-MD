#!/usr/bin/env python3
"""
Screen pRosettaC pose ensembles for PROTAC bridging geometry across all POIs.

For each POI, computes min PROTAC–VHL and PROTAC–kinase distances per pose,
then reports: n_bridging / n_total (both < THRESHOLD) as the bridging fraction.

This fraction = P(productive | ternary complex) — key ML feature.

Output:
  analysis/<POI>/bridging_screen.tsv   — per-pose distances
  analysis/bridging_summary.tsv        — one row per POI (the ML feature table)
"""
import numpy as np
from pathlib import Path
from scipy.spatial.distance import cdist
from multiprocessing import Pool
import sys

ROSETTA_BASE = Path.home() / "PROTAC degradation Prediction/pRosettaC/runs"
ANALYSIS_DIR = Path.home() / "protac-md/analysis"
THRESHOLD = 5.0  # Å — both VHL and kinase must be < this
N_WORKERS = 19   # ~80% of 24 cores

SUMMARY_TSV = ANALYSIS_DIR / 'bridging_summary.tsv'

# POIs to force re-run even if cached results exist (chain-D reruns invalidate old cache)
RERUN_POIS = {'TNIK', 'MAPK14', 'MAP4K5', 'SLK'}

# Config per POI:
#   base       : subdirectory under ROSETTA_BASE
#   src_subdir : directory under the POI folder containing PDBs (default: Patchdock_Results)
#   glob       : filename pattern for combined PDBs (supports subdirs, e.g. cluster*/combined_*.pdb)
#   vhl_max    : chain A residues <= vhl_max are VHL
#   kin_chain  : chain containing kinase atoms
#   kin_min    : residues >= kin_min in kin_chain are kinase
#   category   : experimental label
POIS = {
    # ── High affinity degraders ──────────────────────────────────────────────
    'DDR2':        {'base': 'high_affinity_degraders',        'glob': 'combined_*.pdb', 'vhl_max': 103, 'kin_chain': 'A', 'kin_min': 553,  'category': 'high_degrader'},
    'EPHB2':       {'base': 'high_affinity_degraders',        'glob': 'combined_*.pdb', 'vhl_max': 103, 'kin_chain': 'A', 'kin_min': 615,  'category': 'high_degrader'},
    'MAPK14_full': {'base': 'high_affinity_degraders',        'glob': 'combined_*.pdb', 'vhl_max': 103, 'kin_chain': 'A', 'kin_min': 200,  'category': 'high_degrader'},
    'MET':         {'base': 'high_affinity_degraders',        'glob': 'combined_*.pdb', 'vhl_max': 103, 'kin_chain': 'A', 'kin_min': 1053, 'category': 'high_degrader'},
    'RIPK2_full':  {'base': 'high_affinity_degraders',        'glob': 'combined_*.pdb', 'vhl_max': 103, 'kin_chain': 'A', 'kin_min': 200,  'category': 'high_degrader'},
    # ── High affinity non-degraders ──────────────────────────────────────────
    'ABL1':        {'base': 'high_affinity_no_degradation',   'glob': 'combined_*.pdb', 'vhl_max': 103, 'kin_chain': 'A', 'kin_min': 235,  'category': 'high_non_degrader'},
    'AXL':         {'base': 'high_affinity_no_degradation',   'glob': 'combined_*.pdb', 'vhl_max': 103, 'kin_chain': 'A', 'kin_min': 521,  'category': 'high_non_degrader'},
    'EPHA2':       {'base': 'high_affinity_no_degradation',   'glob': 'combined_*.pdb', 'vhl_max': 103, 'kin_chain': 'A', 'kin_min': 600,  'category': 'high_non_degrader'},
    # chain-D reruns — kinase on chain D, sourced from Results/cluster*/ directories
    'MAP4K5':      {'base': 'high_affinity_no_degradation',   'src_subdir': 'Results', 'glob': 'cluster*/combined_*.pdb', 'vhl_max': 103, 'kin_chain': 'D', 'kin_min': 0, 'category': 'high_non_degrader'},
    'SLK':         {'base': 'high_affinity_no_degradation',   'src_subdir': 'Results', 'glob': 'cluster*/combined_*.pdb', 'vhl_max': 103, 'kin_chain': 'D', 'kin_min': 0, 'category': 'high_non_degrader'},
    # ── Low affinity degraders ───────────────────────────────────────────────
    'CDK17':       {'base': 'low_affinity_degraders',         'glob': 'combined_*.pdb', 'vhl_max': 103, 'kin_chain': 'D', 'kin_min': 1,    'category': 'low_degrader'},
    # chain-D reruns — kinase on chain D, sourced from Results/cluster*/ directories
    'MAPK14':      {'base': 'low_affinity_degraders',         'src_subdir': 'Results', 'glob': 'cluster*/combined_*.pdb', 'vhl_max': 103, 'kin_chain': 'D', 'kin_min': 0, 'category': 'low_degrader'},
    'TNIK':        {'base': 'low_affinity_degraders',         'src_subdir': 'Results', 'glob': 'cluster*/combined_*.pdb', 'vhl_max': 103, 'kin_chain': 'D', 'kin_min': 0, 'category': 'low_degrader'},
    # ── Low affinity non-degraders (new) ────────────────────────────────────
    'CDK2':        {'base': 'low_affinity_no_degradation',    'glob': 'combined_*.pdb', 'vhl_max': 103, 'kin_chain': 'D', 'kin_min': 0,    'category': 'low_non_degrader'},
    'CDK4':        {'base': 'low_affinity_no_degradation',    'glob': 'combined_*.pdb', 'vhl_max': 103, 'kin_chain': 'D', 'kin_min': 0,    'category': 'low_non_degrader'},
    'CDK6':        {'base': 'low_affinity_no_degradation',    'glob': 'combined_*.pdb', 'vhl_max': 103, 'kin_chain': 'D', 'kin_min': 0,    'category': 'low_non_degrader'},
    'DAPK1':       {'base': 'low_affinity_no_degradation',    'glob': 'combined_*.pdb', 'vhl_max': 103, 'kin_chain': 'D', 'kin_min': 0,    'category': 'low_non_degrader'},
    'IRAK1':       {'base': 'low_affinity_no_degradation',    'glob': 'combined_*.pdb', 'vhl_max': 103, 'kin_chain': 'D', 'kin_min': 0,    'category': 'low_non_degrader'},
    'TLK2':        {'base': 'low_affinity_no_degradation',    'glob': 'combined_*.pdb', 'vhl_max': 103, 'kin_chain': 'D', 'kin_min': 0,    'category': 'low_non_degrader'},
    # ── Kinome screen (Phase 2 full pRosettaC runs, renumbered res 200+) ────────
    'sik3_full':    {'abs_dir': '/mnt/d/kinome_results/sik3_full',    'glob': 'combined_*.pdb', 'vhl_max': 103, 'kin_chain': 'A', 'kin_min': 200, 'category': 'kinome_screen'},
    'rps6kb2_full': {'abs_dir': '/mnt/d/kinome_results/rps6kb2_full', 'glob': 'combined_*.pdb', 'vhl_max': 103, 'kin_chain': 'A', 'kin_min': 200, 'category': 'kinome_screen'},
}


def _process_pdb(args):
    """Worker: parse one PDB, return (stem, d_vhl, d_kin) or None."""
    pdb, vhl_max, kin_chain, kin_min = args
    vhl, kin, ptc = [], [], []
    try:
        with open(pdb) as fh:
            for line in fh:
                rec = line[:6].strip()
                if rec not in ('ATOM', 'HETATM'):
                    continue
                chain  = line[21]
                resnum = int(line[22:26])
                aname  = line[12:16].strip()
                coord  = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                if chain == 'A' and resnum <= vhl_max:
                    vhl.append(coord)
                if chain == kin_chain and resnum >= kin_min:
                    kin.append(coord)
                if chain == 'X' and not aname.startswith('H'):
                    ptc.append(coord)
    except Exception:
        return None
    if not vhl or not kin or not ptc:
        return None
    vhl = np.array(vhl); kin = np.array(kin); ptc = np.array(ptc)
    d_vhl = float(cdist(ptc, vhl).min())
    d_kin  = float(cdist(ptc, kin).min())
    return (Path(pdb).stem, d_vhl, d_kin)


# Load existing summary rows to preserve already-processed POIs
existing_rows = {}
if SUMMARY_TSV.exists():
    with open(SUMMARY_TSV) as f:
        f.readline()  # header
        for line in f:
            parts = line.strip().split('\t')
            if parts:
                existing_rows[parts[0]] = line.strip()

summary_rows = []

for poi, cfg in POIS.items():
    out_dir = ANALYSIS_DIR / poi
    out_tsv = out_dir / 'bridging_screen.tsv'

    if poi in existing_rows and out_tsv.exists() and poi not in RERUN_POIS:
        print(f"[{poi}] already processed — skipping")
        summary_rows.append(existing_rows[poi])
        continue

    src_subdir = cfg.get('src_subdir', 'Patchdock_Results')
    if 'abs_dir' in cfg:
        pd_dir = Path(cfg['abs_dir']) / src_subdir
    else:
        pd_dir = ROSETTA_BASE / cfg['base'] / poi / src_subdir
    pdbs   = sorted(pd_dir.glob(cfg['glob']))

    if not pdbs:
        print(f"[{poi}] WARNING: no PDBs found at {pd_dir} (glob: {cfg['glob']})")
        continue

    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"  {poi}: screening {len(pdbs)} poses  "
          f"(VHL: chain A ≤{cfg['vhl_max']}, KIN: chain {cfg['kin_chain']} ≥{cfg['kin_min']})  "
          f"[{N_WORKERS} workers]")
    print(f"{'='*60}", flush=True)

    args_iter = [(str(p), cfg['vhl_max'], cfg['kin_chain'], cfg['kin_min']) for p in pdbs]
    results = []
    done = 0
    with Pool(N_WORKERS) as pool:
        for res in pool.imap_unordered(_process_pdb, args_iter, chunksize=50):
            done += 1
            if res is not None:
                results.append(res)
            if done % 2000 == 0:
                n_bridge_so_far = sum(1 for _, v, k in results if v < THRESHOLD and k < THRESHOLD)
                print(f"  {done}/{len(pdbs)}  ({n_bridge_so_far} bridging so far)", flush=True)

    if not results:
        print(f"  [WARN] no valid poses parsed for {poi}")
        continue

    results.sort(key=lambda r: r[1] + r[2])
    bridging  = [(n, v, k) for n, v, k in results if v < THRESHOLD and k < THRESHOLD]
    n_total   = len(results)
    n_bridge  = len(bridging)
    frac      = n_bridge / n_total * 100

    print(f"  Bridging (<{THRESHOLD}Å both): {n_bridge}/{n_total} = {frac:.3f}%")

    if bridging:
        print(f"  Top 5 bridging poses:")
        print(f"  {'Pose':<35} {'VHL(Å)':>8} {'KIN(Å)':>8}")
        for name, v, k in bridging[:5]:
            print(f"  {name:<35} {v:>8.2f} {k:>8.2f}")

    with open(out_tsv, 'w') as f:
        f.write("pose\td_vhl\td_kin\tbridging\n")
        for name, v, k in results:
            f.write(f"{name}\t{v:.3f}\t{k:.3f}\t{'1' if v<THRESHOLD and k<THRESHOLD else '0'}\n")

    best_vhl = round(bridging[0][1], 3) if bridging else 'NA'
    best_kin = round(bridging[0][2], 3) if bridging else 'NA'
    summary_rows.append(
        f"{poi}\t{cfg['category']}\t{n_total}\t{n_bridge}\t{round(frac, 4)}\t{best_vhl}\t{best_kin}"
    )

# ── Summary table ─────────────────────────────────────────────────────────────
with open(SUMMARY_TSV, 'w') as f:
    f.write("poi\tcategory\tn_total\tn_bridging\tbridge_pct\tbest_vhl_A\tbest_kin_A\n")
    for row in summary_rows:
        f.write(row + '\n')

print(f"\n{'='*60}")
print(f"SUMMARY — bridging fraction per POI (threshold {THRESHOLD} Å)")
print(f"{'='*60}")
print(f"{'POI':<16} {'n_bridging':>12} {'bridge_pct':>12}")
print(f"{'-'*50}")
for row in summary_rows:
    parts = row.split('\t')
    print(f"{parts[0]:<16} {parts[3]:>12} {parts[4]:>11}%")
print(f"\nFull summary: {SUMMARY_TSV}")
