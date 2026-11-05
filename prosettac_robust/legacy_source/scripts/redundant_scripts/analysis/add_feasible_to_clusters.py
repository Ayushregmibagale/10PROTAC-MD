#!/usr/bin/env python3
"""
Add PASS/MARGINAL poses (rsynced from ernie) to the correct Results/clusterN/ directories.

Cluster assignment heuristic: combined_X_Y_0001.pdb shares the same PatchDock global
solution X as other poses already in cluster directories. We find which cluster contains
any combined_X_*_0001.pdb and assign the feasible pose to that cluster.

Usage (run on HOME):
    python add_feasible_to_clusters.py DDR2 EPHB2
"""

import sys
import shutil
from pathlib import Path

RUNS_ROOT = Path(__file__).resolve().parents[2] / 'pRosettaC' / 'runs'

POI_CONFIGS = {
    'DDR2':   {'category': 'high_affinity_degraders'},
    'EPHB2':  {'category': 'high_affinity_degraders'},
    'RIPK2':  {'category': 'high_affinity_degraders'},
    'MAPK14': {'category': 'high_affinity_degraders'},
}


def build_global_to_cluster_map(results_dir):
    """Map global docking ID (X in combined_X_Y_0001.pdb) → cluster directory."""
    global_map = {}
    for cluster_dir in sorted(results_dir.glob('cluster*')):
        if not cluster_dir.is_dir():
            continue
        for pdb in cluster_dir.glob('combined_*_0001.pdb'):
            parts = pdb.stem.split('_')
            global_id = parts[1]
            if global_id not in global_map:
                global_map[global_id] = cluster_dir
    return global_map


def run_poi(poi_name, feasible_file):
    cfg = POI_CONFIGS[poi_name]
    run_dir = RUNS_ROOT / cfg['category'] / poi_name
    results_dir = run_dir / 'Results'
    patchdock_dir = run_dir / 'Patchdock_Results'

    if not feasible_file.exists():
        print(f'{poi_name}: feasible file not found: {feasible_file}')
        return

    lines = [l.strip() for l in feasible_file.read_text().splitlines() if l.strip()]
    if not lines:
        print(f'{poi_name}: feasible file is empty')
        return

    global_map = build_global_to_cluster_map(results_dir)
    print(f'{poi_name}: {len(global_map)} global IDs mapped to {len(set(global_map.values()))} clusters')

    added = 0
    unmapped = 0
    already_present = 0

    for line in lines:
        parts = line.split('\t')
        if len(parts) != 2:
            continue
        verdict, rel_path = parts
        src = patchdock_dir / Path(rel_path).name

        if not src.exists():
            print(f'  WARNING: {src.name} not found in Patchdock_Results — skipped')
            continue

        stem_parts = src.stem.split('_')
        global_id = stem_parts[1]
        cluster_dir = global_map.get(global_id)

        if cluster_dir is None:
            unmapped += 1
            continue

        dst = cluster_dir / src.name
        existed = dst.exists()
        shutil.copy2(src, dst)
        if existed:
            already_present += 1
        else:
            added += 1

    print(f'{poi_name}: added {added} poses | already present {already_present} | unmapped {unmapped}')


def main():
    pois = sys.argv[1:] if len(sys.argv) > 1 else list(POI_CONFIGS)
    home = Path.home()

    for poi in pois:
        if poi not in POI_CONFIGS:
            print(f'Unknown POI: {poi}')
            continue
        feasible_file = home / f'{poi.lower()}_feasible.txt'
        run_poi(poi, feasible_file)


if __name__ == '__main__':
    main()
