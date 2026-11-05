#!/usr/bin/env python3
import os, sys, shutil, subprocess
from pathlib import Path

# Works on any machine — resolves relative to this script's location
_HERE  = Path(__file__).resolve().parent          # scripts/rosetta_utils/
BASE   = str(_HERE.parent.parent)                 # project root
SCRIPT = str(_HERE / 'resume_pipeline.py')

POIS = [
    ('high_affinity_degraders',    'MET'),
    ('high_affinity_degraders',    'DDR2'),
    ('high_affinity_degraders',    'RIPK2'),
    ('high_affinity_degraders',    'EPHB2'),
    ('high_affinity_degraders',    'MAPK14'),
    ('high_affinity_no_degradation', 'AXL'),
    ('high_affinity_no_degradation', 'ABL1'),
    ('high_affinity_no_degradation', 'EPHA2'),
]

env = os.environ.copy()
env['N_WORKERS_RESUME'] = '35'

for subdir, poi in POIS:
    run_dir = f'{BASE}/pRosettaC/runs/{subdir}/{poi}'
    print(f'=== Starting {poi} ===', flush=True)
    os.chdir(run_dir)
    r = subprocess.run([sys.executable, SCRIPT], env=env)
    print(f'=== Done {poi} (rc={r.returncode}) ===', flush=True)

print('All 8 POIs complete.')
