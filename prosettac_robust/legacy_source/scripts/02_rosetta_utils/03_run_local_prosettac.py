#!/usr/bin/env python3
"""
Local pRosettaC runner — no cluster scheduler, no Rosetta binaries, no PatchDock.

Usage (run from the target's run directory):
  cd ~/PROTAC\ degradation\ Prediction/pRosettaC/runs/<category>/<target>
  python ~/PROTAC\ degradation\ Prediction/scripts/02_rosetta_utils/03_run_local_prosettac.py \\
      Protac_params.txt [--n_dock 200] [--n_workers 4]

What it does:
  1. Applies pyrosetta_adapter patches (PyRosetta + Python OpenBabel replace executables)
  2. Delegates to pRosettaC's main.py with ClusterName=LOCAL (runs locally via
     Python's ProcessPoolExecutor instead of PBS/SGE/SLURM)

Prerequisites (all already installed):
  - PyRosetta       (~/miniconda3, Python 3.13)
  - openbabel-wheel (Python OpenBabel)
  - RDKit, scipy, numpy
  - PRosettaC repo  (~/PRosettaC)
"""

import sys
import os
import argparse
from pathlib import Path

# ── paths ─────────────────────────────────────────────────────────────────────
_SCRIPT_DIR   = Path(__file__).parent
_PROJECT_ROOT = _SCRIPT_DIR.parent.parent
_PROSETTAC    = Path.home() / 'PRosettaC'

# ── inject PRosettaC into path so its modules are importable ──────────────────
if str(_PROSETTAC) not in sys.path:
    sys.path.insert(0, str(_PROSETTAC))

# ── apply monkey-patches BEFORE importing pRosettaC ──────────────────────────
sys.path.insert(0, str(_SCRIPT_DIR))
from pyrosetta_adapter import apply_patches
apply_patches()

# ── now import pRosettaC main ─────────────────────────────────────────────────
import main as prosettac_main


def parse_args():
    p = argparse.ArgumentParser(description='Local pRosettaC runner')
    p.add_argument('params_file', help='Protac_params.txt path')
    p.add_argument('--n_dock',    type=int, default=200,
                   help='Max global docking poses to generate (default 200)')
    p.add_argument('--n_workers', type=int, default=None,
                   help='Parallel workers for local docking (default: CPU count - 1)')
    return p.parse_args()


def main():
    args = parse_args()

    params_path = Path(args.params_file)
    if not params_path.exists():
        sys.exit(f'ERROR: {params_path} not found.\n'
                 f'Run this script from the target run directory.')

    # Patch LOCAL cluster worker count if requested
    if args.n_workers:
        from cluster.LOCAL.LOCAL import LOCAL as LocalCluster
        LocalCluster.__init__ = lambda self, _=None: setattr(self, 'n_workers', args.n_workers)

    # Optionally scale the number of global docking poses
    # pRosettaC uses Global=1000 (Full=True) or 500 (Full=False);
    # override via the adapter's global_dock if needed.
    # For local workstations 200 is a practical default.
    import pyrosetta_adapter as _adp
    _orig_global_dock = _adp.global_dock

    def _scaled_global_dock(structs, anchors, min_dist, max_dist,
                             num_results=1000, threshold=2.0):
        return _orig_global_dock(structs, anchors, min_dist, max_dist,
                                 num_results=min(num_results, args.n_dock),
                                 threshold=threshold)

    import utils
    utils.patchdock = _scaled_global_dock

    print(f'[run_local] Starting pRosettaC for: {params_path.resolve()}')
    print(f'[run_local] Max global poses: {args.n_dock}  |  Workers: {args.n_workers or "auto"}')

    prosettac_main.main(sys.argv[0], [str(params_path)])

    # pRosettaC's main.py builds constraint_generation commands with "../abs_path"
    # (broken for absolute heads/linker paths). Re-run constraint_gen + clustering
    # via resume_pipeline which uses correct absolute paths and shifts scores.
    import subprocess
    print('\n[run_local] Re-running constraint_gen + clustering via 06_resume_pipeline...')
    subprocess.run(
        [sys.executable, str(_SCRIPT_DIR / '06_resume_pipeline.py')],
        check=True,
        cwd=os.getcwd(),
    )


if __name__ == '__main__':
    main()
