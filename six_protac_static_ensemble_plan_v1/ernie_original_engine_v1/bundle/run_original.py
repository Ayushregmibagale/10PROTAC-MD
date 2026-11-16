"""Run original main.py; only bound the LOCAL scheduler and propagate its errors.

No replacement docking, restraint, linker or score implementation is imported.
"""
import argparse
import concurrent.futures
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--workers', type=int, required=True)
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    if args.workers < 1:
        raise ValueError('workers must be positive')
    source = Path(cfg['source']).resolve()
    os.environ.update(PATCHDOCK=cfg['patchdock'], ROSETTA_FOL=cfg['rosetta'],
                      OB=cfg['openbabel'], SCRIPTS_FOL=str(source) + '/',
                      OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    sys.path.insert(0, str(source))
    cluster = importlib.import_module('cluster')
    from cluster.LOCAL.LOCAL import LOCAL

    class BoundedLocal(LOCAL):
        def __init__(self):
            super().__init__(n_workers=args.workers)

        def runBatchCommands(self, commands, mem, batch_size=12):
            def run(command):
                result = subprocess.run(command, shell=True, executable='/bin/bash')
                if result.returncode:
                    raise RuntimeError(f'Original subprocess failed ({result.returncode}): {command}')
            with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
                # Iteration forces every exception to reach the original main process.
                list(executor.map(run, commands))
            return ['local_done']

    original_get = cluster.getCluster
    cluster.getCluster = lambda name: BoundedLocal() if name == 'LOCAL' else original_get(name)
    entry = importlib.import_module('main')
    if Path(entry.__file__).resolve() != source / 'main.py':
        raise RuntimeError('Wrong main module imported')
    entry.main('main.py', ['Protac_params.txt'])


if __name__ == '__main__':
    main()
