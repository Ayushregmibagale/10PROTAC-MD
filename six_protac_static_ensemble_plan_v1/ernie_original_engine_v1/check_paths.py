"""Read-only check of the original PatchDock + Rosetta installation, no docking."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def inspect(config):
    checks = []
    def check(name, path, executable=False):
        p = Path(path)
        good = p.is_file() and (not executable or os.access(p, os.X_OK))
        checks.append(dict(name=name, path=str(p), ok=good,
                           sha256=sha(p) if good else None))
        return good
    source = Path(config['source'])
    patch = Path(config['patchdock'])
    rosetta = Path(config['rosetta'])
    for name in ('main.py', 'utils.py', 'rosetta.py', 'protac_lib.py',
                 'constraint_generation.py', 'clustering.py', 'Rosetta/docking.xml',
                 'Rosetta/docking.flags', 'Rosetta/relax.flags',
                 'cluster/LOCAL/LOCAL.py'):
        check('source:' + name, source / name)
    for name in ('patch_dock.Linux', 'buildParams.pl', 'transOutput.pl'):
        check('PatchDock:' + name, patch / name, True)
    check('RosettaScripts', rosetta / 'main/source/bin/rosetta_scripts.default.linuxgccrelease', True)
    check('clean_pdb.py', rosetta / 'tools/protein_tools/scripts/clean_pdb.py', True)
    check('molfile_to_params.py', rosetta / 'main/source/scripts/python/public/molfile_to_params.py')
    check('babel (original code spelling)', Path(config['openbabel']) / 'babel', True)
    check('python3', config['python'], True)
    effective_path = str(Path(config['python']).parent) + os.pathsep + os.environ.get('PATH', '')
    for command in ('python2.7', 'perl', 'bash', 'sed', 'grep', 'cat', 'mv'):
        path = shutil.which(command, path=effective_path)
        checks.append(dict(name='PATH:' + command, path=path, ok=path is not None))
    utils = source / 'utils.py'
    if utils.is_file():
        text = utils.read_text()
        checks.append(dict(name='source calls PatchDock executable',
                           ok='/patch_dock.Linux' in text and '/buildParams.pl' in text))
    main = source / 'main.py'
    if main.is_file():
        text = main.read_text()
        checks.append(dict(name='entrypoint has no PyRosetta adapter import',
                           ok='pyrosetta_adapter' not in text and 'apply_patches' not in text))
    python = Path(config['python'])
    if python.is_file():
        probe = subprocess.run([str(python), '-c',
            'import rdkit,numpy,scipy,sklearn,pymol; print("imports OK")'],
            capture_output=True, text=True, timeout=45,
            env={**os.environ, 'PYTHONDONTWRITEBYTECODE': '1'})
        checks.append(dict(name='original Python imports', ok=probe.returncode == 0,
                           stdout=probe.stdout[-1000:], stderr=probe.stderr[-2000:]))
    return dict(status='PATHS_PASS' if all(c['ok'] for c in checks) else 'PATHS_BLOCKED',
                checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                hostname=socket.gethostname(), checks=checks,
                scheduler_commands={x:shutil.which(x) for x in ('sbatch','qsub','srun')},
                scope='Original engine paths/imports only; no molecular jobs or reference-recovery gate',
                docking_launched=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    args = parser.parse_args()
    report = inspect(json.loads(args.config.read_text()))
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if report['status'] == 'PATHS_PASS' else 2)
