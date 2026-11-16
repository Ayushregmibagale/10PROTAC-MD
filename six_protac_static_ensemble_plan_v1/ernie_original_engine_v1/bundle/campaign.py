"""Plan/run a task manifest with original PatchDock/Rosetta pRosettaC.

Default is a read-only plan. --execute starts tasks only after path checks.
Completed attempts are hash-verified and skipped; failed attempts are preserved.
"""
import argparse
import datetime
import fcntl
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time

from check_paths import inspect, sha


def write(path, data):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def validate_tasks(root, manifest):
    seen = set()
    if not manifest['tasks']:
        raise ValueError('Empty campaign manifest')
    for task in manifest['tasks']:
        key = task['id']
        if key in seen or not key or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in key):
            raise ValueError('Duplicate or unsafe task ID')
        seen.add(key)
        directory = (root / task['input']).resolve()
        directory.relative_to(root.resolve())
        for name, digest in task['files'].items():
            path = (directory / name).resolve()
            path.relative_to(directory)
            if sha(path) != digest:
                raise ValueError(f'Input changed: {key}/{name}')


def outputs(run):
    score = run / 'Patchdock_Results/score.sc'
    if not score.is_file() or not (run / 'Results').is_dir():
        raise RuntimeError('Original process ended without score.sc and Results; inspect log, not a completed success')
    lines = score.read_text().splitlines()
    values = []
    for line in lines:
        fields = line.split()
        if len(fields) > 2 and fields[0] == 'SCORE:':
            try:
                value = float(fields[1])
            except ValueError:
                continue  # Rosetta header.
            if not math.isfinite(value):
                raise RuntimeError('Nonfinite original score')
            values.append(value)
    pdbs = sorted((run / 'Results').rglob('*.pdb'))
    if not values or not pdbs:
        raise RuntimeError('Empty scores/clusters; not a completed successful task')
    return {'Patchdock_Results/score.sc': sha(score),
            **{str(p.relative_to(run)): sha(p) for p in pdbs}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--retry-failed', action='store_true')
    parser.add_argument('--task', action='append', help='Exact task ID; repeat to select several')
    args = parser.parse_args()
    root = args.root.resolve()
    cfg = json.loads((root / 'engine_config.json').read_text())
    if any(any(c.isspace() or c in ';|&$`<>\"\'()' for c in cfg[key])
           for key in ('source','patchdock','rosetta','openbabel','python')):
        raise ValueError('Original shell commands require plain whitespace-free installation paths')
    if cfg['workers'] < 1 or cfg['task_timeout_seconds'] <= 0 or cfg['min_free_bytes'] < 0:
        raise ValueError('Invalid execution limits')
    manifest = json.loads((root / 'TASKS.json').read_text())
    validate_tasks(root, manifest)
    tasks = [t for t in manifest['tasks'] if args.task is None or t['id'] in args.task]
    if not tasks or (args.task and set(args.task)-{t['id'] for t in tasks}):
        raise ValueError('Unknown/empty task selection')
    print(json.dumps({'engine':'Original LondonLab pRosettaC: PatchDock + standalone Rosetta',
        'tasks':len(tasks), 'execute':args.execute, 'sampling':manifest['sampling'],
        'workers':cfg['workers'], 'no_MD':True}, indent=2), flush=True)
    if not args.execute:
        return
    if any(c.isspace() or c in ';|&$`<>\"\'()' for c in str(root)):
        raise ValueError('Original engine requires a plain whitespace-free execution directory')
    for relative, digest in json.loads((root/'BUNDLE_MANIFEST.json').read_text())['files'].items():
        path = (root/relative).resolve()
        path.relative_to(root)
        if sha(path) != digest:
            raise RuntimeError('Bundle changed after preparation: '+relative)
    report = inspect(cfg)
    print(json.dumps(report, indent=2), flush=True)
    if report['status'] != 'PATHS_PASS':
        raise RuntimeError('Original-engine path/dependency check failed; no tasks launched')
    with (root / 'campaign.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for task in tasks:
            if (root / 'STOP').exists():
                raise RuntimeError('STOP sentinel present')
            validate_tasks(root, {'tasks':[task]})
            if shutil.disk_usage(root).free < cfg['min_free_bytes']:
                raise RuntimeError('Insufficient free storage')
            taskroot = root / 'runs' / task['id']
            previous = sorted(taskroot.glob('attempt_*')) if taskroot.exists() else []
            if previous:
                receipt = previous[-1] / 'TASK_STATUS.json'
                status = json.loads(receipt.read_text()) if receipt.exists() else {}
                if status.get('status') == 'COMPLETE':
                    if status.get('input_files') != task['files']:
                        raise RuntimeError('Completed task input identity changed')
                    if status.get('outputs') != outputs(previous[-1]):
                        raise RuntimeError('Completed task outputs changed')
                    print('SKIP verified complete: ' + task['id'], flush=True)
                    continue
                if not args.retry_failed:
                    raise RuntimeError('Incomplete/failed attempt preserved; inspect it, then use --retry-failed')
            run = taskroot / ('attempt_' + str(len(previous)+1).zfill(3))
            shutil.copytree(root / task['input'], run)
            state = dict(task=task['id'], status='RUNNING', input_files=task['files'],
                         engine_checks=report, started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
            write(run / 'TASK_STATUS.json', state)
            active = None
            try:
                env = {**os.environ, 'PATH':str(Path(cfg['python']).parent)+os.pathsep+os.environ['PATH'],
                       'OMP_NUM_THREADS':'1', 'OPENBLAS_NUM_THREADS':'1', 'MKL_NUM_THREADS':'1',
                       'PYTHONDONTWRITEBYTECODE':'1'}
                with (run / 'console.log').open('w') as log:
                    command = [cfg['python'], str(root/'run_original.py'), '--config',
                               str(root/'engine_config.json'), '--workers', str(cfg['workers'])]
                    active = subprocess.Popen(command, cwd=run, env=env, stdout=log,
                                              stderr=subprocess.STDOUT, start_new_session=True)
                    started = time.monotonic()
                    while active.poll() is None:
                        if (root/'STOP').exists() or time.monotonic()-started > cfg['task_timeout_seconds']:
                            raise RuntimeError('STOP requested or task time limit reached')
                        if shutil.disk_usage(root).free < cfg['min_free_bytes']:
                            raise RuntimeError('Free-storage limit reached')
                        time.sleep(2)
                    if active.returncode:
                        raise RuntimeError(f'Original engine exit code {active.returncode}; see console.log')
                state.update(status='COMPLETE', outputs=outputs(run))
            except BaseException as exc:
                state.update(status='FAILED', error=str(exc))
                raise
            finally:
                if active is not None and active.poll() is None:
                    os.killpg(active.pid, signal.SIGTERM)
                    try:
                        active.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(active.pid, signal.SIGKILL)
                        active.wait()
                state['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
                write(run/'TASK_STATUS.json', state)
                print(json.dumps({'task':task['id'], 'status':state['status']}), flush=True)


if __name__ == '__main__':
    main()
