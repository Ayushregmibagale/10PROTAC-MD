#!/usr/bin/env python3
"""Run one explicitly requested local benchmark and record terminal validation.

This is a benchmark harness, not an Ernie submission mechanism. Existing attempts
cannot be reused. It never upgrades scientific-release status or forgives missing
outputs. Use an engine copy that reports failures rather than placeholder scores.
"""
import argparse
import fcntl
import json
import os
import signal
import subprocess
import sys
import time
import uuid
import threading
from contextlib import contextmanager
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import atomic_json, sha256, validate_outputs, ValidationError, verify_execution_snapshot

@contextmanager
def interruption_signals():
    previous={}
    if threading.current_thread() is threading.main_thread():
        def interrupted(signum,frame):raise KeyboardInterrupt('Received signal '+str(signum))
        for signum in [signal.SIGTERM,signal.SIGHUP]:previous[signum]=signal.signal(signum,interrupted)
    try:yield
    finally:
        for signum,handler in previous.items():signal.signal(signum,handler)

def record(run,command,timeout):
    run=Path(run).resolve()
    if not command or timeout<=0: raise ValidationError('Command and positive timeout required')
    with (run/'.attempt.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if (run/'attempt.json').exists() or (run/'completion.json').exists() or (run/'Patchdock_Results').exists():
            raise ValidationError('Attempt directory already used; create a fresh isolated copy')
        snapshot=run/'source_snapshot'
        if snapshot.exists():
            try:
                manifest=json.loads((snapshot/'source_manifest.json').read_text())
            except (OSError,ValueError) as exc:
                raise ValidationError('Execution source manifest required before starting a benchmark') from exc
            verify_execution_snapshot(run,{'engine_sha256':manifest.get('source_sha256')})
        selection=json.loads((run/'selection.json').read_text())
        inputs={p.name:sha256(p) for p in run.iterdir() if p.is_file() and p.name!='.attempt.lock'}
        receipt={'schema_version':1,'attempt_id':str(uuid.uuid4()),'command':command,
                 'recorder_sha256':sha256(Path(__file__)),
                 'input_sha256':inputs,'selection_sha256':sha256(run/'selection.json'),
                 'started_unix':time.time(),'status':'RUNNING','exit_code':None}
        atomic_json(run/'attempt.json',receipt)
        status='FAILED'; reason=''; code=None
        try:
            with (run/'engine.log').open('w') as log, interruption_signals():
                process=subprocess.Popen(command,cwd=run,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
                receipt.update(process_id=process.pid,process_group_id=process.pid)
                atomic_json(run/'attempt.json',receipt)
                try: code=process.wait(timeout=timeout)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid,signal.SIGTERM)
                    try: process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid,signal.SIGKILL); process.wait()
                    raise ValidationError('Benchmark timeout; process group terminated')
                except BaseException:
                    # Ctrl-C must not leave expensive workers orphaned.
                    os.killpg(process.pid,signal.SIGTERM)
                    try:process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid,signal.SIGKILL);process.wait()
                    raise
            if code!=0: raise ValidationError(f'Engine exit code {code}')
            for name,digest in inputs.items():
                if sha256(run/name)!=digest: raise ValidationError('Input mutated during run: '+name)
            rows=validate_outputs(run,selection)
            receipt['engine_status_sha256']=sha256(run/'engine_status.json')
            receipt['score_sha256']=sha256(run/'Patchdock_Results/score.sc')
            receipt['pose_sha256']={r['description']:r['pdb_sha256'] for r in rows}
            receipt['n_scored_valid']=len(rows); status='COMPLETE' if rows else 'COMPLETE_NO_VALID_POSES'
        except KeyboardInterrupt:
            status='INTERRUPTED';reason='Interrupted; child process group terminated'
            code=process.returncode if 'process' in locals() else None
        except (ValidationError,OSError,ValueError) as exc: reason=str(exc)
        receipt.update(status=status,reason=reason,exit_code=code,finished_unix=time.time())
        atomic_json(run/'completion.json',receipt)
        return receipt

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-dir',type=Path,required=True)
    p.add_argument('--timeout-seconds',type=float,required=True)
    p.add_argument('command',nargs=argparse.REMAINDER)
    a=p.parse_args(); command=a.command[1:] if a.command[:1]==['--'] else a.command
    result=record(a.run_dir,command,a.timeout_seconds)
    print(json.dumps(result,indent=2)); raise SystemExit(0 if result['status'] in {'COMPLETE','COMPLETE_NO_VALID_POSES'} else 1)
