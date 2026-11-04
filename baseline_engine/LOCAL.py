"""
Local multiprocessing cluster for pRosettaC.
Replaces PBS/SGE/SLURM when running on a single workstation.
"""
import os
import json,hashlib,shlex
from pathlib import Path
import subprocess
import concurrent.futures
from ..Cluster import Cluster


class LOCAL(Cluster):
    def __init__(self, n_workers=None):
        self.n_workers = n_workers or max(1, int(os.environ.get('PROTAC_WORKERS', '2')))

    def getRunningJobIDs(self):
        return []

    def submit(self, job_file):
        subprocess.run(['bash', job_file], check=True)
        return job_file

    def writeJobFile(self, job_file, commands, mem):
        with open(job_file, 'w') as f:
            f.write('#!/bin/bash\nset -e\n')
            for cmd in commands:
                f.write(cmd + '\n')
        os.chmod(job_file, 0o755)

    def runBatchCommands(self, commands, mem, batch_size=12):
        """Run all commands in parallel using a thread pool (threads share CWD and avoid
        pickle issues that arise with ProcessPoolExecutor for local functions)."""
        def _run(cmd):
            root=Path(os.environ.get('PROTAC_ATTEMPT_ROOT',os.getcwd())).resolve()
            words=shlex.split(cmd)
            normalized=[]
            for index,word in enumerate(words[1:]):
                path=Path(word)
                if index==0: normalized.append(path.name)
                elif path.is_absolute():
                    try: normalized.append(str(path.resolve().relative_to(root)))
                    except ValueError: normalized.append(path.name)
                else: normalized.append(word)
            task_key=hashlib.sha256(json.dumps(normalized).encode()).hexdigest()
            env=dict(os.environ,PROTAC_TASK_KEY=task_key)
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True,env=env)
            logs=Path('LOCAL_logs');logs.mkdir(exist_ok=True)
            (logs/(task_key+'.log')).write_text(result.stdout+result.stderr)
            (logs/(task_key+'.json')).write_text(json.dumps({'command':cmd,'normalized_task':normalized,'task_key':task_key,'exit_code':result.returncode},indent=2)+'\n')
            if result.returncode != 0:
                print(f'[LOCAL] Command failed (rc={result.returncode}): {cmd}\n{result.stderr[-500:]}')
            result.check_returncode()
            return result.returncode

        with concurrent.futures.ThreadPoolExecutor(max_workers=self.n_workers) as ex:
            futures = {ex.submit(_run, cmd): cmd for cmd in commands}
            for fut in concurrent.futures.as_completed(futures):
                try:
                    fut.result()
                except Exception:
                    for pending in futures: pending.cancel()
                    raise
        return ['local_done']

    def wait(self, job_ids, timeout=-1):
        pass  # Already waited in runBatchCommands

    def jobsRunning(self, job_ids):
        return False
