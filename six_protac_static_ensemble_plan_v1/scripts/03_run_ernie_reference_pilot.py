"""Bounded four-core static reference pilot; no kinome production or MD.

Runs only frozen inputs and source hashes, with a lock, CPU affinity, group
termination, timeout, memory/storage limits and durable terminal receipts.
"""
import argparse,fcntl,hashlib,json,os,signal,subprocess,sys,time,shutil
from pathlib import Path

def atomic(path,data):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');os.replace(tmp,path)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(root,manifest):
    for rel,digest in manifest['files'].items():
        if sha(root/rel)!=digest:raise RuntimeError('Frozen file changed: '+rel)
def memory_bytes(pgid):
    total=0
    for proc in Path('/proc').iterdir():
        if not proc.name.isdigit():continue
        try:
            fields=(proc/'stat').read_text().rsplit(')',1)[1].split()
            if int(fields[2])==pgid:total+=int((proc/'statm').read_text().split()[1])*os.sysconf('SC_PAGE_SIZE')
        except (OSError,ValueError,IndexError):pass
    return total
def disk_bytes(root):
    total=0
    for p in root.rglob('*'):
        try:
            if p.is_file():total+=p.stat().st_size
        except OSError:pass
    return total
def stop(proc):
    try:os.killpg(proc.pid,signal.SIGTERM)
    except ProcessLookupError:return
    try:proc.wait(timeout=10)
    except subprocess.TimeoutExpired:
        try:os.killpg(proc.pid,signal.SIGKILL)
        except ProcessLookupError:pass
        proc.wait(timeout=10)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args()
    root=a.root.resolve();lock=(root/'pilot.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if (root/'STATUS.json').exists():raise RuntimeError('Existing attempt: inspect receipts; no blind relaunch')
    cfg=json.loads((root/'PILOT_CONFIG.json').read_text());manifest=json.loads((root/'BUNDLE_MANIFEST.json').read_text());verify(root,manifest)
    permitted=sorted(os.sched_getaffinity(0));cpus=permitted[:cfg['cpu_limit']]
    if len(cpus)!=cfg['cpu_limit']:raise RuntimeError('Insufficient available CPUs')
    os.sched_setaffinity(0,cpus);os.nice(10)
    env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','PYTHONDONTWRITEBYTECODE':'1'}
    state={'status':'RUNNING','scope':'static reference/resource pilot only','MD':False,'pid':os.getpid(),
           'started_unix':time.time(),'cpus':cpus,'configuration':cfg,'tasks':[]}
    atomic(root/'STATUS.json',state);active=None
    try:
        (root/'runs').mkdir();(root/'logs').mkdir();(root/'receipts').mkdir()
        for task in cfg['tasks']:
            verify(root,manifest)
            if (root/'STOP').exists():raise RuntimeError('STOP sentinel present')
            if shutil.disk_usage(root).free<cfg['min_free_bytes']:raise RuntimeError('Insufficient free storage')
            run=root/'runs'/task['name'];shutil.copytree(root/'inputs'/task['reference'],run)
            command=[sys.executable,str(root/'workflow/baseline_engine/run.py'),'--run-dir',str(run),
                '--seed',str(task['seed']),'--protocol',task['protocol'],'--n-dock',str(task['n_dock']),
                '--n-local',str(task['n_local']),'--conformers',str(task['conformers']),
                '--workers',str(cfg['cpu_limit']),'--max-global-trials',str(task['max_global_trials'])]
            entry={'name':task['name'],'command':command,'started_unix':time.time(),'status':'RUNNING'}
            state['tasks'].append(entry);atomic(root/'STATUS.json',state)
            with (root/'logs'/f'{task["name"]}.log').open('w') as log:
                active=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
                entry['pid']=active.pid;entry['process_group']=active.pid;atomic(root/'STATUS.json',state)
                peakrss=0;peakdisk=0;reason=None
                while active.poll() is None:
                    rss=memory_bytes(active.pid);used=disk_bytes(root/'runs');peakrss=max(peakrss,rss);peakdisk=max(peakdisk,used)
                    if (root/'STOP').exists():reason='STOP_REQUESTED'
                    elif time.time()-entry['started_unix']>task['timeout_seconds']:reason='TIME_LIMIT'
                    elif rss>cfg['max_group_RSS_bytes']:reason='MEMORY_LIMIT'
                    elif used>cfg['max_output_bytes'] or shutil.disk_usage(root).free<cfg['min_free_bytes']:reason='STORAGE_LIMIT'
                    if reason:stop(active);break
                    entry.update(elapsed_seconds=time.time()-entry['started_unix'],peak_group_RSS_bytes=peakrss,peak_output_bytes=peakdisk)
                    state['heartbeat_unix']=time.time();atomic(root/'STATUS.json',state);time.sleep(cfg.get('poll_seconds',15))
                entry.update(exit_code=active.wait(),finished_unix=time.time(),stop_reason=reason)
                engine=run/'engine_status.json'
                if engine.exists():
                    report=json.loads(engine.read_text());entry['engine_status']=report['status'];entry['stages']=report.get('stages',{});entry['engine_status_sha256']=sha(engine)
                entry['status']='COMPLETE' if entry['exit_code']==0 and reason is None and entry.get('engine_status') in ('COMPLETE','COMPLETE_NO_VALID_POSES') else 'FAILED_OR_STOPPED'
                active=None;atomic(root/'receipts'/f'{task["name"]}.json',entry);atomic(root/'STATUS.json',state)
                if entry['status']!='COMPLETE':raise RuntimeError('Task failed/stopped; retained all artifacts: '+task['name'])
            verify(root,manifest)
        state['status']='REFERENCE_PILOT_COMPLETE_NEEDS_RECOVERY_ASSESSMENT'
    except BaseException as exc:
        if active is not None:
            stop(active)
            entry.update(status='FAILED_OR_STOPPED',error=repr(exc),exit_code=active.returncode,finished_unix=time.time())
            atomic(root/'receipts'/f'{entry["name"]}.json',entry)
        state.update(status='FAILED_OR_STOPPED',error=repr(exc));raise
    finally:
        state['finished_unix']=time.time();atomic(root/'STATUS.json',state)

if __name__=='__main__':main()
