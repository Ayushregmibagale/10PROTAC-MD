"""Monitor retained static extraction, then execute the frozen paired GPU analysis.
No docking/MD; explicit failures stop the chain. Persistent JSON reports actual state.
"""
from pathlib import Path
import os,sys,json,time,subprocess,hashlib,datetime,shutil,traceback
ROOT=Path(__file__).resolve().parents[1]
CANONICAL=Path('/home/ayushregmibagale/10PROTAC-MD/SJF8240_reanalysis')
sys.path.insert(0,str(ROOT/'scripts/lib'))
from protac_qc import atomic_json
PYTHON='/home/ayushregmibagale/miniconda3/envs/boltz2/bin/python'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def status(phase,**extra):
 static=list((ROOT/'analysis/static').glob('*/receipt.json'))
 folds={r:len(list((ROOT/'analysis/gnn'/r/'seed_42').glob('fold_*.json'))) for r in ['elob','vhl']}
 arms=sum((ROOT/'analysis/gnn'/r/'seed_42/completion.json').exists() for r in folds)
 record={'updated_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'phase':phase,'pid':os.getpid(),'active_work_root':str(ROOT),'static_targets_complete':len(static),'static_targets_expected':19,'dynamic_cases_complete':174,'dynamic_cases_expected':174,'folds_complete':folds,'completed_model_arms':arms,'full_reanalysis_complete':phase=='PAIRED_REANALYSIS_COMPLETE',**extra}
 atomic_json(ROOT/'JOB_STATUS.json',record);atomic_json(ROOT/'STATUS.json',record)
 atomic_json(CANONICAL/'STATUS.json',record)
 print(json.dumps(record),flush=True)
def sync():
 """Save only this derivative, preserving replaced artifacts in dated backups."""
 backup=CANONICAL/'recovery_backups_20261111'
 for source in ROOT.rglob('*'):
  relative=source.relative_to(ROOT)
  if not source.is_file() or any(x in {'runtime','__pycache__','.pytest_cache'} for x in relative.parts):continue
  target=CANONICAL/relative
  if target.exists() and sha(source)==sha(target):continue
  if target.exists():
   archived=backup/relative
   if not archived.exists():archived.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(target,archived)
  target.parent.mkdir(parents=True,exist_ok=True)
  temporary=target.with_name('.'+target.name+'.recovery_pending');shutil.copy2(source,temporary);os.replace(temporary,target)
def execute(script,args,log,manifest):
 if any(sha(ROOT/p)!=v for p,v in manifest.items()):raise RuntimeError('Frozen execution source changed')
 env=os.environ.copy();env['PYTHONPATH']=str(ROOT/'runtime');env['PYTHONUNBUFFERED']='1'
 with (ROOT/'logs'/log).open('x') as output:
  executable='/home/ayushregmibagale/miniconda3/bin/python3' if script=='01_build_static_graphs.py' else PYTHON
  process=subprocess.Popen([executable,str(ROOT/'scripts'/script),*args],env=env,stdout=output,stderr=subprocess.STDOUT)
  while process.poll() is None:
   status(script,child_pid=process.pid,log='logs/'+log);time.sleep(30)
  if process.returncode:raise RuntimeError(script+' failed with exit '+str(process.returncode)+'; see logs/'+log)
 if any(sha(ROOT/p)!=v for p,v in manifest.items()):raise RuntimeError('Source changed during execution')
def main():
 # Exclusive lock prevents accidental second model chains.
 import fcntl
 lock=(ROOT/'logs/continuation.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 manifest={str(p.relative_to(ROOT)):sha(p) for folder in ['scripts','config'] for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts}
 atomic_json(ROOT/'logs/continuation_source_manifest_20261111.json',manifest)
 sync();execute('01_build_static_graphs.py',['--workers','4'],'static_build_resume_v3_20261111.log',manifest)
 sync();execute('05_check_model_runtime.py',[],'gpu_runtime_resume_20261111.log',manifest)
 for role in ['elob','vhl']:
  execute('04_train_dual_lopo.py',['--role',role,'--device','cuda'],'train_'+role+'_20261111.log',manifest);sync()
 execute('06_make_comparison.py',[],'paired_comparison_20261111.log',manifest)
 status('PAIRED_REANALYSIS_COMPLETE');sync()
if __name__=='__main__':
 try:main()
 except Exception as error:
  status('STOPPED_ON_ERROR',error=repr(error),traceback=traceback.format_exc());raise
