"""Monitored corrected-label MD-only execution with fail-closed stage gates."""
from pathlib import Path
import sys,os,json,time,subprocess,datetime,fcntl
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'))
from protac_qc import atomic_json
CANONICAL=Path('/home/ayushregmibagale/10PROTAC-MD/SJF8240_reanalysis/md_weighted_v1')
PYTHON='/home/ayushregmibagale/miniconda3/envs/boltz2/bin/python'
CPU_PYTHON='/home/ayushregmibagale/miniconda3/bin/python3'
def status(phase,**details):
 record={'updated_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'phase':phase,'pid':os.getpid(),'work_root':str(ROOT),'cohort_size':8,'positive_count':6,'negative_count':2,'EPHA2_degraded':True,'folds_complete':{role:len(list((ROOT/'analysis'/role).glob('fold_*.json'))) for role in ['elob','vhl']},'full_reanalysis_complete':phase=='COMPLETE','active_job':phase not in ['COMPLETE','STOPPED_ON_ERROR'],**details};atomic_json(ROOT/'STATUS.json',record)
 if CANONICAL.exists():
  atomic_json(CANONICAL/'STATUS.json',record)
  parent=CANONICAL.parent/'STATUS.json';old=json.loads(parent.read_text())
  old.update({'updated_UTC':record['updated_UTC'],'phase':'CORRECTED_MD_ONLY_'+phase,'current_label_reanalysis_complete':record['full_reanalysis_complete'],'full_reanalysis_complete':record['full_reanalysis_complete'],'active_job':record['active_job'],'current_run':'md_weighted_v1','current_run_work_root':str(ROOT),'folds_complete':record['folds_complete'],'completed_model_arms':sum((ROOT/'analysis'/r/'completion.json').exists() for r in ['elob','vhl']),'current_findings':'md_weighted_v1/RESULTS.md' if phase=='COMPLETE' else 'CURRENT_LABELS.md'})
  atomic_json(parent,old)
 print(json.dumps(record),flush=True)
def execute(script,args,log,cpu=False):
 env=os.environ.copy();
 if cpu:env.pop('PYTHONPATH',None)
 else:env['PYTHONPATH']='/mnt/c/Users/Ayush Regmi Bagale/SJF8240_reanalysis_work/runtime'
 env['PYTHONUNBUFFERED']='1'
 with (ROOT/'logs'/log).open('x') as stream:
  process=subprocess.Popen([CPU_PYTHON if cpu else PYTHON,str(ROOT/'scripts'/script),*args],env=env,stdout=stream,stderr=subprocess.STDOUT)
  while process.poll() is None:status(script,child_pid=process.pid,arguments=args,log='logs/'+log);time.sleep(10)
  if process.returncode:raise RuntimeError(script+' failed with exit '+str(process.returncode)+'; see logs/'+log)
def main():
 lock=(ROOT/'logs/pipeline.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 execute('preflight.py',[],'preflight.log')
 for role in ['elob','vhl']:execute('train.py',['--role',role,'--device','cuda'],'train_'+role+'.log')
 execute('verify_and_report.py',[],'independent_review.log',cpu=True);status('COMPLETE',independent_artifact_review_complete=True,independent_biological_validation=False,results='RESULTS.md')
if __name__=='__main__':
 try:main()
 except Exception as exc:status('STOPPED_ON_ERROR',error=repr(exc));raise
