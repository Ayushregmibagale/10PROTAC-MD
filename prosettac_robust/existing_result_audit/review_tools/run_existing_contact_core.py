from pathlib import Path
import concurrent.futures,subprocess,sys,json,time,hashlib
run=Path('/tmp/10protac_finish/prosettac_robust/existing_result_audit/contact_core_v1');plan=json.loads((run/'plan.json').read_text());snapshot=run/'source_snapshot';script=snapshot/'existing_result_audit/review_tools/recheck_existing_trajectory_contacts.py'
start=time.time();outcomes=[]
def work(case):
 output=run/'poses'/case['poi']/case['pose_id'];output.parent.mkdir(parents=True,exist_ok=True);log=output.with_suffix('.log');command=[sys.executable,str(script),'--project',plan['project'],'--poi',case['poi'],'--pose-id',case['pose_id'],'--time-end-ns',str(plan['time_end_ns']),'--output',str(output)]
 with log.open('w') as f:result=subprocess.run(command,text=True,stdout=f,stderr=subprocess.STDOUT,timeout=600)
 return dict(case,exit_code=result.returncode,output=str(output),log=str(log))
with concurrent.futures.ThreadPoolExecutor(max_workers=plan['workers']) as executor:
 for future in concurrent.futures.as_completed([executor.submit(work,case) for case in plan['cases']]):
  try:outcomes.append(future.result())
  except Exception as exc:outcomes.append({'exit_code':-1,'error':str(exc)})
  if len(outcomes)%10==0:print(json.dumps({'completed':len(outcomes),'failed':sum(x['exit_code']!=0 for x in outcomes),'elapsed_s':round(time.time()-start,1)}),flush=True)
manifest=json.loads((snapshot/'source_manifest.json').read_text());assert all(hashlib.sha256((snapshot/name).read_bytes()).hexdigest()==digest for name,digest in manifest['files_sha256'].items())
summary={'status':'CONTACT_REANALYSIS_COMPLETE' if all(x['exit_code']==0 for x in outcomes) else 'CONTACT_REANALYSIS_REVIEW_REQUIRED','cases':len(outcomes),'failed':sum(x['exit_code']!=0 for x in outcomes),'elapsed_s':time.time()-start,'outcomes':outcomes,'source_manifest_sha256':hashlib.sha256((snapshot/'source_manifest.json').read_bytes()).hexdigest(),'plan_sha256':hashlib.sha256((run/'plan.json').read_bytes()).hexdigest(),'original_files_modified':False,'new_MD_or_docking_executed':False};(run/'completion.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='outcomes'},indent=2),flush=True)
