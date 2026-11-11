"""One-time pretraining integrity amendment after the current extraction exits.
Preserve successful graphs, every original input, stopped source and frozen protocol.
"""
from pathlib import Path
import json,sys,time,hashlib,shutil,os,subprocess,datetime
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts/lib'))
from protac_qc import atomic_json
PLAN=ROOT/'recovery/input_validity_plan_v4.json'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 import fcntl
 lock=(ROOT/'recovery/activation.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 plan=json.loads(PLAN.read_text());start=time.monotonic()
 while True:
  state=json.loads((ROOT/'JOB_STATUS.json').read_text())
  if state['phase']=='STOPPED_ON_ERROR':break
  if state['phase']!='01_build_static_graphs.py':raise RuntimeError('Unexpected old phase; do not alter an active training run')
  if time.monotonic()-start>24*3600:raise RuntimeError('Old extraction did not exit within 24 hours')
  time.sleep(30)
 if '01_build_static_graphs.py failed with exit' not in state.get('error',''):raise RuntimeError('Old failure does not match the inspected extraction failure')
 if any((ROOT/'analysis/gnn').rglob('fold_*.json')):raise RuntimeError('Never amend the protocol after model results')
 if sha(ROOT/'config/protocol.json')!=plan['original_protocol_sha256']:raise RuntimeError('Original frozen protocol drift')
 if any(sha(ROOT/p)!=digest for p,digest in plan['candidate_sha256'].items()):raise RuntimeError('Reviewed candidate source drift')
 tests=(ROOT/'logs/input_validity_v4_tests_20261111.txt').read_text()
 if '14 passed' not in tests:raise RuntimeError('Fourteen candidate tests must pass')
 archive=ROOT/'logs/pre_input_validity_amendment_20261111';archive.mkdir(exist_ok=False)
 for relative in ['scripts/01_build_static_graphs.py','scripts/04_train_dual_lopo.py','scripts/07_continue_reanalysis.py','scripts/06_make_comparison.py','config/protocol.json','JOB_STATUS.json']:
  destination=archive/relative;destination.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/relative,destination)
 protocol=json.loads((ROOT/'config/protocol.json').read_text());protocol['static_graph']['input_validity_policy']=plan['input_validity_policy']
 for name in ['01_build_static_graphs.py','04_train_dual_lopo.py','06_make_comparison.py','07_continue_reanalysis.py']:shutil.copy2(ROOT/'candidate_scripts_v4'/name,ROOT/'scripts'/name)
 test=(ROOT/'candidate_tests_v4/test_static_integrity.py').read_text().replace("ROOT/'candidate_scripts_v4/01_build_static_graphs.py'","ROOT/'scripts/01_build_static_graphs.py'")
 (ROOT/'tests/test_static_integrity.py').write_text(test)
 atomic_json(ROOT/'config/protocol.json',protocol)
 receipt={'status':'PRETRAINING_INTEGRITY_AMENDMENT_APPLIED','UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'old_status':state,'original_protocol_sha256':plan['original_protocol_sha256'],'amended_protocol_sha256':sha(ROOT/'config/protocol.json'),'original_files_modified':False,'new_docking_or_MD':False,'training_results_seen':False,'plan_sha256':sha(PLAN)}
 atomic_json(ROOT/'logs/input_validity_activation_v4.json',receipt)
 log=ROOT/'logs/continuation_v4_20261111.log'
 with log.open('x') as output:
  result=subprocess.run(['/home/ayushregmibagale/miniconda3/bin/python3',str(ROOT/'scripts/07_continue_reanalysis.py')],stdout=output,stderr=subprocess.STDOUT)
 if result.returncode:raise RuntimeError('Amended chain stopped; see '+str(log))
if __name__=='__main__':
 try:main()
 except Exception as exc:
  atomic_json(ROOT/'recovery/activation_error_v4.json',{'error':repr(exc),'UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()});raise
