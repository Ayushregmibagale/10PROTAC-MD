"""Inspect existing compiled production settings; never run grompp or mdrun.

GROMACS dump -om output is an inspection artifact, not a reusable input MDP.
"""
from pathlib import Path
import concurrent.futures,collections,json,subprocess,sys,re
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import sha256,atomic_json,ValidationError
project=Path('/home/ayushregmibagale/protac-md');out=ROOT/'existing_result_audit/compiled_md_protocol_v1';gmx=Path('/home/ayushregmibagale/.conda/envs/gromacs/bin/gmx');roster=ROOT/'sjf8240_review/md_graph_pose_lineage.json'
if out.exists():raise ValidationError('Fresh output required')
source_hash=sha256(Path(__file__));gmx_hash=sha256(gmx);roster_hash=sha256(roster);cases=json.loads(roster.read_text())
if len(cases)!=174 or len({(r['poi'],r['pose_id']) for r in cases})!=174:raise ValidationError('Expected 174 unique saved MD graph cases')
out.mkdir(parents=True)
keys=['integrator','dt','nsteps','nstxout-compressed','rcoulomb','rvdw','rlist','coulombtype','coulomb-modifier','vdw-type','vdw-modifier','DispCorr','cutoff-scheme','tcoupl','ref-t','tau-t','pcoupl','pcoupltype','ref-p','tau-p','pbc','continuation','constraint-algorithm','lincs-order','lincs-iter']
def worker(case):
 poi=case['poi'];pose=case['pose_id'];tpr=project/'random_screen'/poi/pose/'prod/prod.tpr';before=sha256(tpr);folder=out/poi/pose;folder.mkdir(parents=True);mdp=folder/'compiled_input_inspection.mdp';log=folder/'dump.log';command=[str(gmx),'dump','-s',str(tpr),'-om',str(mdp)]
 with log.open('w') as f:run=subprocess.run(command,stdout=subprocess.DEVNULL,stderr=f,text=True,timeout=60)
 if run.returncode:raise ValidationError('GROMACS compiled-input inspection failed: '+str(log))
 settings={}
 for line in mdp.read_text().splitlines():
  if '=' not in line:continue
  name,value=line.split('=',1);name=name.strip().replace('_','-');value=' '.join(value.split())
  if name in settings:raise ValidationError('Duplicate compiled setting '+name)
  settings[name]=value
 selected={k:settings.get(k) for k in keys}
 if any(selected[k] is None for k in ['dt','nsteps','nstxout-compressed','rcoulomb','rvdw']):raise ValidationError('Missing compiled protocol field')
 row={'poi':poi,'pose_id':pose,'compiled_settings':selected,'nominal_duration_ns':float(settings['dt'])*int(settings['nsteps'])/1000,'stored_frame_interval_ps':float(settings['dt'])*int(settings['nstxout-compressed']),'tpr':str(tpr),'tpr_sha256':before,'inspection_sha256':sha256(mdp),'log_sha256':sha256(log),'command':command}
 if sha256(tpr)!=before:raise ValidationError('TPR changed during read-only inspection')
 return row
rows=[];failures=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 futures={pool.submit(worker,r):r for r in cases}
 for f in concurrent.futures.as_completed(futures):
  try:rows.append(f.result())
  except Exception as exc:failures.append(dict(poi=futures[f]['poi'],pose_id=futures[f]['pose_id'],error=str(exc)))
  if (len(rows)+len(failures))%30==0:print(json.dumps({'completed':len(rows),'failed':len(failures)}),flush=True)
if sha256(gmx)!=gmx_hash or sha256(roster)!=roster_hash or sha256(Path(__file__))!=source_hash:raise ValidationError('Audit source changed')
rows.sort(key=lambda r:(r['poi'],r['pose_id']));atomic_json(out/'per_pose_protocol.json',rows)
classes=collections.Counter(json.dumps(r['compiled_settings'],sort_keys=True) for r in rows)
summary={'status':'COMPILED_PROTOCOL_INSPECTION_COMPLETE' if not failures else 'INCOMPLETE_COMPILED_PROTOCOL_INSPECTION','requested':len(cases),'inspected':len(rows),'failures':failures,'protocol_classes':[{'cases':n,'settings':json.loads(k)} for k,n in classes.items()],'nominal_duration_ns_counts':dict(collections.Counter(str(r['nominal_duration_ns']) for r in rows)),'stored_frame_interval_ps_counts':dict(collections.Counter(str(r['stored_frame_interval_ps']) for r in rows)),'source_sha256':source_hash,'gmx_sha256':gmx_hash,'roster_sha256':roster_hash,'records_sha256':sha256(out/'per_pose_protocol.json'),'new_MD_or_model_training_executed':False,'original_files_modified':False,'limitations':['Compiled input settings do not establish completed steps, convergence or independent replicas.','Generated inspection MDP files cannot be reused as grompp inputs.','Absent template-level constraint selection cannot be reconstructed merely from the inputrec; bonded interactions require separate inspection.']}
atomic_json(out/'summary.json',summary);print(json.dumps({k:v for k,v in summary.items() if k not in ['protocol_classes','limitations']},indent=2))
