"""Build a held, immutable six-compound Ernie pilot draft; cannot submit/run jobs.

The controller refuses execution until a separately qualified release is built.
This builder deliberately always writes a HOLD release gate.
"""
import argparse,hashlib,json,re,shutil
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser(description=__doc__)
 for k in ('qualification','staged','output'):p.add_argument('--'+k,type=Path,required=True)
 a=p.parse_args()
 for k in ('qualification','staged','output'):setattr(a,k,getattr(a,k).resolve())
 if any(c.isspace() for c in str(a.output)):raise ValueError('Use a whitespace-free bundle directory')
 a.output.mkdir(parents=True,exist_ok=False)
 frozen=json.loads((a.qualification/'BUNDLE_MANIFEST.json').read_text())
 for name,digest in frozen['files'].items():
  if sha(a.qualification/name)!=digest:raise ValueError('Qualification source changed')
 shutil.copytree(a.qualification/'workflow',a.output/'workflow',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
 manifest=json.loads((a.staged/'INPUT_MANIFEST.json').read_text());tasks=[]
 for row in manifest['controls']:
  name=row['compound']+'_'+row['target']
  if not re.fullmatch(r'[A-Za-z0-9_]+',name):raise ValueError('Unsafe task identifier')
  source=a.staged/row['directory'];source.resolve().relative_to(a.staged)
  for f,digest in row['input_file_sha256'].items():
   if sha(source/f)!=digest:raise ValueError('Staged input changed')
  rel=row['compound']+'/'+row['target'];shutil.copytree(source,a.output/'inputs'/rel)
  tasks.append({'name':'global_'+name+'_seed42','reference':rel,'seed':42,'protocol':'global','n_dock':50,'n_local':10,'conformers':20,'max_global_trials':20000,'timeout_seconds':14400,'role':row['role']})
 cfg={'purpose':'Six-compound static technical/resource pilot only, not kinome production or biological validation','cpu_limit':4,'poll_seconds':15,'max_group_RSS_bytes':24*1024**3,'max_output_bytes':24*1024**3,'min_free_bytes':50*1024**3,'tasks':tasks}
 (a.output/'PILOT_CONFIG.json').write_text(json.dumps(cfg,indent=2)+'\n')
 controller=(a.qualification/'run_pilot.py').read_text()
 target="    root=a.root.resolve();lock=(root/'pilot.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)"
 replacement="""    root=a.root.resolve()
    gate=json.loads((root/'RELEASE_GATE.json').read_text())
    if gate.get('status')!='QUALIFIED_FOR_SUBMISSION':
        raise RuntimeError('HOLD: draft inputs/protocol are not qualified; no jobs launched')
    lock=(root/'pilot.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)"""
 if target not in controller:raise ValueError('Unexpected controller version')
 controller=controller.replace(target,replacement).replace("'scope':'static reference/resource pilot only'","'scope':'six-compound technical/resource pilot only'").replace('REFERENCE_PILOT_COMPLETE_NEEDS_RECOVERY_ASSESSMENT','SIX_COMPOUND_TECHNICAL_PILOT_COMPLETE_NEEDS_ASSESSMENT')
 (a.output/'run_pilot.py').write_text(controller)
 gate={'status':'HOLD','production_ready':False,'reason':['Reference recovery qualification pending','Per-input steric/state qualification pending','MS143 PTM-aware post-QC integration pending','Ernie resource/scheduling mode and final submission approval pending'],'note':'No full kinome campaign or MD; builder cannot authorize release'}
 (a.output/'RELEASE_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
 plan={'remote_root':'/home/s1638308/six_protac_ernie_pilot_v1','runtime':'/home/s1638308/miniconda3/envs/protac/bin/python3','controller_argv':['/home/s1638308/miniconda3/envs/protac/bin/python3','/home/s1638308/six_protac_ernie_pilot_v1/run_pilot.py','--root','/home/s1638308/six_protac_ernie_pilot_v1'],'status':'DRAFT_HELD_NOT_SUBMITTED','tasks':len(tasks),'qualification_bundle':str(a.qualification),'qualification_manifest_sha256':sha(a.qualification/'BUNDLE_MANIFEST.json'),'staged_input_manifest_sha256':sha(a.staged/'INPUT_MANIFEST.json'),'max_sequential_task_wall_hours':24,'max_simultaneous_cores':4,'scheduler_policy_note':'Read-only PATH discovery does not establish allocation or scheduling policy','approval':'Present a qualified concrete package to user before remote submission'}
 (a.output/'SUBMISSION_PLAN.json').write_text(json.dumps(plan,indent=2)+'\n')
 m={'scope':'Held six-compound technical/resource draft','files':{str(f.relative_to(a.output)):sha(f) for f in sorted(a.output.rglob('*')) if f.is_file()}}
 (a.output/'BUNDLE_MANIFEST.json').write_text(json.dumps(m,indent=2)+'\n');print(json.dumps(plan,indent=2))
if __name__=='__main__':main()
