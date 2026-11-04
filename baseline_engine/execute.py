#!/usr/bin/env python3
"""Execute the existing baseline stages with explicit sampling/resource budgets.

Native-frame mode is a diagnostic control. Global mode retains the existing
adapter search. Invoke through run.py for durable output QC and stage receipts.
"""
import argparse,os,sys,json,hashlib,shutil,time,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--run',type=Path,required=True);p.add_argument('--seed',type=int,required=True)
p.add_argument('--protocol',choices=['global','native_frame'],default='global');p.add_argument('--n-dock',type=int,default=200);p.add_argument('--workers',type=int,default=2);p.add_argument('--max-global-trials',type=int,default=10000000);p.add_argument('--local-poses',type=int,default=50);p.add_argument('--conformers',type=int,default=2);a=p.parse_args()
if not 1<=a.local_poses<=50 or not 1<=a.conformers<=100 or a.n_dock<1 or a.workers<1 or a.max_global_trials<1:raise SystemExit('Invalid baseline budgets: local poses 1–50; conformers 1–100; positive global/workers/trial counts')
base=Path(__file__).resolve().parent;run=a.run.resolve();run.mkdir(exist_ok=False,parents=True)
for f in a.input.iterdir():
 if f.is_file() and f.name not in {'attempt.json','completion.json','engine_status.json','engine.log','.attempt.lock'}:shutil.copy2(f,run/f.name)
os.chdir(run)
os.environ.pop('PROTAC_DIAGNOSTIC',None)
os.environ.update(PROTAC_PROSETTAC=str(base/'runtime/PRosettaC'),PROTAC_ATTEMPT_ROOT=str(run),PROTAC_SEED=str(a.seed),PROTAC_WORKERS=str(a.workers),PROTAC_CONFORMERS=str(a.conformers),PROTAC_MAX_GLOBAL_TRIALS=str(a.max_global_trials))
ad=base/'runtime/adapter';pr=base/'runtime/PRosettaC';sys.path[:0]=[str(ad),str(pr)]
import pyrosetta_adapter as adapter
adapter.apply_patches()
import utils,main
sys.path.insert(0,str(base.parent/'scripts/lib'))
from protac_qc import parse_scores
# Original global stage is replaced ONLY for this explicitly native diagnostic.
def native_frame(structs,anchors,min_dist,max_dist,num_results=1000,threshold=2):
 out=Path('Patchdock_Results');out.mkdir(exist_ok=False)
 with (out/'pd.1.pdb').open('w') as f:
  for i,path in enumerate(structs):
   for line in Path(path).read_text().splitlines():
    if line.startswith(('END','MODEL')):continue
    if i==1 and line.startswith('HETATM') and line[17:20].strip()=='PT1':line=line[:21]+'Y'+line[22:]
    f.write(line+'\n')
 return 1
if a.protocol=='native_frame':
    utils.patchdock=native_frame
else:
    original_global=utils.patchdock
    def global_sample(structs,anchors,min_dist,max_dist,num_results=1000,threshold=2):
        return original_global(structs,anchors,min_dist,max_dist,min(num_results,a.n_dock),threshold)
    utils.patchdock=global_sample
original_local=adapter.local_docking
import rosetta
rosetta.local_docking=lambda struct,ca,cb,pa,pb,nstruct=50:original_local(struct,ca,cb,pa,pb,a.local_poses)
intent={'status':'RUNNING','scope':a.protocol+' baseline execution; scientific validation remains separate',
 'requested_local_poses':a.local_poses,'conformers_per_local_pose':a.conformers,
 'seed':a.seed,'started_unix':time.time(),'baseline_runtime_manifest_sha256':hashlib.sha256((base/'runtime_manifest.json').read_bytes()).hexdigest()}
(run/'diagnostic.json').write_text(json.dumps(intent,indent=2)+'\n')
shutil.copy2(base/'runtime_manifest.json',run/'runtime_manifest.json')
try:
    try:main.main('baseline_execution',['Protac_params.txt'])
    except SystemExit as early:
        ledger=run/'global_sampling.json'
        if early.code not in (None,0) or not ledger.exists() or json.loads(ledger.read_text())['succeeded']!=0:raise RuntimeError('Baseline terminated before verified stage completion')
    sampling=sorted((run/'Patchdock_Results').glob('confs_*.sdf.sampling.json'))
    local=sorted((run/'Patchdock_Results').glob('pd.*_docking_????.pdb'))
    if len(sampling)!=len(local) or (not local and not (run/'global_sampling.json').exists()):raise RuntimeError('Missing per-local-pose linker accounting')
    reports=[json.loads(p.read_text()) for p in sampling]
    count=sum(d['generated_conformers']>0 for d in reports)
    intent['linker_sampling_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sampling}
    intent.update(local_pose_count=len(local),local_poses_with_feasible_linkers=count,
                  generated_full_conformers=sum(d['generated_conformers'] for d in reports),biological_label=None)
    score=run/'Patchdock_Results/score.sc'
    if count==0:
        if score.exists() and parse_scores(score,allow_empty=True):raise RuntimeError('Scored output exists despite no feasible linker')
        if not score.exists():score.write_text('SEQUENCE: \nSCORE: total_score description\n')
        intent['status']='NO_FEASIBLE_LINKER_POSES_IN_DIAGNOSTIC'
    else:
        rows=parse_scores(score)
        if len(rows)!=count or not (run/'result_summary.txt').exists():raise RuntimeError('Incomplete scoring/clustering stage')
        intent['status']='MECHANICAL_PIPELINE_COMPLETE_REQUIRES_OUTPUT_QC'
except BaseException as e:
    intent.update(status='FAILED',error=repr(e));raise
finally:
    intent['finished_unix']=time.time();(run/'diagnostic.json').write_text(json.dumps(intent,indent=2)+'\n')
