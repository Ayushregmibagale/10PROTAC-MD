#!/usr/bin/env python3
"""Bounded native-frame integration diagnostic of the existing baseline pipeline.

This deliberately supplies the native global frame and caps local/conformer work.
It is not an independent docking recovery benchmark or a production launcher.
"""
import argparse,os,sys,json,hashlib,shutil,time,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--run',type=Path,required=True);p.add_argument('--seed',type=int,required=True)
p.add_argument('--local-poses',type=int,default=1);p.add_argument('--conformers',type=int,default=2);a=p.parse_args()
if not 1<=a.local_poses<=10 or not 1<=a.conformers<=100:raise SystemExit('Diagnostic budgets must be 1–10 local poses and 1–100 conformers')
base=Path(__file__).resolve().parent;run=a.run.resolve();run.mkdir(exist_ok=False,parents=True)
for f in a.input.iterdir():
 if f.is_file():shutil.copy2(f,run/f.name)
os.chdir(run)
os.environ.update(PROTAC_PROSETTAC=str(base/'runtime/PRosettaC'),PROTAC_ATTEMPT_ROOT=str(run),PROTAC_SEED=str(a.seed),PROTAC_WORKERS='1',PROTAC_DIAGNOSTIC='1',PROTAC_DIAGNOSTIC_CONFORMERS=str(a.conformers))
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
utils.patchdock=native_frame
original_local=adapter.local_docking
import rosetta
rosetta.local_docking=lambda struct,ca,cb,pa,pb,nstruct=50:original_local(struct,ca,cb,pa,pb,a.local_poses)
intent={'status':'RUNNING','scope':'native-frame integration diagnostic; not an independent global docking benchmark',
 'requested_local_poses':a.local_poses,'conformers_per_local_pose':a.conformers,
 'seed':a.seed,'started_unix':time.time(),'baseline_runtime_manifest_sha256':hashlib.sha256((base/'runtime_manifest.json').read_bytes()).hexdigest()}
(run/'diagnostic.json').write_text(json.dumps(intent,indent=2)+'\n')
shutil.copy2(base/'runtime_manifest.json',run/'runtime_manifest.json')
try:
    main.main('baseline_native_diagnostic',['Protac_params.txt'])
    sampling=sorted((run/'Patchdock_Results').glob('confs_*.sdf.sampling.json'))
    local=sorted((run/'Patchdock_Results').glob('pd.*_docking_????.pdb'))
    if len(sampling)!=len(local) or not local:raise RuntimeError('Missing per-local-pose linker accounting')
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
