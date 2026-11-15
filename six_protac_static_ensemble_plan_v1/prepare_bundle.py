from pathlib import Path
import shutil,json,hashlib,tarfile,argparse
p=argparse.ArgumentParser(description='Build a new immutable private Ernie pilot bundle; refuses existing output')
p.add_argument('--project',type=Path,default=Path('/home/ayushregmibagale/10PROTAC-MD'))
p.add_argument('--inputs',type=Path,required=True,help='Frozen <reference> input directories')
p.add_argument('--output',type=Path,required=True,help='New bundle directory')
a=p.parse_args()
src=a.project/'prosettac_robust/workflow'
base=Path(__file__).resolve().parent;bundle=a.output.resolve()
bundle.mkdir(exist_ok=False)
def copytree(s,d):shutil.copytree(s,d,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
copytree(src/'baseline_engine',bundle/'workflow/baseline_engine')
copytree(src/'scripts/lib',bundle/'workflow/scripts/lib')
(bundle/'workflow/scripts/01_setup').mkdir(parents=True)
for name in ('02_batch_prosettac.py','03_record_attempt.py'):
 shutil.copy2(src/'scripts/01_setup'/name,bundle/'workflow/scripts/01_setup'/name)
for ref in ('5T35_MZ1','7PI4_GSK215'):
 copytree(a.inputs/ref,bundle/'inputs'/ref)
 # Evaluation structures are unnecessary for generation; do not ship them as inputs.
 (bundle/'inputs'/ref/'native_complex.pdb').unlink(missing_ok=True)
shutil.copy2(base/'scripts/03_run_ernie_reference_pilot.py',bundle/'run_pilot.py')
cfg={'purpose':'bounded reference/resource benchmark; no kinome production and no MD',
 'cpu_limit':4,'poll_seconds':15,'max_group_RSS_bytes':24*1024**3,'max_output_bytes':10*1024**3,
 'min_free_bytes':50*1024**3,'tasks':[]}
cfg['tasks'].append({'name':'native_7PI4_GSK215_seed42','reference':'7PI4_GSK215','seed':42,'protocol':'native_frame',
 'n_dock':1,'n_local':3,'conformers':20,'max_global_trials':20000,'timeout_seconds':1800})
for ref in ('7PI4_GSK215','5T35_MZ1'):
 cfg['tasks'].append({'name':f'global_{ref}_seed42','reference':ref,'seed':42,'protocol':'global',
  'n_dock':50,'n_local':10,'conformers':20,'max_global_trials':20000,'timeout_seconds':8*3600})
(bundle/'PILOT_CONFIG.json').write_text(json.dumps(cfg,indent=2)+'\n')
manifest={'source_project':'10PROTAC-MD','source_commit':'f5f3efc','licensed_binaries_bundled':False,
 'evaluation_native_coordinates_in_inputs':False,
 'files':{str(p.relative_to(bundle)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(bundle.rglob('*')) if p.is_file()}}
(bundle/'BUNDLE_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
archive=bundle.with_suffix('.tar.gz')
with tarfile.open(archive,'w:gz') as tar:tar.add(bundle,arcname=bundle.name)
print(json.dumps({'archive':str(archive),'bytes':archive.stat().st_size,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
 'files':len(manifest['files']),'max_workers':4,'MD':False},indent=2))
