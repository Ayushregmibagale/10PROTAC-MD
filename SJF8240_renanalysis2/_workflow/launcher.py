"""Ordered launchers: explicit paths, reviewed backends, no implicit model fitting."""
from pathlib import Path
import argparse,json,subprocess,sys,os,hashlib,ast
ROOT=Path(__file__).resolve().parents[1];GNN=ROOT/'_workflow/gnn';E2=ROOT/'_workflow/e2_full_frame';FEATURES=ROOT/'_workflow/lysine_features'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def write_json(p,value):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_name(p.name+'.pending');tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');os.replace(tmp,p)
def verify_sources():
 expected=json.loads((ROOT/'config/software_manifest.json').read_text())
 for rel,digest in expected.items():require((ROOT/rel).is_file() and sha(ROOT/rel)==digest,'Packaged source/config changed: '+rel)
 return len(expected)
def source_project():
 c=json.loads((ROOT/'config/source_project.json').read_text());source=Path(c['source_project']);require(source.is_dir(),'Original data/project folder is required: '+str(source));return source
def call(script,args,cpu=False):
 verify_sources();env=os.environ.copy()
 if cpu:env.pop('PYTHONPATH',None)
 command=[sys.executable,str(script),*map(str,args)];result=subprocess.run(command,env=env)
 require(result.returncode==0,'Step failed: '+str(script)+' (exit '+str(result.returncode)+')')
def check_original_bindings():
 source=source_project();bindings=json.loads((ROOT/'config/original_input_bindings.json').read_text())
 for rel,digest in bindings.items():require(sha(source/rel)==digest,'Original manifest/label source changed: '+rel)
 for rel,key in [('corrected_md_only_v1/config/input_manifest.json','entries'),('mechanistic_v3/config/temporal_input_manifest.json','graphs')]:
  base=source/rel.split('/')[0];manifest=json.loads((source/rel).read_text())
  for row in manifest[key]:require(sha(base/row['path'])==row['sha256'],'Original graph drift: '+row['path'])
def e2_cached_preflight():
 for p in (E2/'analysis/poses').glob('*/*/receipt.json'):
  record=json.loads(p.read_text());require(sha(p.parent/'packet.npz')==record['packet_sha256'],'Cached E2 packet changed')
  for name,digest in record['input_sha256'].items():require(sha(Path(name))==digest,'Cached E2 source changed: '+name)
def main(number):
 step={1:1,2:2,3:12,4:13,5:3,6:4,7:5,8:6,9:7,10:8,11:9,12:10,13:11}[number]
 descriptions={1:'Prepare SJF8240 input configuration and verify original graph bindings.',2:'Optionally rebuild MD contact graphs from retained trajectories.',3:'Choose existing or rebuilt graphs and prepare a GNN study.',4:'Validate every prepared GNN study before fitting.',5:'Train grouped cross-validation GNNs and save every fold.',6:'Independently replay the GNN fold weights on CPU.',7:'Compute the full201-frame per-lysine E2 geometry.',8:'Compare E2 conventions and compute protein-level statistics.',9:'Verify full-frame E2 geometry and rehash all inputs.',10:'Collect completed results into one summary.',11:'Check this ordered workflow source and configuration.',12:'Optionally extract the SJF8240 lysine/SASA geometry at41 frames per pose.',13:'Build29-feature four-window GNN graphs from extracted lysine packets.'}
 p=argparse.ArgumentParser(description=descriptions[step]);p.add_argument('--dry-run',action='store_true',help='Show the selected paths without running the analysis.')
 if step==2:p.add_argument('--manifest',type=Path,default=ROOT/'config/trajectories.json');p.add_argument('--output',type=Path,default=ROOT/'outputs/md_contact_graphs')
 if step==3:
  p.add_argument('--model',choices=['original_md','lysine_temporal'],default='original_md');p.add_argument('--graph-source',choices=['existing','rebuilt'],default='existing')
 if step in [5,6]:
  p.add_argument('--model',choices=['original_md','lysine_temporal'],default='original_md');p.add_argument('--role',choices=['vhl','elob','both'],default='both')
 if step==5:p.add_argument('--device',choices=['cpu','cuda'],default='cpu')
 if step in [7,12]:p.add_argument('--workers',type=int,default=2)
 a=p.parse_args();require(sys.flags.optimize==0,"Use normal Python; -O disables validation in the preserved feature builder");verify_sources()
 if a.dry_run:print(json.dumps({'step':number,'description':descriptions[step],'root':str(ROOT),'options':{k:str(v) for k,v in vars(a).items()}}));return
 for group,names in [('e2_full_frame',['analysis','results','validation','logs']),('lysine_features',['analysis','inputs','results'])]:
  for name in names:(ROOT/'outputs'/group/name).mkdir(parents=True,exist_ok=True)
 if step==1:
  check_original_bindings();source=source_project();data=json.loads((source/'folder_audit_20261112/real_extraction_cases.json').read_text());write_json(ROOT/'config/trajectories.json',data);write_json(ROOT/'outputs/preparation.json',{'status':'PREPARED','compound_id':'SJF8240','cases':len(data['cases']),'labels':'six current-labelled positives; ABL1/AXL negative; SIK3/RPS6KB2 unknown','source_project':str(source)});print('Inputs prepared; no trajectory extraction or GNN fitting run.')
 elif step==2:
  require(a.manifest.is_file(),'Run01 first or supply --manifest');call(GNN/'build_dynamic.py',['--manifest',a.manifest,'--output',a.output])
 elif step==3:
  check_original_bindings();staging=ROOT/'outputs/prepared_input_manifests';call(GNN/'make_sjf8240_examples.py',['--project-root',source_project(),'--output',staging]);c=json.loads((staging/(a.model+'.json')).read_text())
  if a.graph_source=='rebuilt':
   folder=ROOT/'outputs/md_contact_graphs' if a.model=='original_md' else FEATURES
   manifest=None
   if a.model=='original_md':
    receipt=json.loads((folder/'completion.json').read_text());require(receipt['status']=='COMPLETE','Incomplete contact graph rebuild')
   else:
    manifest=json.loads((FEATURES/'config/temporal_input_manifest.json').read_text())
   c['targets']=[r for r in c['targets'] if r['label'] is not None]
   for row in c['targets']:
    for role in c['roles']:
     if a.model=='original_md':
      graph=folder/'targets'/row['target_id']/role/'graph.npz';r=json.loads(graph.with_suffix('.json').read_text());require(sha(graph)==r['graph_sha256'] and r['compound_id']=='SJF8240' and r['role']==role and r['target_id']==row['target_id'],'Rebuilt graph receipt mismatch')
     else:
      match=[r for r in manifest['graphs'] if r['target']==row['target_id'] and r['role']==role];require(len(match)==1,'Missing/duplicate temporal graph');graph=FEATURES/match[0]['path'];require(sha(graph)==match[0]['sha256'],'Temporal graph hash mismatch')
     row['graphs'][role].update(path=str(graph.resolve()),sha256=sha(graph))
  c['clean_folder_graph_source']=a.graph_source;path=ROOT/'config/studies'/(a.model+'.json');write_json(path,c);call(GNN/'validate_study.py',['--study',path]);print('Prepared',path)
 elif step==4:
  paths=sorted((ROOT/'config/studies').glob('*.json'));require(paths,'Run05 for each desired model first.');batch={'expected_compound_ids':['SJF8240'],'studies':[{'path':str(path),'sha256':sha(path)} for path in paths]};write_json(ROOT/'config/batch.json',batch);call(GNN/'validate_batch.py',['--batch',ROOT/'config/batch.json','--output',ROOT/'outputs/preflight.json'])
 elif step in [5,6]:
  study=ROOT/'config/studies'/(a.model+'.json');require(study.exists(),'Run05 then06 first');call(GNN/'validate_study.py',['--study',study]);roles=['vhl','elob'] if a.role=='both' else [a.role]
  for role in roles:
   args=['--study',study,'--output',ROOT/'outputs/gnn'/a.model/role,'--role',role,'--device',a.device if step==5 else 'cpu']
   if step==6:args+=['--verify']
   call(GNN/'run_study.py',args,cpu=step==6)
 elif step==7:
  require(a.workers>0,'Positive worker count required');e2_cached_preflight();call(E2/'scripts/extract.py',['--workers',a.workers])
 elif step==8:
  require((E2/'analysis/completion.json').exists(),'Complete09 first');e2_cached_preflight();call(E2/'scripts/analyze.py',[])
 elif step==9:
  require((E2/'analysis/completion.json').exists(),'Complete09 first');call(E2/'scripts/verify_replays.py',[]);call(E2/'scripts/check_inputs.py',[])
 elif step==10:
  lines=['# New results in SJF8240_renanalysis2','', 'Only runs made in this folder are listed below. Original findings are in reference_results/ and are not presented as new runs.',''];found=0
  for path in sorted((ROOT/'outputs/gnn').glob('*/*/completion.json')):
   c=json.loads(path.read_text());v=path.parent/'cpu_verification.json';lines+=['## GNN: '+str(path.parent.relative_to(ROOT/'outputs/gnn')),'',json.dumps(c['metrics'],indent=2),'','CPU verification: '+('passed' if v.exists() and json.loads(v.read_text())['status']=='CPU_VERIFIED' else 'pending'),''];found+=1
  e2=E2/'RESULTS.md'
  if e2.exists():
   require((E2/'validation/independent_unwrap_replay.json').exists() and (E2/'validation/current_input_recheck.json').exists(),'Complete11 before exporting the E2 report');lines+=['## Full-frame E2 result','',e2.read_text()];found+=1
  if not found:lines+=['No new analyses have been run. See reference_results/ for preserved earlier findings.']
  dest=ROOT/'outputs/SUMMARY.md';dest.write_text('\n'.join(lines)+'\n');print(dest)
 elif step==12:
  check_original_bindings();require(a.workers>0,'Positive worker count required')
  for p in (FEATURES/'analysis/lysines/poses').glob('*/*/receipt.json'):
   record=json.loads(p.read_text());require(sha(p.parent/'packet.npz')==record['packet_sha256'],'Cached lysine packet changed')
   for name,digest in record['source_sha256'].items():require(sha(Path(name))==digest,'Cached lysine source changed')
  call(FEATURES/'scripts/extract_lysines.py',['--workers',a.workers])
 elif step==13:
  check_original_bindings();require((FEATURES/'analysis/lysines/completion.json').exists(),'Complete03 before04');call(FEATURES/'scripts/build_temporal.py',[])
 elif step==11:
  provenance=json.loads((ROOT/'config/backend_provenance.json').read_text())
  for rel,row in provenance.items():require(sha(ROOT/rel)==row['sha256'],'Copied backend differs from reviewed source')
  for path in ROOT.rglob('*.py'):
   if '__pycache__' not in path.parts:ast.parse(path.read_text())
  # Frozen E2 source manifests remain valid after copying and routing output directories.
  for name in ['execution_manifest.json','reporting_manifest.json']:
   for rel,digest in json.loads((E2/'config'/name).read_text()).items():require(sha(E2/rel)==digest,'Frozen E2 backend drift '+rel)
  write_json(ROOT/'outputs/workflow_check.json',{'status':'PASS','reviewed_backend_files':len(provenance),'immutable_source_files':verify_sources(),'original_scientific_results_modified':False});print('Ordered workflow source checks passed.')
if __name__=='__main__':raise SystemExit('Use the numbered launchers in the parent directory.')
