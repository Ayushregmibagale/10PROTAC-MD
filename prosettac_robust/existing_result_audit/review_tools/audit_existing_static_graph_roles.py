"""Read one retained pose per configured core target; never execute old scripts."""
from pathlib import Path
import ast,json,sys,importlib.util
from Bio.SeqUtils import seq1
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import pdb_atoms,sha256,atomic_json
spec=importlib.util.spec_from_file_location('static_role_identity',Path(__file__).with_name('audit_existing_protein_roles.py'));roles=importlib.util.module_from_spec(spec);spec.loader.exec_module(roles)
original=Path('/home/ayushregmibagale/protac-md/scripts/02_analysis/04a_build_ensemble_graphs.py');tree=ast.parse(original.read_text());node=next(x.value for x in tree.body if isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='POIS_CONFIG' for t in x.targets))
cohort=json.loads((ROOT/'workflow/config/sjf8240_publication_cohort.json').read_text())['protein_groups'];cache=ROOT/'workflow/sequence_provenance/uniprot_batches/sjf8240_reference_cohort.json';entries={r['primaryAccession']:r['sequence']['value'] for r in json.loads(cache.read_text())['results'] if r['primaryAccession'] in {'P40337','Q15370','Q15369'}}
base=Path('/mnt/d/PROTAC degradation Prediction/pRosettaC/runs');records=[]
for key,value in zip(node.keys,node.values):
 poi=ast.literal_eval(key)
 if poi not in cohort:continue
 cfg=ast.literal_eval(value);directory=base/cfg['base']/poi/cfg.get('src_subdir','Patchdock_Results');pose=next(directory.rglob('combined_*.pdb'),None)
 if pose is None:records.append({'poi':poi,'status':'NO_REPRESENTATIVE_FOUND','directory':str(directory)});continue
 before=sha256(pose);atoms=pdb_atoms(pose);ca=[x for x in atoms if x['record']=='ATOM' and x['name']=='CA'];selected=[x for x in ca if x['chain']=='A' and 1<=x['resnum']<=103];sequence=''.join(seq1(x['resname'],undef_code='X') for x in selected)
 try:found=roles.identify(sequence,entries);status='LEGACY_SELECTED_COMPONENT_IDENTIFIED'
 except ValueError as exc:found={'error':str(exc)};status='REVIEW_REQUIRED'
 chains={}
 for chain in sorted({x['chain'] for x in ca}):
  seq=''.join(seq1(x['resname'],undef_code='X') for x in ca if x['chain']==chain)
  try:chains[chain]=roles.identify(seq,entries)
  except ValueError:pass
 assert sha256(pose)==before
 records.append({'poi':poi,'status':status,'representative_selection':'First available configured retained combined PDB; chosen for sequence audit without labels/geometry optimization','pose':str(pose),'pdb_sha256':before,'legacy_static_VHL_CA_count':len(selected),'legacy_static_selected_protein':found,'independently_identified_E3_chains':chains})
out=ROOT/'existing_result_audit/static_role_representatives_v1';out.mkdir(exist_ok=False);atomic_json(out/'records.json',records);summary={'status':'REPRESENTATIVE_STATIC_ROLE_AUDIT','targets_requested':len(cohort),'targets_found':sum(x['status']=='LEGACY_SELECTED_COMPONENT_IDENTIFIED' for x in records),'legacy_selected_UID_counts':{uid:sum(x.get('legacy_static_selected_protein',{}).get('uniprot')==uid for x in records) for uid in entries},'records':len(records),'source_script_sha256':sha256(original),'audit_source_sha256':sha256(Path(__file__)),'sequence_cache_sha256':sha256(cache),'records_sha256':sha256(out/'records.json'),'original_files_modified':False,'new_PRosettaC_or_MD_executed':False,'limitation':'One retained structure per target verifies these representatives; this is not an exhaustive audit of every retained pose.'};atomic_json(out/'summary.json',summary);print(json.dumps(summary,indent=2))
