"""Inventory exact retained source structures; no original writes or execution."""
from pathlib import Path
import ast,json,os,hashlib
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).resolve().parents[1]
ORIGINAL=Path('/home/ayushregmibagale/protac-md')
AUDIT=Path('/home/ayushregmibagale/10PROTAC-MD/prosettac_robust')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def inventory(item):
 poi,cfg=item
 folder=Path(cfg['path']) if 'path' in cfg else Path('/mnt/d/PROTAC degradation Prediction/pRosettaC/runs')/cfg['base']/poi/cfg.get('src_subdir','Patchdock_Results')
 files=sorted(folder.glob('cluster*/combined_*.pdb') if cfg.get('src_subdir')=='Results' else folder.glob('combined_*.pdb'))
 if not files:raise ValueError('No retained static poses: '+str(folder))
 dest=ROOT/'config/pose_lists';dest.mkdir(exist_ok=True)
 p=dest/(poi+'.txt');p.write_text('\n'.join(map(str,files))+'\n')
 record={'poi':poi,'source_directory':str(folder),'poses':len(files),'first_pose':str(files[0]),'first_pose_sha256':sha(files[0]),'pose_list':str(p.relative_to(ROOT)),'pose_list_sha256':sha(p),'kin_chain':cfg['kin_chain'],'kin_min':cfg['kin_min'],'historical_label':None if cfg['degrader'] is None else int(cfg['degrader']),'training_inclusion':poi not in ['CDK17','SIK3_full','RPS6KB2_full'],'structure_identity_exception':'Retained CDK17-labelled structure is CDK14/O94921; excluded from main evaluation' if poi=='CDK17' else None}
 print(json.dumps({'poi':poi,'poses':len(files),'folder':str(folder)}),flush=True);return record
if __name__=='__main__':
 source=ORIGINAL/'scripts/02_analysis/04a_build_ensemble_graphs.py';tree=ast.parse(source.read_text());node=next(x.value for x in tree.body if isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='POIS_CONFIG' for t in x.targets))
 cfg=eval(compile(ast.Expression(node),str(source),'eval'),{'Path':Path,'__builtins__':{}})
 labels=json.loads((AUDIT/'workflow/config/sjf8240_publication_cohort.json').read_text())['historical_labels']
 selected=[(p,cfg[p]) for p in list(labels)+['SIK3_full','RPS6KB2_full']]
 with ThreadPoolExecutor(max_workers=2) as pool:records=list(pool.map(inventory,selected))
 (ROOT/'config/source_inventory.json').write_text(json.dumps({'status':'SOURCE_INVENTORY_COMPLETE','records':records,'original_config_sha256':sha(source),'new_MD_or_docking':False},indent=2)+'\n')
