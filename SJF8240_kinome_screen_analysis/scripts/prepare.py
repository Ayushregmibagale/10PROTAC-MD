from pathlib import Path
import os,re,json,hashlib
from collections import defaultdict
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];D=Path('/mnt/d/kinome_results');K=Path('/home/ayushregmibagale/kinome_scan');CAN=Path('/home/ayushregmibagale/10PROTAC-MD/SJF8240_reanalysis')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
meta=pd.read_csv(K/'data/kinase_list.tsv',sep='\t');genes=dict(zip(meta.uniprot_id,meta.gene));refs=json.loads((CAN/'config/role_reference.json').read_text());labels=json.loads((CAN/'config/current_md_labels.json').read_text())['labels'];known={refs['targets'][target]:label for target,label in labels.items()};inventory=[];pools=[]
candidates=[(p.name,'initial_clustered',p) for p in sorted(D.iterdir()) if p.is_dir() and re.fullmatch('[A-Z0-9]{6,10}',p.name)]
for name,uid in [('sik3_full','Q9Y2K2'),('rps6kb2_full','Q9UBS0')]:candidates.append((uid,name,D/name))
for uid in ['Q13164','Q6DT37','Q9UBS0','Q9Y5S2']:
 physical=D/'top5_full'/uid
 if not (physical/'Results').exists():physical=K/'data/top5_full'/uid
 candidates.append((uid,'top5_full',physical))
for uid,pool,directory in candidates:
 poses=[];results=directory/'Results'
 if results.exists():
  for entry in os.scandir(results):
   if entry.is_dir():
    for f in os.scandir(entry.path):
     if f.is_file() and f.name.startswith('combined_') and f.name.endswith('.pdb'):poses.append(f.path)
   elif entry.is_file() and entry.name.startswith('combined_') and entry.name.endswith('.pdb'):poses.append(entry.path)
 # Results root and cluster directories can contain exact copies of the same pose.
 aliases=defaultdict(list)
 for name in poses:aliases[Path(name).name].append(name)
 unique=[];duplicate_copies=0
 for name,paths in aliases.items():
  if len(paths)>1:
   digests={sha(p) for p in paths}
   if len(digests)!=1:raise ValueError('Conflicting duplicate pose name '+uid+' '+name)
   duplicate_copies+=len(paths)-1
  unique.append(sorted(paths,key=lambda p:(0 if Path(p).parent.name.startswith('cluster') else 1,p))[0])
 poses=unique
 inventory.append({'duplicate_copies_removed':duplicate_copies,'uid':uid,'gene':genes.get(uid,uid),'pool':pool,'directory':str(directory),'clustered_pose_count':len(poses),'available':bool(poses)})
 if poses:pools.append({'uid':uid,'gene':genes.get(uid,uid),'pool':pool,'directory':str(directory),'poses':sorted(poses),'key':uid+'__'+pool})
protocol={'compound':'SJF8240','contact_A':5.,'pocket_definition_A':5.,'severe_clash_A':1.5,'min_ligand_contact_atoms':3,'anchor_rmsd_A':1.,'vhl_backbone_rmsd_A':1.,'full_e3_reference':'/home/ayushregmibagale/protac-md/random_screen/DDR2/pose001/structure/complex_protein.pdb','role_reference':str(CAN/'config/role_reference.json'),'kinase_structure_root':str(K/'data/structures'),'known_labels':known,'ranking':'Eligible dual pocket-proximity, fit and severe-clash screens first; balanced ligand-contact count descending; maximum pocket distance ascending; deterministic path tie-break. No learned degradation probability. Separate all-poses and one-per-kinase lists; prospective list excludes8 current-labelled targets.','scope':'Only clustered Results/combined*.pdb supplied by available runs. Unclustered Patchdock_Results are not pooled with refined clustered representatives. All incomplete/missing runs recorded; not zeros. Full-VHL rigid reference fit explicitly model-based and not a new docking/MD or proof of productive binding.','pocket_caveat':'Pocket regions derived from supplied PT0/PT1 reference coordinates; proximity of any PROTAC heavy atom, not atom-mapped verification of the two cognate heads.','post_score_tuning':False,'new_MD_or_GNN':False}
protocol['common_input_sha256']={str(p):sha(p) for p in [Path(protocol['full_e3_reference']),Path(protocol['role_reference']),K/'data/kinase_list.tsv',CAN/'config/current_md_labels.json']}
(ROOT/'config/protocol.json').write_text(json.dumps(protocol,indent=2)+'\n');(ROOT/'config/pools.json').write_text(json.dumps(pools,indent=2)+'\n');pd.DataFrame(inventory).to_csv(ROOT/'results/inventory.tsv',sep='\t',index=False);print('Accession dirs:',sum(r['pool']=='initial_clustered' for r in inventory),'available pools:',len(pools),'clustered PDBs:',sum(len(p['poses']) for p in pools));print('Missing/empty:',[r['uid']+' '+r['pool'] for r in inventory if not r['available']])
