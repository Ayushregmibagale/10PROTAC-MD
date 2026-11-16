"""Stage six technical-control inputs with immutable source receipts; not release.

Does not run docking, submit jobs, alter historical inputs or invent biological
labels. Output is explicitly held until engine and per-input qualification pass.
"""
import argparse,hashlib,json,shutil,datetime
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser(description=__doc__)
 for k in ('project','study','output'):p.add_argument('--'+k,type=Path,required=True)
 a=p.parse_args()
 for k in ('project','study','output'):setattr(a,k,getattr(a,k).resolve())
 a.output.mkdir(parents=True,exist_ok=False);rows=[]
 controls=[('GSK215','PTK2_Q05397','native ternary reference / technical control'),('SJF8240','MET_P08581','parent-head technical control'),('SJF1521','EGFR_P00533','modeled conjugated-head technical control'),('JNJ1013','IRAK1_P51617','modeled conjugated-head technical control'),('SJF_delta','MET_P08581','parent-head technical control only; NOT a claimed degradation positive')]
 for drug,target,role in controls:
  source=a.project/'prepared_inputs'/drug/target;dest=a.output/'inputs'/drug/target;dest.mkdir(parents=True)
  names=['full.smi','vhl_head.sdf','poi_head.sdf','protein.pdb','vbc.pdb','Protac_params.txt','selection.json','preparation.json']
  provenance={}
  for name in names:
   f=source/name
   if not f.is_file():raise ValueError('Missing source: '+str(f))
   shutil.copy2(f,dest/name);provenance[name]={'source':str(f),'sha256':sha(f)}
  rows.append({'compound':drug,'target':target,'role':role,'directory':str(dest.relative_to(a.output)),'sources':provenance,'production_ready':False})
 drug='MS143';target='AKT1_P31749';dest=a.output/'inputs'/drug/target;dest.mkdir(parents=True)
 source=a.study/'MS143/strict_review_v2';mapping={'full.smi':a.study.parent/'molecules/MS143/full.smi','vhl_head.sdf':source/'vhl_head.sdf','poi_head.sdf':source/'poi_modeled_head.sdf','protein.pdb':source/'AKT1_parent_protein.pdb','vbc.pdb':a.project/'prepared_inputs/SJF8240/MET_P08581/vbc.pdb','protein_modification_spec.json':a.study/'MS143/protein_modification_spec.json','preparation.json':source/'PREPARATION.json'}
 prov={}
 for name,f in mapping.items():
  if not f.is_file():raise ValueError('Missing MS143 source: '+str(f))
  shutil.copy2(f,dest/name);prov[name]={'source':str(f),'sha256':sha(f)}
 prep=json.loads((source/'PREPARATION.json').read_text());anchors=[prep['vhl_anchor_one_based'],prep['poi_anchor_one_based']]
 (dest/'Protac_params.txt').write_text('Structures: vbc.pdb protein.pdb\nChains: ABC D\nHeads: vhl_head.sdf poi_head.sdf\nAnchor atoms: '+' '.join(map(str,anchors))+'\nProtac: full.smi\nFull: True\nClusterName: LOCAL\n')
 selection=json.loads((a.project/'prepared_inputs/SJF8240/MET_P08581/selection.json').read_text());(dest/'selection.json').write_text(json.dumps(selection,indent=2)+'\n')
 rows.append({'compound':drug,'target':target,'role':'modeled AKT technical control, phosphoThr308 retained','directory':str(dest.relative_to(a.output)),'sources':prov,'production_ready':False,'extra_blocker':'Explicit PTM-aware geometry selection and end-to-end docking preservation required; ATOM-only POI selection is insufficient'})
 for row in rows:
  dest=a.output/row['directory'];m=Chem.MolFromSmiles((dest/'full.smi').read_text().strip());registered=Chem.MolFromSmiles((a.study.parent/'molecules'/row['compound']/'full.smi').read_text().strip())
  if m is None or registered is None or Chem.MolToInchiKey(m)!=Chem.MolToInchiKey(registered):raise ValueError('Supplier identity mismatch: '+row['compound'])
  head_counts={}
  for name in ['vhl_head.sdf','poi_head.sdf']:
   h=Chem.SDMolSupplier(str(dest/name),removeHs=True)[0]
   if h is None or h.GetNumConformers()!=1:raise ValueError('Invalid bound head')
   head_counts[name]=h.GetNumAtoms()
  row.update(full_molecule_inchikey=Chem.MolToInchiKey(m),formula=rdMolDescriptors.CalcMolFormula(m),head_heavy_atoms=head_counts,input_file_sha256={f.name:sha(f) for f in dest.iterdir() if f.is_file()},qualification_status='STAGED_NOT_QUALIFIED')
 report={'status':'SIX_COMPOUND_TECHNICAL_INPUTS_STAGED_ON_HOLD','created_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'MD':False,'submission_performed':False,'production_ready':False,'scope':'Six chemistry/resource technical controls; not full kinase cohort, biological validation or proof of binding','global_blockers':['Docking reference qualification incomplete','Input pocket sterics/state review pending','MS143 full pipeline modification selection/preservation pending','Sampling/resources to be frozen after qualified pilot'],'controls':rows,'source_script_sha256':sha(Path(__file__))}
 (a.output/'INPUT_MANIFEST.json').write_text(json.dumps(report,indent=2)+'\n');(a.output/'DO_NOT_SUBMIT.json').write_text(json.dumps({'production_ready':False,'reason':'Unqualified staged technical inputs; read INPUT_MANIFEST.json'},indent=2)+'\n')
 print(json.dumps({'status':report['status'],'compounds':[r['compound'] for r in rows],'output':str(a.output)},indent=2))
if __name__=='__main__':main()
