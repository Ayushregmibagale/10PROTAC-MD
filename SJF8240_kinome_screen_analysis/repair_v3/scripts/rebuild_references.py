"""Regenerate verified bound-parent reference heads; preserve original chemistry/files.
Correct supplier identity is a separate candidate branch, not historical relabeling.
"""
from pathlib import Path
import os,sys,json,importlib.util,concurrent.futures,argparse
sys.dont_write_bytecode=True;os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,pandas as pd
from rdkit import Chem
BASE=Path(os.environ.get('KINOME_ANALYSIS_ROOT',str(Path(__file__).resolve().parents[2])));ROOT=Path(__file__).resolve().parents[1];W=Path('/home/ayushregmibagale/10PROTAC-MD/prosettac_robust/workflow');REF=W/'reference_structures/validated/SJF8240'
sys.path.insert(0,str(BASE/'engagement_v2/scripts'));import check_engagement as c
sys.path.insert(0,str(W/'scripts/lib'))
spec=importlib.util.spec_from_file_location('reference_placement',W/'scripts/01_setup/01_place_warhead_batch.py');placement=importlib.util.module_from_spec(spec);spec.loader.exec_module(placement)
from protac_qc import pdb_atoms

def save_sdf(path,mol,xyz=None):
 m=Chem.Mol(mol)
 if xyz is not None:
  for i,p in enumerate(xyz):m.GetConformer().SetAtomPosition(i,p)
 writer=Chem.SDWriter(str(path));writer.write(m);writer.close()

def validate_head(m):
 xyz=m.GetConformer().GetPositions();lengths=[np.linalg.norm(xyz[b.GetBeginAtomIdx()]-xyz[b.GetEndAtomIdx()]) for b in m.GetBonds()];c.require(min(lengths)>=1. and max(lengths)<=2.1,'Invalid bound-head bond geometry');return {'atoms':m.GetNumAtoms(),'min_bond_A':float(min(lengths)),'max_bond_A':float(max(lengths)),'inchikey':Chem.MolToInchiKey(m)}

def rebuild(pool):
 b=c.load_backend(BASE);cfg=json.loads((BASE/'config/protocol.json').read_text());out=ROOT/'supplier_reference_candidates'/pool['key'];out.mkdir(exist_ok=True,parents=True);c.require(not (out/'receipt.json').exists(),'Reference candidate already exists; use archived new run');result={'pool':pool['key'],'gene':pool['gene'],'uid':pool['uid'],'status':'REJECTED','production_ready':False,'native_ternary':False,'reference_kind':'separately deposited bound-parent heads; transfer into archived kinase model is a prediction','supplier_vs_historical':'Separate supplier-matched SJF8240 candidate inputs; not the molecule used in archived poses/MD','alignment_attempts':[]}
 try:
  context=b.model_context(pool,cfg);sequences=json.loads(Path(cfg['role_reference']).read_text())['sequences'];native=b.parse(REF/'proteins.pdb');native_eloc=next(g for g in native['groups'] if b.role_identity(g,'Q15369',sequences)>=.95);vhl=Chem.MolFromMolFile(str(REF/'vhl_head.sdf'),removeHs=True);poi=Chem.MolFromMolFile(str(REF/'poi_head.sdf'),removeHs=True);result['source_head_geometry']={'vhl':validate_head(vhl),'poi':validate_head(poi)}
  vxyz,rmsd=b.fit(b.ca_xyz(native_eloc),b.ca_xyz(context['roles']['eloc']),vhl.GetConformer().GetPositions());c.require(rmsd<=1.,'ELOC reference alignment failed');save_sdf(out/'vhl_head.sdf',vhl,vxyz);result['vhl_anchor_RMSD_A']=rmsd
  # Transform all full native E3 atoms through that same anchor; chain roles explicit.
  proteinlines=[];chainmap={'A':'C','B':'A','C':'B'}
  for line in (REF/'proteins.pdb').read_text().splitlines():
   if line.startswith('ATOM  '):
    point=np.array([[float(line[30:38]),float(line[38:46]),float(line[46:54])]]);transformed,_=b.fit(b.ca_xyz(native_eloc),b.ca_xyz(context['roles']['eloc']),point);x,y,z=transformed[0];proteinlines.append(line[:21]+chainmap[line[21]]+line[22:30]+f'{x:8.3f}{y:8.3f}{z:8.3f}'+line[54:])
   elif line.startswith('TER'):proteinlines.append('TER')
  (out/'vbc.pdb').write_text('\n'.join(proteinlines)+'\nEND\n');d=Path(pool['directory']);targetsource=d/'Side1_0001_A.pdb'
  if not targetsource.exists():targetsource=d/'protein.pdb'
  targetatoms=[a for a in pdb_atoms(targetsource) if a['record']=='ATOM'];targetca=[a for a in targetatoms if a['name']=='CA'];referenceca=[a for a in pdb_atoms(REF/'protein.pdb') if a['name']=='CA'];refcoords=poi.GetConformer().GetPositions();fit=None
  for method in ['sequence','ce']:
   try:r,t,meta=placement.align_domain(targetca,referenceca,refcoords,method);fit=(r,t,meta);result['alignment_attempts'].append({'method':method,'status':'PASS'});break
   except Exception as exc:result['alignment_attempts'].append({'method':method,'status':'REJECTED','reason':str(exc)})
  c.require(fit is not None,'No trustworthy kinase/domain pocket transfer');r,t,meta=fit;coords=(refcoords-t)@r;save_sdf(out/'poi_head.sdf',poi,coords);result['kinase_alignment']=meta
  (out/'protein.pdb').write_text('\n'.join(a['line'][:21]+'D'+a['line'][22:] for a in targetatoms)+'\nTER\nEND\n');registry=json.loads((W/'config/molecules.json').read_text())['compounds']['SJF8240'];full=Chem.MolFromSmiles(registry['smiles']);c.require(Chem.MolToInchiKey(full)==registry['inchikey'],'Supplier identity mismatch');c.require(full.HasSubstructMatch(vhl,useChirality=True) and full.HasSubstructMatch(poi,useChirality=True),'Head does not match supplier identity');(out/'full.smi').write_text(registry['smiles']+'\n');heads=json.loads((REF/'heads.json').read_text());(out/'Protac_params.txt').write_text('Structures: vbc.pdb protein.pdb\nChains: ABC D\nHeads: vhl_head.sdf poi_head.sdf\nAnchor atoms: '+str(heads['vhl']['anchor_one_based'])+' '+str(heads['poi']['anchor_one_based'])+'\nProtac: full.smi\nFull: True\nClusterName: LOCAL\n')
  tar=np.array([a['xyz'] for a in targetatoms if a['element'] not in ['H','D']]);els=[a['element'].upper() for a in targetatoms if a['element'] not in ['H','D']];headels=[a.GetSymbol().upper() for a in poi.GetAtoms()];ov=c.overlap(coords,headels,tar,els);result['poi_head_overlap']=ov;result['molecule_inchikey']=registry['inchikey'];result['status']='REFERENCE_CANDIDATE' if ov['max_overlap_A']<=.8 else 'REJECTED_HEAD_CLASH';result['reason']='Not approved as bound pose or production input; no ternary docking performed'
 except Exception as exc:result['reason']=str(exc)
 result['source_sha256']={str(p):c.sha(p) for p in [REF/'vhl_head.sdf',REF/'poi_head.sdf',REF/'proteins.pdb',REF/'protein.pdb',REF/'heads.json',W/'config/molecules.json',W/'scripts/01_setup/01_place_warhead_batch.py',Path(__file__)]};result['files_sha256']={p.name:c.sha(p) for p in out.iterdir() if p.is_file() and p.name!='receipt.json'};(out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');return result

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=4);ap.add_argument('--limit',type=int);args=ap.parse_args();receipts=json.loads((BASE/'engagement_v2/head_mapping_receipts.json').read_text());keys={r['pool'] for r in receipts};pools=[p for p in json.loads((BASE/'config/pools.json').read_text()) if p['key'] in keys];pools=pools[:args.limit] if args.limit else pools;results=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as ex:
  for i,r in enumerate(ex.map(rebuild,pools)):results.append(r);print(i+1,len(pools),r['pool'],r['status'],flush=True)
 (ROOT/'reference_rebuild_receipts.json').write_text(json.dumps(results,indent=2)+'\n');pd.DataFrame([{'pool':r['pool'],'gene':r['gene'],'uid':r['uid'],'status':r['status'],'reason':r['reason'],'vhl_reference_written':bool((ROOT/'supplier_reference_candidates'/r['pool']/'vhl_head.sdf').exists())} for r in results]).to_csv(ROOT/'reference_rebuild_summary.tsv',sep='\t',index=False)
if __name__=='__main__':main()
