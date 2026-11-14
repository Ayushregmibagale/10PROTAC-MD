"""Chemical head mapping and heavy-atom overlap audit; no new docking/MD."""
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import sys,json,hashlib,argparse,concurrent.futures,importlib.util
import numpy as np,pandas as pd
from scipy.spatial.distance import cdist
from scipy.spatial import cKDTree
from rdkit import Chem
RADII={'C':1.70,'N':1.55,'O':1.52,'F':1.47,'S':1.80,'P':1.80,'CL':1.75,'BR':1.85,'I':1.98}
PARAM=Path('/home/ayushregmibagale/kinome_scan/data/kinase_pRosettaC/Q16832/Patchdock_Results/PT_1_14.params')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def require(ok,msg):
 if not ok:raise ValueError(msg)
def molecule(path):
 rw=Chem.RWMol();names=[];bonds=[]
 for line in Path(path).read_text().splitlines():
  s=line.split()
  if s and s[0]=='ATOM':
   name=s[1];element=''.join(c for c in name if c.isalpha());rw.AddAtom(Chem.Atom(element));names.append(name)
  elif s and s[0]=='BOND_TYPE':bonds.append(s[1:])
 for a,b,o in bonds:rw.AddBond(names.index(a),names.index(b),{'1':Chem.BondType.SINGLE,'2':Chem.BondType.DOUBLE,'3':Chem.BondType.TRIPLE,'ar':Chem.BondType.AROMATIC}[o])
 m=rw.GetMol();Chem.SanitizeMol(m);return m,names,bonds
def atoms(path,chain=None):
 out=[]
 for l in Path(path).read_text().splitlines():
  if not l.startswith(('ATOM  ','HETATM')) or l[16:17] not in [' ','A']:continue
  if chain is not None and l[21]!=chain:continue
  name=l[12:16].strip();el=l[76:78].strip().upper() or name.lstrip('0123456789')[0]
  if el in ['H','D']:continue
  require(el in RADII,'Unsupported element '+el);out.append((name,el,np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])))
 return out
def bond_audit(names,bonds,records):
 coords={a[0]:a[2] for a in records};require(len(coords)==len(records),'Duplicate ligand atom name');require(set(coords)==set(names),'Ligand atom identity differs from topology')
 lengths=[float(np.linalg.norm(coords[a]-coords[b])) for a,b,_ in bonds]
 return {'bond_min_A':min(lengths),'bond_max_A':max(lengths),'bad_bonds':sum(d<1. or d>2.1 for d in lengths)},np.array([coords[n] for n in names])
def load_backend(root):
 spec=importlib.util.spec_from_file_location('static_backend',root/'scripts/rank_poses.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod

def overlap(a,ea,b,eb):
 ds=cdist(a,b);v=np.array([RADII[x] for x in ea])[:,None]+np.array([RADII[x] for x in eb])[None,:]-ds
 return {'max_overlap_A':float(v.max()),'pairs_gt04':int((v>.4).sum()),'pairs_gt06':int((v>.6).sum()),'pairs_gt08':int((v>.8).sum()),'min_distance_A':float(ds.min())}
def protein_elements(g):return [name.lstrip('0123456789')[0] for name in g['names']]
def head_map(pool,full):
 m,n,b=full;directory=Path(pool['directory']);heads={};evidence={}
 for key,label in [('PT0','VHL'),('PT1','POI')]:
  p=directory/(key+'.params');hm,hn,hb=molecule(p);matches=m.GetSubstructMatches(hm,uniquify=False);sets={tuple(sorted(x)) for x in matches};require(len(sets)==1,'Ambiguous/absent chemical head map '+key)
  ix=next(iter(sets));heads[label]=list(ix)
  check,xyz=bond_audit(hn,hb,atoms(directory/(key+'_0001.pdb')))
  evidence[label]={'atom_count':len(ix),'mapped_full_atom_names':[n[i] for i in ix],'reference_geometry':check,'params_sha256':sha(p),'reference_sha256':sha(directory/(key+'_0001.pdb')),'reference_engagement_verifiable':check['bad_bonds']==0}
 require(not set(heads['VHL'])&set(heads['POI']),'Overlapping head atom maps')
 return heads,evidence

def worker(task):
 root,pool,records=task;root=Path(root);backend=load_backend(root);cfg=json.loads((root/'config/protocol.json').read_text());full=molecule(PARAM);heads,evidence=head_map(pool,full);context=backend.model_context(pool,cfg);refs=json.loads(Path(cfg['role_reference']).read_text())['sequences'];results=[];errors=[]
 for row in records:
  try:
   path=Path(row['path']);require(sha(path)==row['pdb_sha256'],'Raw source hash mismatch');d=backend.parse(path);ligrecords=atoms(path,'X');bondaudit,lig=bond_audit(full[1],full[2],ligrecords);els=[next(a[1] for a in ligrecords if a[0]==name) for name in full[1]]
   g={role:backend.select_by_sequence(d['groups'],r['sequence']) for role,r in context['roles'].items() if role!='vhl'};poi=backend.select_by_sequence(d['groups'],context['target']['sequence']);vhl=context['full_vhl'];vfull,_=backend.fit(context['anchor_reference'],backend.ca_xyz(g['eloc']),vhl['xyz']);vfrag=[x for x in d['groups'] if len(x['ca'])>=20 and backend.role_identity(x,'P40337',refs)>=.95];require(len(vfrag)==1,'VHL identity')
   ev=protein_elements(vhl);ep=protein_elements(poi);e3=np.vstack([vfull,g['elob']['xyz'],g['eloc']['xyz']]);ee=ev+protein_elements(g['elob'])+protein_elements(g['eloc'])
   out={**row,**bondaudit,'head_VHL_atoms':len(heads['VHL']),'head_POI_atoms':len(heads['POI']),'VHL_reference_valid':evidence['VHL']['reference_engagement_verifiable'],'POI_reference_valid':evidence['POI']['reference_engagement_verifiable']}
   for name,a,ea,b,eb in [('ligand_VHL',lig,els,vfull,ev),('ligand_POI',lig,els,poi['xyz'],ep),('ligand_scaffold',lig,els,e3,ee),('E3_POI',e3,ee,poi['xyz'],ep),('observedE3_POI',np.vstack([vfrag[0]['xyz'],g['elob']['xyz'],g['eloc']['xyz']]),protein_elements(vfrag[0])+protein_elements(g['elob'])+protein_elements(g['eloc']),poi['xyz'],ep)]:
    out.update({name+'_'+k:v for k,v in overlap(a,ea,b,eb).items()})
   out['maximum_intercomponent_overlap_A']=max(out['ligand_scaffold_max_overlap_A'],out['ligand_POI_max_overlap_A'],out['E3_POI_max_overlap_A'])
   # Internal nonbonded check excludes directly bonded and 1-3 neighbors.
   graphdist=Chem.GetDistanceMatrix(full[0]);internal=np.array([RADII[e] for e in els])[:,None]+np.array([RADII[e] for e in els])[None,:]-cdist(lig,lig);mask=np.triu(graphdist>2,1);out['ligand_internal_nonbonded_max_overlap_A']=float(internal[mask].max());out['ligand_internal_pairs_gt08']=int((internal[mask]>.8).sum())
   dv=cdist(lig[heads['VHL']],vfull).min(axis=1);dk=cdist(lig[heads['POI']],poi['xyz']).min(axis=1)
   out['correct_VHL_head_min_A']=float(dv.min());out['correct_POI_head_min_A']=float(dk.min())
   for cutoff,tag in [(4.,'4A'),(4.5,'4p5A'),(5.,'5A')]:
    nv=int((dv<cutoff).sum());nk=int((dk<cutoff).sum());out['VHL_head_contacts_'+tag]=nv;out['POI_head_contacts_'+tag]=nk;out['balanced_head_contacts_'+tag]=min(nv,nk);out['heads_contact_'+tag]=bool(min(nv,nk)>=3)
    for tol,label in [(.4,'strict04'),(.6,'moderate06'),(.8,'lenient08')]:out['geometry_pass_'+tag+'_'+label]=bool(min(nv,nk)>=3 and out['maximum_intercomponent_overlap_A']<=tol and out['bad_bonds']==0 and out['ligand_internal_pairs_gt08']==0)
   out['engagement_confirmable']=bool(out['VHL_reference_valid'] and out['POI_reference_valid']);results.append(out)
  except Exception as e:errors.append({'path':row['path'],'error':str(e)})
 return results,errors,{'pool':pool['key'],'heads':evidence,'context_sha256':context['source_sha256']}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--analysis-root',type=Path,default=Path(__file__).resolve().parents[2]);ap.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--workers',type=int,default=6);args=ap.parse_args();root=args.analysis_root;out=args.output;out.mkdir(exist_ok=True,parents=True)
 table=root/'results/all_static_ranked_poses.tsv';df=pd.read_csv(table,sep='\t');df=df[df.passes_full_screen.eq(True)].drop_duplicates('pdb_sha256');pools=json.loads((root/'config/pools.json').read_text());tasks=[]
 for (uid,poolname),frame in df.groupby(['uid','pool']):tasks.append((str(root),next(p for p in pools if p['uid']==uid and p['pool']==poolname),frame.to_dict('records')))
 protocol={'scope':'All 723 unique original QC-passing poses; compare previous distinct-kinase top10. Not a fresh search of rejected pools or discarded poses.','cutoffs_A':[4,4.5,5],'minimum_correct_head_contact_atoms':3,'heavy_atom_vdw_radii_A':RADII,'overlap_thresholds_A':[.4,.6,.8],'primary_diagnostic':'4A correct-head contacts plus maximum intercomponent heavy-atom overlap <=0.6A; ligand covalent geometry 1.0-2.1A and no internal nonbonded overlaps >0.8A excluding 1-2/1-3 pairs. Operational thresholds, not MolProbity or binding validation.','reference_policy':'Head chemical mappings must be unique atom sets; invalid reference-head bond lengths prevent validated pocket-engagement/RMSD claims. Correct-head contacts alone do not establish pocket engagement.','source_sha256':{'original_table':sha(table),'full_ligand_params':sha(PARAM),'script':sha(__file__),'static_backend':sha(root/'scripts/rank_poses.py'),'config':sha(root/'config/protocol.json'),'pool_inventory':sha(root/'config/pools.json')}}
 (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n');rows=[];errors=[];contexts=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as executor:
  for i,(rr,ee,cc) in enumerate(executor.map(worker,tasks)):
   rows+=rr;errors+=ee;contexts.append(cc);print(f'{i+1}/{len(tasks)} pools; {len(rows)} poses checked; {len(errors)} errors',flush=True)
 d=pd.DataFrame(rows);d.to_csv(out/'pose_checks.tsv',sep='\t',index=False);(out/'rejections.json').write_text(json.dumps(errors,indent=2)+'\n');(out/'head_mapping_receipts.json').write_text(json.dumps(contexts,indent=2)+'\n');require(len(d)+len(errors)==len(df),'Pose accounting mismatch')
 old=pd.read_csv(root/'results/top10_static_distinct_prospective_kinases.tsv',sep='\t');comparison=d[d.path.isin(old.path)].copy();comparison['previous_rank']=comparison.path.map(dict(zip(old.path,range(1,len(old)+1))));comparison.sort_values('previous_rank').to_csv(out/'previous_top10_comparison.tsv',sep='\t',index=False)
 stats=[]
 for cutoff,tag in [(4.,'4A'),(4.5,'4p5A'),(5.,'5A')]:
  for tol,label in [(.4,'strict04'),(.6,'moderate06'),(.8,'lenient08')]:
   field='geometry_pass_'+tag+'_'+label;eligible=d[d[field]&d.source_label.astype(str).eq('unknown')].copy();eligible=eligible.sort_values(['balanced_head_contacts_'+tag,'maximum_intercomponent_overlap_A','path'],ascending=[False,True,True]).drop_duplicates('uid');eligible.head(10).to_csv(out/f'top10_geometry_{tag}_{label}.tsv',sep='\t',index=False);stats.append({'cutoff_A':cutoff,'maximum_overlap_A':tol,'passing_poses_all':int(d[field].sum()),'prospective_kinases':len(eligible),'previous_top10_passing':int(comparison[field].sum())})
 pd.DataFrame(stats).to_csv(out/'cutoff_sensitivity.tsv',sep='\t',index=False)
 receipt={'complete':True,'expected_poses':len(df),'checked_poses':len(d),'rejected_poses':len(errors),'source_sha256':protocol['source_sha256'],'outputs_sha256':{p.name:sha(p) for p in out.glob('*.tsv')}};(out/'completion.json').write_text(json.dumps(receipt,indent=2)+'\n');print(pd.DataFrame(stats).to_string(index=False),flush=True)
if __name__=='__main__':main()
