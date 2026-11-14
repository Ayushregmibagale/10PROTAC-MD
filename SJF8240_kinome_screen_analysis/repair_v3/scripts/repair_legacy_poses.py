"""Bounded local geometry repairs of original top10; preserves archived graph/3D stereo.
The historical molecule is NOT supplier SJF8240. No docking or MD is run.
"""
from pathlib import Path
import os,sys,json,importlib.util,hashlib,argparse,concurrent.futures,time
os.environ['OPENBLAS_NUM_THREADS']='1';sys.dont_write_bytecode=True
import numpy as np,pandas as pd
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from rdkit import Chem
from rdkit.Chem import AllChem
BASE=Path(os.environ.get('KINOME_ANALYSIS_ROOT',str(Path(__file__).resolve().parents[2])))
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE/'engagement_v2/scripts'));import check_engagement as c

def make_complex(path,context,backend,cfg):
 d=backend.parse(path);g={role:backend.select_by_sequence(d['groups'],r['sequence']) for role,r in context['roles'].items() if role!='vhl'};poi=backend.select_by_sequence(d['groups'],context['target']['sequence']);vhl=context['full_vhl'];vfull,_=backend.fit(context['anchor_reference'],backend.ca_xyz(g['eloc']),vhl['xyz']);prot=np.vstack([vfull,g['elob']['xyz'],g['eloc']['xyz'],poi['xyz']]);el=c.protein_elements(vhl)+c.protein_elements(g['elob'])+c.protein_elements(g['eloc'])+c.protein_elements(poi)
 return d,g,poi,vhl,vfull,prot,el

def process(task):
 row,pool,attempts=task;start=time.monotonic();backend=c.load_backend(BASE);cfg=json.loads((BASE/'config/protocol.json').read_text());context=backend.model_context(pool,cfg);d,g,poi,vhl,vfull,protein,pels=make_complex(Path(row['path']),context,backend,cfg);c.require(c.sha(row['path'])==row['pdb_sha256'],'Source drift');params,names,bonds=c.molecule(c.PARAM);_,xyz=c.bond_audit(names,bonds,c.atoms(row['path'],'X'))
 # Derive stereo from actual archived coordinates; do not silently substitute intended/supplier stereo.
 mol=Chem.Mol(params);conf=Chem.Conformer(len(names))
 for i,v in enumerate(xyz):conf.SetAtomPosition(i,v)
 mol.AddConformer(conf);Chem.AssignAtomChiralTagsFromStructure(mol,replaceExistingTags=True);Chem.AssignStereochemistry(mol,cleanIt=True,force=True);original_key=Chem.MolToInchiKey(mol);nh=len(names);els=[a.GetSymbol().upper() for a in mol.GetAtoms()];radii=np.array([c.RADII[e] for e in els]);graphdist=Chem.GetDistanceMatrix(mol);i_internal,j_internal=np.where(np.triu(graphdist>2,1));rs_internal=radii[i_internal]+radii[j_internal]
 hmol=Chem.AddHs(mol,addCoords=True);props=AllChem.MMFFGetMoleculeProperties(hmol,mmffVariant='MMFF94s');c.require(props is not None,'MMFF parameters unavailable');ff=AllChem.MMFFGetMoleculeForceField(hmol,props,nonBondedThresh=100.);c.require(ff is not None,'MMFF unavailable');init=hmol.GetConformer().GetPositions();nall=len(init);ptree=cKDTree(protein);prad=np.array([c.RADII[e] for e in pels]);polar=np.isin(els,['N','O'])[:,None];headmaps,evidence=c.head_map(pool,(params,names,bonds));headids=set(headmaps['VHL']+headmaps['POI']);restraint=np.array([5. if i in headids else .5 for i in range(nh)]);attempt_rows=[];best=None
 def objective(flat):
  x=flat.reshape(nall,3);heavy=x[:nh];energy=float(ff.CalcEnergy(flat.tolist()));grad=np.array(ff.CalcGrad(flat.tolist())).reshape(nall,3);delta=heavy-xyz;energy+=float(np.sum(restraint[:,None]*delta**2));grad[:nh]+=2*restraint[:,None]*delta
  hits=ptree.query_ball_point(heavy,3.9);ii=np.repeat(np.arange(nh),[len(h) for h in hits]);jj=np.array([j for h in hits for j in h],dtype=int)
  if len(ii):
   diff=heavy[ii]-protein[jj];dist=np.maximum(np.linalg.norm(diff,axis=1),1e-8);allowed=np.where(np.isin(np.array(els)[ii],['N','O'])&np.isin(np.array(pels)[jj],['N','O']),.6,.4);over=radii[ii]+prad[jj]-dist-allowed;active=over>0;ii,jj,diff,dist,over=ii[active],jj[active],diff[active],dist[active],over[active];energy+=float(100.*np.sum(over**2));np.add.at(grad,ii,-200.*over[:,None]*diff/dist[:,None])
  diff=heavy[i_internal]-heavy[j_internal];dist=np.maximum(np.linalg.norm(diff,axis=1),1e-8);over=rs_internal-dist-.6;active=over>0
  if active.any():
   ov=over[active];dv=diff[active];dd=dist[active];energy+=float(100.*np.sum(ov**2));gg=-200.*ov[:,None]*dv/dd[:,None];np.add.at(grad,i_internal[active],gg);np.add.at(grad,j_internal[active],-gg)
  return energy,grad.ravel()
 initialE,initialG=objective(init.ravel());finite_difference=[]
 for index in [0,17,80,nh*3+2]:
  flat=init.ravel().copy();step=1e-5;flat[index]+=step;plus=objective(flat)[0];flat[index]-=2*step;minus=objective(flat)[0];numeric=(plus-minus)/(2*step);error=abs(numeric-initialG[index]);c.require(error<1e-2*(1+abs(numeric)),'Objective gradient verification failed');finite_difference.append({'coordinate':index,'analytic':float(initialG[index]),'numeric':float(numeric),'absolute_error':float(error)})
 rng=np.random.default_rng(42);bounds=[]
 for i,point in enumerate(init):bounds.extend([(float(v-3),float(v+3)) if i<nh else (None,None) for v in point])
 for attempt in range(attempts):
  startx=init.copy()
  if attempt:startx[:nh]+=rng.normal(0,.12,(nh,3))
  result=minimize(objective,startx.ravel(),jac=True,method='L-BFGS-B',bounds=bounds,options={'maxiter':1200,'maxls':40,'ftol':1e-10,'gtol':1e-3});xx=result.x.reshape(nall,3);trial=Chem.Mol(hmol)
  for i,p in enumerate(xx):trial.GetConformer().SetAtomPosition(i,p)
  heavytrial=Chem.RemoveHs(trial);stereo=Chem.Mol(heavytrial);Chem.AssignAtomChiralTagsFromStructure(stereo,replaceExistingTags=True);key=Chem.MolToInchiKey(stereo);x=xx[:nh];external=c.overlap(x,els,protein,pels);internal=rs_internal-np.linalg.norm(x[i_internal]-x[j_internal],axis=1);bondlength=np.array([np.linalg.norm(x[names.index(a)]-x[names.index(b)]) for a,b,_ in bonds]);gradient=float(np.max(np.abs(objective(result.x)[1])));maxmove=float(np.linalg.norm(x-xyz,axis=1).max());record={'attempt':attempt,'optimizer_success':bool(result.success),'optimizer_message':str(result.message),'iterations':int(result.nit),'final_objective':float(result.fun),'gradient_max_abs':gradient,'coordinate_stereo_inchikey':key,'stereo_preserved':key==original_key,'max_displacement_A':maxmove,'external_max_overlap_A':external['max_overlap_A'],'internal_max_overlap_A':float(internal.max()),'bond_min_A':float(bondlength.min()),'bond_max_A':float(bondlength.max())};record['geometry_ok']=bool(key==original_key and external['max_overlap_A']<=.8 and internal.max()<=.8 and bondlength.min()>=1. and bondlength.max()<=2.1);record['convergence_ok']=bool(result.success or gradient<.01);attempt_rows.append(record)
  rank=(not(record['geometry_ok'] and record['convergence_ok']),max(external['max_overlap_A'],float(internal.max())),float(result.fun))
  if best is None or rank<best[0]:best=(rank,record,heavytrial,x.copy())
  if record['geometry_ok'] and record['convergence_ok']:break
 _,record,finalmol,newxyz=best;out=ROOT/'legacy_repairs'/row['gene'];c.require(not (out/'receipt.json').exists(),'Repair output already exists; archive before rerunning');out.mkdir(parents=True,exist_ok=True);writer=Chem.SDWriter(str(out/'repaired_historical_ligand.sdf'));writer.write(finalmol);writer.close();(out/'atom_name_map.json').write_text(json.dumps({str(i):n for i,n in enumerate(names)},indent=2)+'\n')
 lines=[]
 for line in Path(row['path']).read_text().splitlines():
  if line.startswith(('ATOM  ','HETATM')) and line[21]=='X':
   i=names.index(line[12:16].strip());x,y,z=newxyz[i];line=line[:30]+f'{x:8.3f}{y:8.3f}{z:8.3f}'+line[54:]
  if line.startswith(('ATOM  ','HETATM','TER','END')):lines.append(line)
 (out/'repaired_observed_complex.pdb').write_text('\n'.join(lines)+'\n');resultrow={**row,**record,'archived_coordinate_inchikey':original_key,'historical_graph_inchikey':Chem.MolToInchiKey(params),'supplier_equivalent':False,'pocket_engagement_validated':False,'new_MD':False,'seconds':time.monotonic()-start,'repair_method':'MMFF94s plus declared positional restraints and soft heavy-atom overlap penalties; proteins fixed; maximum heavy-coordinate shift 3A per Cartesian axis','initial_objective':initialE,'gradient_verification':finite_difference,'attempts':attempt_rows,'output_dir':str(out),'source_sha256':{'original_PDB':c.sha(row['path']),'full_params':c.sha(c.PARAM),'script':c.sha(__file__)}}
 for cutoff,tag in [(4.,'4A'),(4.5,'4p5A'),(5.,'5A')]:
  dv=cKDTree(vfull).query(newxyz[headmaps['VHL']])[0];dk=cKDTree(poi['xyz']).query(newxyz[headmaps['POI']])[0];resultrow['balanced_correct_head_contacts_'+tag]=int(min((dv<cutoff).sum(),(dk<cutoff).sum()))
 resultrow['structural_diagnostic_only']=True;resultrow['files_sha256']={p.name:c.sha(p) for p in out.iterdir() if p.is_file() and p.name!='receipt.json'};(out/'receipt.json').write_text(json.dumps(resultrow,indent=2)+'\n');return resultrow

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--limit',type=int,default=10);ap.add_argument('--attempts',type=int,default=2);ap.add_argument('--workers',type=int,default=4);args=ap.parse_args();c.require(1<=args.limit<=10 and 1<=args.attempts<=3,'Bounded repair allows original top10, up to3 starts each');table=pd.read_csv(BASE/'results/top10_static_distinct_prospective_kinases.tsv',sep='\t').head(args.limit);pools=json.loads((BASE/'config/pools.json').read_text());tasks=[(r,next(p for p in pools if p['uid']==r['uid'] and p['pool']==r['pool']),args.attempts) for r in table.to_dict('records')];results=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as ex:
  for r in ex.map(process,tasks):results.append(r);print(r['gene'],r['geometry_ok'],r['convergence_ok'],r['external_max_overlap_A'],r['internal_max_overlap_A'],flush=True)
 (ROOT/'legacy_repair_receipts.json').write_text(json.dumps(results,indent=2)+'\n');pd.DataFrame([{k:v for k,v in r.items() if not isinstance(v,(list,dict))} for r in results]).to_csv(ROOT/'legacy_repair_summary.tsv',sep='\t',index=False)
if __name__=='__main__':main()
