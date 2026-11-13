"""Four fixed geometry/mask combinations on identical full201-frame MD grids."""
from pathlib import Path
import os
os.environ['OMP_NUM_THREADS']='1'
import sys,json,time,argparse,concurrent.futures,hashlib
import numpy as np,pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from Bio.SeqUtils import seq1
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'))
from protac_qc import atomic_json,sha256,pdb_atoms
from read_only_trajectory import ReadOnlyXTCReader
from e2_geometry import WholeProteinTree,fit_reference_points
from coherent_images import assemble,ligand_connected_images
CAN=Path('/home/ayushregmibagale/10PROTAC-MD/SJF8240_reanalysis');OLD=Path('/home/ayushregmibagale/protac-md')
def require(ok,msg):
 if not ok:raise ValueError(msg)
def verify_source():
 f=json.loads((ROOT/'config/execution_manifest.json').read_text())
 for rel,h in f.items():require(sha256(ROOT/rel)==h,'Frozen source drift '+rel)
 return sha256(ROOT/'config/execution_manifest.json')
def legacy_reference(poi):
 paths=[OLD/'simulations'/poi/'structure/rank1/complex_protein.pdb'] if poi!='EPHB2' else []
 paths+=sorted((OLD/'random_screen'/poi).glob('pose*/structure/complex_protein.pdb'))
 for path in paths:
  if not path.exists():continue
  atoms=pdb_atoms(path);rows=[a for a in atoms if a['record']=='ATOM' and a['chain']=='B' and a['name']=='CA']
  if len(rows)==85:return np.array([a['xyz'] for a in rows],float),path
 raise ValueError('No original ElonginC reference '+poi)
def paired_delta(distances,bridge,cutoff=50):
 d=np.asarray(distances);b=np.asarray(bridge,bool);require(d.ndim==2 and b.shape==(len(d),) and np.isfinite(d).all(),'Invalid arrays')
 if not b.any() or not (~b).any():return None
 return (d[b]<cutoff).mean(0)-(d[~b]<cutoff).mean(0)
def work(case,independent=False,output=None):
 poi,pose=case;out=(ROOT/'analysis/poses' if output is None else Path(output))/poi/pose;receipt=out/'receipt.json';frozen=verify_source()
 if receipt.exists():
  r=json.loads(receipt.read_text());require(r['source_manifest_sha256']==frozen and sha256(out/'packet.npz')==r['packet_sha256'],'Cached packet/source drift');return r
 start=time.monotonic();base=OLD/'random_screen'/poi/pose;tpr=base/'prod/prod.tpr';xtc=OLD/'analysis/random_screen'/poi/pose/'prod_nojump.xtc';pdb=base/'structure/complex_protein.pdb';metrics=OLD/'analysis/random_screen'/poi/pose/'metrics.csv';contacts=CAN/'../prosettac_robust/existing_result_audit/contact_core_v1/poses'/poi/pose/'contacts_A.tsv';prior=CAN/'analysis/e2/ligand_connected_e2_cohort_v1/poses'/poi/pose/'frames.tsv';ensemblefile=OLD/'reference_structures/e2_c85_ensemble.json';legacyref,legacyrefpath=legacy_reference(poi)
 files=[tpr,xtc,pdb,metrics,contacts,prior,ensemblefile,legacyrefpath];sources={str(p.resolve()):sha256(p) for p in files};u=mda.Universe(str(tpr),str(xtc),format=ReadOnlyXTCReader);whole=WholeProteinTree(u.select_atoms('protein'));ligall=u.select_atoms('resname UNL');ligwhole=WholeProteinTree(ligall);lig=ligall[ligall.masses>1.5];ca={c:u.select_atoms('moltype Protein_chain_'+c+' and name CA') for c in 'BCD'};nz=u.select_atoms('moltype Protein_chain_D and resname LYS and name NZ');vh=u.select_atoms('moltype Protein_chain_C');vh=vh[vh.masses>1.5];po=u.select_atoms('moltype Protein_chain_D');po=po[po.masses>1.5]
 require(len(lig)==79 and len(nz)>0 and len(ca['B'])==85,'Selection differs from reviewed compound')
 # Bind exact topology and role assignment to the independent sequence audit.
 roles=json.loads((CAN/'../prosettac_robust/existing_result_audit/protein_roles_core_complete/per_pose_roles.json').read_text());row=next(x for x in roles if x['poi']==poi and x['pose_id']==pose);require(row['source_sha256']['tpr']==sources[str(tpr)],'Topology identity drift');expected=json.loads((ROOT/'config/role_reference.json').read_text())['targets'][poi]
 for chain,uid in [('B','Q15369'),('C','P40337'),('D',expected)]:require(row['proteins']['Protein_chain_'+chain]['uniprot']==uid,'Wrong molecular role')
 atoms=pdb_atoms(pdb);reference={c:np.array([a['xyz'] for a in atoms if a['record']=='ATOM' and a['chain']==c and a['name']=='CA'],float) for c in 'BCD'}
 for c in 'BCD':require(reference[c].shape==ca[c].positions.shape,'Reference atom correspondence')
 arrays=[np.asarray(v,float) for v in json.loads(ensemblefile.read_text()).values()];require(all(np.array_equal(a,arrays[0]) for a in arrays),'E2 duplication changed');points_ref=arrays[0]
 old=pd.read_csv(prior,sep='\t');contact=pd.read_csv(contacts,sep='\t');m=pd.read_csv(metrics,comment='#');grid=np.arange(201)/10
 require(len(old)==len(contact)==201 and np.allclose(old.time_ns,grid,atol=5.1e-5,rtol=0) and np.allclose(contact.time_ns,grid,atol=5.1e-5,rtol=0),'Reviewed frame grid mismatch')
 m=m[(m.time_ns>=-5.1e-5)&(m.time_ns<=20+5.1e-5)].copy()
 if m.time_ns.duplicated().any():
  for _,g in m.groupby('time_ns'):require(np.allclose(g[['d_protac_vhl_A','d_protac_kinase_A']],g[['d_protac_vhl_A','d_protac_kinase_A']].iloc[0],atol=0,rtol=0),'Conflicting duplicate metrics')
  m=m.drop_duplicates('time_ns')
 require(len(m)==201 and np.allclose(m.time_ns,grid,atol=5.1e-5,rtol=0),'Metrics grid mismatch')
 legacy_mask=((m.d_protac_vhl_A<5)&(m.d_protac_kinase_A<5)).to_numpy();true_mask=((contact.true_VHL_heavy_A<5)&(contact.POI_heavy_A<5)).to_numpy();raw=[];corrected=[];times=[];errors=[];contact_errors=[]
 if independent:
  from MDAnalysis.transformations import unwrap
  reference_unwrap=unwrap(u.select_atoms('protein'))
 for index in range(201):
  ts=u.trajectory[int(old.frame.iloc[index])];require(abs(ts.time/1000-grid[index])<5.1e-5,'Trajectory time mismatch')
  if independent:reference_unwrap(ts)
  else:whole(ts)
  ligwhole(ts);anchor=ca['B'].positions.copy();v=ca['C'].positions.copy();p=ca['D'].positions.copy();nzpos=nz.positions.copy()
  oldpoints,_,_=fit_reference_points(legacyref,anchor,points_ref);raw.append(distance_array(nzpos,oldpoints).min(1))
  im=assemble(anchor,v,p,reference['B'],reference['C'],reference['D'],ts.dimensions);lc=ligand_connected_images(vh.positions,lig.positions,po.positions,im['vhl_shift_A'],ts.dimensions);points,_,_=fit_reference_points(reference['B'],anchor,points_ref);d=distance_array(nzpos+lc['poi_shift_A'],points).min(1);corrected.append(d);times.append(ts.time/1000);errors.append(abs(float(d.min())-old.ligand_image_min_E2_A.iloc[index]));contact_errors.append(max(abs(lc['VHL_ligand_min_A']-contact.true_VHL_heavy_A.iloc[index]),abs(lc['ligand_POI_min_A']-contact.POI_heavy_A.iloc[index])))
 u.trajectory.close();raw=np.asarray(raw);corrected=np.asarray(corrected);require(max(errors)<.002 and max(contact_errors)<.002,'Previously audited geometry/contact replay failed');require(np.isfinite(raw).all() and np.isfinite(corrected).all(),'Nonfinite distances')
 packet={'times_ns':np.asarray(times),'lys_resids':nz.resids.copy(),'legacy_min_e2_A':raw,'corrected_min_e2_A':corrected,'legacy_mask':legacy_mask,'true_vhl_mask':true_mask}
 out.mkdir(parents=True,exist_ok=True);tmp=out/'packet.pending.npz';np.savez_compressed(tmp,**packet);os.replace(tmp,out/'packet.npz');require(all(sha256(Path(p))==h for p,h in sources.items()),'Input changed during computation');require(verify_source()==frozen,'Source changed during computation')
 r={'poi':poi,'pose':pose,'frames':201,'lysines':len(nz),'source_manifest_sha256':frozen,'input_sha256':sources,'packet_sha256':sha256(out/'packet.npz'),'max_prior_distance_error_A':max(errors),'max_prior_contact_error_A':max(contact_errors),'legacy_bridge_frames':int(legacy_mask.sum()),'true_VHL_bridge_frames':int(true_mask.sum()),'legacy_paired':bool(legacy_mask.any() and (~legacy_mask).any()),'true_VHL_paired':bool(true_mask.any() and (~true_mask).any()),'independent_MDA_unwrap':independent,'seconds':time.monotonic()-start};atomic_json(receipt,r);return r

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=2);parser.add_argument('--limit',type=int);a=parser.parse_args();labels=json.loads((ROOT/'config/protocol.json').read_text())['labels'];records=json.loads((CAN/'analysis/dynamic/completion.json').read_text())['records'];cases=[(r['poi'],r['pose_id']) for r in records if r['poi'] in labels];require(len(cases)==141,'Confirmed cohort roster mismatch')
 if a.limit:cases=cases[:a.limit]
 rows=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=a.workers) as pool:
  for row in pool.map(work,cases):rows.append(row);atomic_json(ROOT/'STATUS.json',{'phase':'EXTRACTING','complete':len(rows),'expected':len(cases),'active_job':True});print(json.dumps({'complete':len(rows),'expected':len(cases),'poi':row['poi'],'pose':row['pose']}),flush=True)
 atomic_json(ROOT/'analysis'/('benchmark.json' if a.limit else 'completion.json'),{'status':'FULL_FRAME_EXTRACTION_COMPLETE','cases':len(rows),'frames':sum(r['frames'] for r in rows),'records':rows})
if __name__=='__main__':main()
