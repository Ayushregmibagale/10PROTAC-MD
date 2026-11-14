"""Fixed 0.5 ns sampling, heavy-atom NZ SASA and ligand-connected E2 geometry."""
from pathlib import Path
import os
os.environ['OMP_NUM_THREADS']='1'
import sys,json,hashlib,time,argparse,concurrent.futures
import numpy as np
import MDAnalysis as mda
import mdtraj as md
from MDAnalysis.lib.distances import distance_array
from Bio.SeqUtils import seq1
ROOT=Path(__file__).resolve().parents[1]
CAN=Path('/home/ayushregmibagale/10PROTAC-MD');OLD=Path('/home/ayushregmibagale/protac-md')
sys.path.insert(0,str(CAN/'prosettac_robust/workflow/scripts/lib'))
from protac_qc import pdb_atoms,atomic_json,sha256
from e2_geometry import WholeProteinTree,fit_reference_points
from coherent_images import assemble,ligand_connected_images,image_to_reference
from read_only_trajectory import ReadOnlyXTCReader
AA='ACDEFGHIKLMNPQRSTVWY';ALIASES={'HID':'H','HIE':'H','HIP':'H','HSD':'H','HSE':'H','HSP':'H','CYX':'C','MSE':'M'}
def require(ok,msg):
 if not ok:raise ValueError(msg)
def topology(group):
 top=md.Topology();chain=top.add_chain();prev=None
 for a in group:
  if a.resindex!=prev:res=top.add_residue(a.resname,chain,resSeq=int(a.resid));prev=a.resindex
  masses=np.array([12.011,14.007,15.999,32.06,30.974,18.998,35.45,79.904]);symbols=['C','N','O','S','P','F','Cl','Br'];k=int(np.argmin(abs(masses-a.mass)));require(abs(masses[k]-a.mass)<.3,'Unknown element mass')
  top.add_atom(a.name,md.element.get_by_symbol(symbols[k]),res)
 return top

def work(case):
 poi,pose=case;out=ROOT/'analysis/lysines/poses'/poi/pose;receipt=out/'receipt.json'
 if receipt.exists():
  r=json.loads(receipt.read_text());require(sha256(out/'packet.npz')==r['packet_sha256'],'Cached packet drift');require(r['script_sha256']==sha256(__file__),'Extraction source drift');return r
 start=time.monotonic();base=OLD/'random_screen'/poi/pose;tpr=base/'prod/prod.tpr';xtc=OLD/'analysis/random_screen'/poi/pose/'prod_nojump.xtc';ref=base/'structure/complex_protein.pdb'
 contact=CAN/'prosettac_robust/existing_result_audit/contact_core_v1/poses'/poi/pose/'contacts_A.tsv';prior=CAN/'SJF8240_reanalysis/analysis/e2/ligand_connected_e2_cohort_v1/poses'/poi/pose/'frames.tsv';ensemble_path=OLD/'reference_structures/e2_c85_ensemble.json'
 sources={str(p):sha256(p) for p in [tpr,xtc,ref,contact,prior,ensemble_path]};u=mda.Universe(str(tpr),str(xtc),format=ReadOnlyXTCReader)
 protein=u.select_atoms('protein');whole=WholeProteinTree(protein);ligand=u.select_atoms('resname UNL');lw=WholeProteinTree(ligand);lig=ligand[ligand.masses>1.5]
 ca={c:u.select_atoms('moltype Protein_chain_'+c+' and name CA') for c in 'ABCD'};heavy={c:u.select_atoms('moltype Protein_chain_'+c) for c in 'ABCD'};heavy={c:g[g.masses>1.5] for c,g in heavy.items()}
 allheavy=heavy['A']+heavy['B']+heavy['C']+heavy['D']+lig;nz=u.select_atoms('moltype Protein_chain_D and resname LYS and name NZ')
 atoms=pdb_atoms(ref);reference={c:np.asarray([a['xyz'] for a in atoms if a['record']=='ATOM' and a['chain']==c and a['name']=='CA']) for c in 'ABCD'}
 for c in 'ABCD':require(len(ca[c])==len(reference[c]),'Reference CA count')
 ensemble=np.asarray(next(iter(json.loads(ensemble_path.read_text()).values())),dtype=float)
 indices=np.array([int(np.flatnonzero(allheavy.indices==i)[0]) for i in nz.indices]);poiindices=np.array([int(np.flatnonzero(heavy['D'].indices==i)[0]) for i in nz.indices]);lysmap=np.array([int(np.flatnonzero(ca['D'].resindices==i)[0]) for i in nz.resindices]);require(len(set(lysmap))==len(nz),'Lysine residue mapping')
 contacts=np.genfromtxt(contact,delimiter='\t',names=True);old=np.genfromtxt(prior,delimiter='\t',names=True);sample=np.arange(0,201,5);coords=[];isolated=[];dist=[];reach=[];bridges=[];times=[];errors=[]
 window_edges={role:{k:[] for k in ['left','right']} for role in ['elob','vhl']}
 for index in sample:
  ts=u.trajectory[int(old['frame'][index])];require(abs(ts.time/1000-index*.1)<5.1e-5,'Sampling time mismatch');whole(ts);lw(ts);box=ts.dimensions
  im=assemble(ca['B'].positions,ca['C'].positions,ca['D'].positions,reference['B'],reference['C'],reference['D'],box);lc=ligand_connected_images(heavy['C'].positions,lig.positions,heavy['D'].positions,im['vhl_shift_A'],box)
  predA,_,_=fit_reference_points(reference['B'],ca['B'].positions,reference['A']);_,ashift,_=image_to_reference(ca['A'].positions,predA,box)
  xyz=np.concatenate([heavy['A'].positions+ashift,heavy['B'].positions,heavy['C'].positions+im['vhl_shift_A'],heavy['D'].positions+lc['poi_shift_A'],lig.positions+lc['ligand_shift_A']]);coords.append(xyz/10);isolated.append((heavy['D'].positions+lc['poi_shift_A'])/10)
  points,_,_=fit_reference_points(reference['B'],ca['B'].positions,ensemble);d=distance_array(nz.positions+lc['poi_shift_A'],points);dist.append(d.min(1));reach.append((d<80).mean(1));bridges.append(bool(contacts['true_VHL_heavy_A'][index]<5 and contacts['POI_heavy_A'][index]<5));times.append(ts.time/1000)
  errors.append(abs(d.min()-old['ligand_image_min_E2_A'][index]));require(max(abs(lc['VHL_ligand_min_A']-contacts['true_VHL_heavy_A'][index]),abs(lc['ligand_POI_min_A']-contacts['POI_heavy_A'][index]))<.002,'Saved contact mismatch')
  right=distance_array(lig.positions,ca['D'].positions,box=box).min(0)
  for role,c in [('elob','A'),('vhl','C')]:window_edges[role]['left'].append(distance_array(lig.positions,ca[c].positions,box=box).min(0));window_edges[role]['right'].append(right)
 require(max(errors)<.002,'E2 descriptor mismatch');u.trajectory.close()
 c_sasa=md.shrake_rupley(md.Trajectory(np.asarray(coords,dtype=np.float32),topology(allheavy)),n_sphere_points=960,atom_indices=indices)[:,indices]*100
 p_sasa=md.shrake_rupley(md.Trajectory(np.asarray(isolated,dtype=np.float32),topology(heavy['D'])),n_sphere_points=960,atom_indices=poiindices)[:,poiindices]*100
 require(np.isfinite(c_sasa).all() and (c_sasa>=0).all() and np.min(p_sasa-c_sasa)>-.001,'Invalid or inconsistent SASA')
 arrays={'times':np.array(times),'lys_resids':nz.resids,'lys_node_indices':lysmap,'complex_sasa_A2':c_sasa,'isolated_sasa_A2':p_sasa,'burial_A2':np.maximum(p_sasa-c_sasa,0),'min_e2_A':np.array(dist),'reach80_fraction':np.array(reach),'bridge':np.array(bridges)}
 for c in 'ACD':arrays['resids_'+c]=ca[c].resids;arrays['resnames_'+c]=np.array(ca[c].resnames,dtype='U4')
 for role,vals in window_edges.items():
  for key,val in vals.items():arrays[role+'_'+key]=np.array(val)
 out.mkdir(parents=True,exist_ok=True);np.savez_compressed(out/'packet.npz',**arrays)
 require(all(sha256(Path(p))==h for p,h in sources.items()),'Source input changed')
 r={'poi':poi,'pose':pose,'frames':41,'lysines':len(nz),'source_sha256':sources,'script_sha256':sha256(__file__),'packet_sha256':sha256(out/'packet.npz'),'max_e2_replay_error_A':max(errors),'seconds':time.monotonic()-start,'mdtraj':md.__version__,'sample_step_ns':.5,'SASA_probe_A':1.4,'sphere_points':960,'heavy_atom_occluders_only':True};atomic_json(receipt,r);print(json.dumps(r),flush=True);return r

def main():
 p=argparse.ArgumentParser();p.add_argument('--limit',type=int);p.add_argument('--workers',type=int,default=2);a=p.parse_args();records=json.loads((CAN/'SJF8240_reanalysis/analysis/dynamic/completion.json').read_text())['records'];labels=json.loads((ROOT/'config/protocol.json').read_text())['labels'];cases=[(r['poi'],r['pose_id']) for r in records if r['poi'] in labels];require(len(cases)==141,'Core roster mismatch')
 if a.limit:cases=cases[:a.limit]
 with concurrent.futures.ProcessPoolExecutor(max_workers=a.workers) as pool:rows=list(pool.map(work,cases))
 atomic_json(ROOT/'analysis/lysines'/('benchmark.json' if a.limit else 'completion.json'),{'cases':len(rows),'records':rows,'sampling':'41 frames at 0.5 ns','new_MD':False})
if __name__=='__main__':main()
