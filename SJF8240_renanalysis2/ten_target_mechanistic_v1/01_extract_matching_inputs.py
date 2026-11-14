"""Extract the existing feature definitions from retained trajectories, without training."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,hashlib,time,argparse,concurrent.futures,importlib.util
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'lib'))
import numpy as np,MDAnalysis as mda,mdtraj as md
from MDAnalysis.lib.distances import distance_array
from scipy.spatial import cKDTree
from Bio.SeqUtils import seq1
from e2_geometry import WholeProteinTree,fit_reference_points
from coherent_images import assemble,ligand_connected_images,image_to_reference
from read_only_trajectory import ReadOnlyXTCReader
from protac_qc import pdb_atoms
BASE=Path('/home/ayushregmibagale/10PROTAC-MD');CAN=BASE/'SJF8240_reanalysis';NATIVE=BASE/'SJF8240_renanalysis2/native_catalytic_validation_v2';MECH=CAN/'mechanistic_v3'
def load(name,p):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
native=load('native_snapshot',ROOT/'source_snapshots/native_extractor.py');legacy=load('legacy_snapshot',ROOT/'source_snapshots/lysine_extractor.py')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):
 p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(x,indent=2,allow_nan=False));tmp.replace(p)
def require(x,msg):
 if not x:raise RuntimeError(msg)
def freeze():
 for p,h in json.loads((ROOT/'execution_manifest.json').read_text()).items():require(sha(p)==h,'Frozen source drift '+p)
 return sha(ROOT/'execution_manifest.json')
def work(task):
 target,source,pose=task;frozen=freeze();out=ROOT/'extracted'/target/pose;out.mkdir(parents=True,exist_ok=True);receipt=out/'receipt.json'
 if receipt.exists():
  r=json.loads(receipt.read_text());require(r['manifest_sha256']==frozen and all(sha(out/k)==v for k,v in r['outputs'].items()),'Cached drift');return r
 start=time.monotonic();oldreceipt=CAN/'analysis/dynamic/poses'/source/pose/'receipt.json';old=json.loads(oldreceipt.read_text());inputs={str(oldreceipt):sha(oldreceipt),**old['source_sha256']}
 for p,h in inputs.items():require(sha(p)==h,'Trajectory drift '+p)
 tpr=next(Path(p) for p in inputs if p.endswith('.tpr'));xtc=next(Path(p) for p in inputs if p.endswith('.xtc'));pdb=tpr.parents[1]/'structure/complex_protein.pdb';inputs[str(pdb)]=sha(pdb)
 evidence=Path(old['sequence_identity_evidence']);inputs[str(evidence)]=sha(evidence);roles=json.loads(evidence.read_text());role=next(x for x in roles if x['poi']==source and x['pose_id']==pose)
 if isinstance(role['proteins'],dict):
  require(role['source_sha256']['tpr']==sha(tpr),'TPR not bound');proteins=role['proteins']
 else:
  require(role['tpr_sha256']==sha(tpr),'TPR not bound');proteins={r['moltype']:r for r in role['proteins']}
 expected=json.loads((CAN/'config/role_reference.json').read_text())['targets'][target]
 for c,uid in [('A','Q15370'),('B','Q15369'),('C','P40337'),('D',expected)]:require(proteins['Protein_chain_'+c]['uniprot']==uid,'Wrong role '+c)
 u=mda.Universe(str(tpr),str(xtc),format=ReadOnlyXTCReader);whole=WholeProteinTree(u.select_atoms('protein'));ligall=u.select_atoms('resname UNL');lw=WholeProteinTree(ligall);lig=ligall[ligall.masses>1.5];require(len(lig)==79,'Different ligand atom count')
 ca={c:u.select_atoms('moltype Protein_chain_'+c+' and name CA') for c in 'ABCD'};heavy={c:u.select_atoms('moltype Protein_chain_'+c+' and mass 1.5:1000') for c in 'ABCD'};nz=u.select_atoms('moltype Protein_chain_D and resname LYS and name NZ');require(len(nz)>0,'No lysines');ce=u.atoms[[next(a.index for a in r.atoms if a.name=='CE') for r in nz.residues]];atoms=pdb_atoms(pdb);ref={c:np.array([a['xyz'] for a in atoms if a['record']=='ATOM' and a['chain']==c and a['name']=='CA']) for c in 'ABCD'}
 for c in 'ABCD':require(ref[c].shape==ca[c].positions.shape,'CA mismatch')
 ensemble_path=Path('/home/ayushregmibagale/protac-md/reference_structures/e2_c85_ensemble.json');inputs[str(ensemble_path)]=sha(ensemble_path);ensemble=np.asarray(next(iter(json.loads(ensemble_path.read_text()).values())),float)
 specs=[]
 for name in ['8R5H_observed_SY8','8RX0_G76_fit75','8RX0_E2_aligned_SY8']:
  path=NATIVE/'references'/(name+'.npz');inputs[str(path)]=sha(path);z=np.load(path);xyz=z['xyz'];radii=np.array([native.RAD[e] for e in z['elements']]);tree=cKDTree(xyz)
  for anchor,c in [('VHL','C'),('ELOC','B')]:
   ii,jj=native.matched(str(z[anchor+'_sequence']),''.join(seq1(r) for r in ca[c].resnames));specs.append({'name':name+'__'+anchor,'chain':c,'ii':ii,'jj':jj,'ref':z[anchor+'_CA'][ii],'xyz':xyz,'radii':radii,'tree':tree,'donor':z['donor_C']})
 allheavy=heavy['A']+heavy['B']+heavy['C']+heavy['D']+lig;nzids=np.array([np.flatnonzero(allheavy.indices==i)[0] for i in nz.indices]);poiids=np.array([np.flatnonzero(heavy['D'].indices==i)[0] for i in nz.indices]);lysmap=np.array([np.flatnonzero(ca['D'].resindices==i)[0] for i in nz.resindices]);require(len(np.unique(lysmap))==len(nz),'Ambiguous lysine mapping')
 rad=np.array([native.RAD[a.name.lstrip('0123456789')[0]] for a in allheavy]);targetrad=np.array([native.RAD[a.name.lstrip('0123456789')[0]] for a in heavy['D']]);nf=201;nl=len(nz);dist=np.zeros((6,nf,1,nl),np.float32);cos=np.zeros_like(dist);exposure=np.zeros((nf,nl),np.float32);clash=np.zeros((6,nf),np.float32);fits=np.zeros_like(clash);min_e2=[];reach50=[];reach80=[];bridges=[];contacts=[];times=[];coords=[];isolated=[];base_dist=[];base_reach=[];base_bridge=[];edges={role:{k:[] for k in ['left','right']} for role in ['elob','vhl']};contacterrors=[];timestamps=np.array([ts.time for ts in u.trajectory]);wanted=np.arange(nf)*100.;indices=[]
 for t in wanted:
  ids=np.flatnonzero(np.abs(timestamps-t)<.051);require(len(ids)==1,'Missing or duplicate frame time');indices.append(int(ids[0]))
 for index,frame in enumerate(indices):
  ts=u.trajectory[frame];whole(ts);lw(ts);box=ts.dimensions;im=assemble(ca['B'].positions,ca['C'].positions,ca['D'].positions,ref['B'],ref['C'],ref['D'],box);lc=ligand_connected_images(heavy['C'].positions,lig.positions,heavy['D'].positions,im['vhl_shift_A'],box);predA,_,_=fit_reference_points(ref['B'],ca['B'].positions,ref['A']);_,ashift,_=image_to_reference(ca['A'].positions,predA,box);shifts={'A':ashift,'B':np.zeros(3),'C':im['vhl_shift_A'],'D':lc['poi_shift_A']};targetxyz=heavy['D'].positions+shifts['D'];nzxyz=nz.positions+shifts['D'];cexyz=ce.positions+shifts['D'];occluders=np.vstack([heavy[c].positions+shifts[c] for c in 'ABCD']+[lig.positions+lc['ligand_shift_A']]);exposure[index]=native.exposed_fraction(nzxyz,nz.indices,occluders,allheavy.indices,rad)
  independent=[distance_array(heavy[c].positions,lig.positions,box=box).min() for c in ['C','D']];error=max(abs(independent[0]-lc['VHL_ligand_min_A']),abs(independent[1]-lc['ligand_POI_min_A']));require(error<.002,'PBC contact mismatch');contacterrors.append(error);bridge=bool(lc['VHL_ligand_min_A']<5 and lc['ligand_POI_min_A']<5);bridges.append(bridge);contacts.append([lc['VHL_ligand_min_A'],lc['ligand_POI_min_A']]);times.append(ts.time/1000)
  points,_,_=fit_reference_points(ref['B'],ca['B'].positions,ensemble);d=distance_array(nzxyz,points);min_e2.append(d.min(1));reach50.append((d<50).mean(1));reach80.append((d<80).mean(1))
  for si,spec in enumerate(specs):
   mobile=ca[spec['chain']].positions[spec['jj']]+shifts[spec['chain']];_,rmsd,rot=fit_reference_points(spec['ref'],mobile,spec['ref']);fits[si,index]=rmsd;localtarget=(targetxyz-mobile.mean(0))@rot+spec['ref'].mean(0);localnz=(nzxyz-mobile.mean(0))@rot+spec['ref'].mean(0);localce=(cexyz-mobile.mean(0))@rot+spec['ref'].mean(0);clash[si,index]=native.overlap_positive(localtarget,targetrad,spec['tree'],spec['xyz'],spec['radii'])[0];vec=spec['donor'][:,None,:]-localnz[None,:,:];dd=np.linalg.norm(vec,axis=2);dist[si,index]=dd;cos[si,index]=np.einsum('eli,li->el',vec,localnz-localce)/np.maximum(dd*np.linalg.norm(localnz-localce,axis=1)[None,:],1e-12)
  if index%5==0:
   coords.append(occluders/10);isolated.append(targetxyz/10);base_dist.append(d.min(1));base_reach.append((d<80).mean(1));base_bridge.append(bridge)
   right=distance_array(lig.positions,ca['D'].positions,box=box).min(0)
   for role,c in [('elob','A'),('vhl','C')]:edges[role]['left'].append(distance_array(lig.positions,ca[c].positions,box=box).min(0));edges[role]['right'].append(right)
 u.trajectory.close();require(np.allclose(times,np.arange(nf)/10,atol=5.1e-5),'Timestamp drift')
 c_sasa=md.shrake_rupley(md.Trajectory(np.asarray(coords,dtype=np.float32),legacy.topology(allheavy)),n_sphere_points=960,atom_indices=nzids)[:,nzids]*100;p_sasa=md.shrake_rupley(md.Trajectory(np.asarray(isolated,dtype=np.float32),legacy.topology(heavy['D'])),n_sphere_points=960,atom_indices=poiids)[:,poiids]*100;require(np.min(p_sasa-c_sasa)>-.001,'Negative SASA burial')
 arrays={'times':np.array(times)[::5],'lys_resids':nz.resids,'lys_node_indices':lysmap,'complex_sasa_A2':c_sasa,'isolated_sasa_A2':p_sasa,'burial_A2':np.maximum(p_sasa-c_sasa,0),'min_e2_A':np.array(base_dist),'reach80_fraction':np.array(base_reach),'bridge':np.array(base_bridge)}
 for c in 'ACD':arrays['resids_'+c]=ca[c].resids;arrays['resnames_'+c]=np.array(ca[c].resnames,dtype='U4')
 for role,vals in edges.items():
  for k,v in vals.items():arrays[role+'_'+k]=np.array(v)
 np.savez_compressed(out/'base_packet.npz',**arrays);np.savez_compressed(out/'native_packet.npz',times_ns=np.array(times),frame_indices=indices,lys_resids=nz.resids,bridge=bridges,reach50=(np.array(min_e2)<50),reach50_fraction=np.array(reach50),reach80_fraction=np.array(reach80),min_e2_A=np.array(min_e2),nz_exposed_fraction=exposure,donor_distance_A=dist,sidechain_outward_cosine=cos,target_scaffold_max_overlap_A=clash,anchor_RMSD_A=fits,spec_names=np.array([s['name'] for s in specs]),contacts_A=np.array(contacts))
 for p,h in inputs.items():require(sha(p)==h,'Input changed during extraction')
 r={'target':target,'source_target':source,'pose':pose,'frames':nf,'lysines':nl,'seconds':time.monotonic()-start,'max_contact_replay_error_A':max(contacterrors),'source_sha256':inputs,'manifest_sha256':frozen,'outputs':{n:sha(out/n) for n in ['base_packet.npz','native_packet.npz']}};save(receipt,r);return r

def cases():
 dynamic=json.loads((CAN/'analysis/dynamic/completion.json').read_text())['records'];unknown=[(r['poi']+'_full',r['poi'],r['pose_id']) for r in dynamic if r['poi'] in ['SIK3','RPS6KB2']];labels=json.loads((NATIVE/'config/labels.json').read_text())['labels'];controls=[(t,t,sorted(r['pose_id'] for r in dynamic if r['poi']==t)[0]) for t in labels];return sorted(unknown)+sorted(controls)
def main():
 p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);p.add_argument('--limit',type=int);a=p.parse_args();todo=cases();todo=todo[:a.limit] if a.limit else todo;rows=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=a.workers) as ex:
  for r in ex.map(work,todo):rows.append(r);save(ROOT/'STATUS.json',{'phase':'extracting_existing_trajectories','completed':len(rows),'expected':len(todo)});print(len(rows),r['target'],r['pose'],round(r['seconds'],1),flush=True)
 save(ROOT/('probe_completion.json' if a.limit else 'extraction_completion.json'),rows)
if __name__=='__main__':main()
