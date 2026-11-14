"""Project existing MD onto explicit CRL2 references; retain all denominators."""
from pathlib import Path
import os,sys,json,time,argparse,concurrent.futures,hashlib
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1';sys.dont_write_bytecode=True
import numpy as np,pandas as pd
import MDAnalysis as mda
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import seq1
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'))
from read_only_trajectory import ReadOnlyXTCReader
from e2_geometry import WholeProteinTree,fit_reference_points
from coherent_images import assemble,ligand_connected_images,image_to_reference
from protac_qc import sha256,pdb_atoms,atomic_json
CAN=Path('/home/ayushregmibagale/10PROTAC-MD/SJF8240_reanalysis');OLD=Path('/home/ayushregmibagale/protac-md');PRIOR=CAN/'full_frame_e2_delta50_v1'
RAD={'C':1.70,'N':1.55,'O':1.52,'S':1.80,'F':1.47,'P':1.80}
def require(ok,msg):
 if not ok:raise ValueError(msg)
def matched(reference,query):
 align=PairwiseAligner(mode='local',match_score=2,mismatch_score=-3,open_gap_score=-5,extend_gap_score=-.5).align(reference,query)[0];pairs=[(i,j) for (a,b),(c,d) in zip(*align.aligned) for i,j in zip(range(a,b),range(c,d))];require(len(pairs)>=60 and sum(reference[i]==query[j] for i,j in pairs)/len(pairs)>=.95,'Anchor sequence correspondence failed');return np.array(pairs,int).T

def sphere_points(n=180):
 i=np.arange(n);z=1-2*(i+.5)/n;phi=i*np.pi*(3-np.sqrt(5));rr=np.sqrt(1-z*z);return np.c_[rr*np.cos(phi),rr*np.sin(phi),z]
SPHERE=sphere_points()
def exposed_fraction(nz,nzids,xyz,ids,radii):
 tree=cKDTree(xyz);output=[]
 for point,own in zip(nz,nzids):
  near=np.array(tree.query_ball_point(point,6.5),int);near=near[ids[near]!=own];surface=point+SPHERE*(1.55+1.4)
  if len(near):occluded=np.any(np.sum((surface[:,None,:]-xyz[near][None,:,:])**2,axis=2)<(radii[near]+1.4)**2,axis=1);output.append(float((~occluded).mean()))
  else:output.append(1.)
 return np.array(output)

def overlap_positive(x,radii,tree,refxyz,refradii):
 hits=tree.query_ball_point(x,3.8);ii=np.repeat(np.arange(len(x)),[len(h) for h in hits]);jj=np.array([j for h in hits for j in h],int)
 if not len(ii):return 0.,-1
 vv=radii[ii]+refradii[jj]-np.linalg.norm(x[ii]-refxyz[jj],axis=1);k=int(np.argmax(vv));return max(0.,float(vv[k])),int(jj[k])

def verify_source():
 f=json.loads((ROOT/'config/execution_manifest.json').read_text())
 for rel,h in f.items():require(sha256(ROOT/rel)==h,'Frozen source drift '+rel)
 return sha256(ROOT/'config/execution_manifest.json')

def work(case):
 poi,pose=case;out=ROOT/'results/poses'/poi/pose;out.mkdir(parents=True,exist_ok=True);frozen=verify_source();receipt=out/'receipt.json'
 if receipt.exists():
  r=json.loads(receipt.read_text());require(r['source_manifest_sha256']==frozen and sha256(out/'packet.npz')==r['packet_sha256'],'Cached source/result drift');return r
 start=time.monotonic();prior=PRIOR/'analysis/poses'/poi/pose;pr=json.loads((prior/'receipt.json').read_text());tpr=OLD/'random_screen'/poi/pose/'prod/prod.tpr';xtc=OLD/'analysis/random_screen'/poi/pose/'prod_nojump.xtc';pdb=OLD/'random_screen'/poi/pose/'structure/complex_protein.pdb';sources={str(p):sha256(p) for p in [tpr,xtc,pdb,prior/'packet.npz',prior/'receipt.json']}
 for p in [tpr,xtc,pdb]:require(sources[str(p)]==pr['input_sha256'][str(p)],'Audited MD source drift')
 roles=json.loads((CAN/'../prosettac_robust/existing_result_audit/protein_roles_core_complete/per_pose_roles.json').read_text());role=next(r for r in roles if r['poi']==poi and r['pose_id']==pose);require(role['source_sha256']['tpr']==sources[str(tpr)],'Role audit not bound to topology')
 expected=json.loads((CAN/'config/role_reference.json').read_text())['targets'][poi]
 for chain,uid in [('B','Q15369'),('C','P40337'),('D',expected)]:require(role['proteins']['Protein_chain_'+chain]['uniprot']==uid,'Wrong protein identity')
 previous=np.load(prior/'packet.npz');u=mda.Universe(str(tpr),str(xtc),format=ReadOnlyXTCReader);protein=u.select_atoms('protein');whole=WholeProteinTree(protein);ligall=u.select_atoms('resname UNL');ligwhole=WholeProteinTree(ligall);lig=ligall[ligall.masses>1.5];groups={c:u.select_atoms('moltype Protein_chain_'+c) for c in 'ABCD'};heavy={c:g[g.masses>1.5] for c,g in groups.items()};ca={c:u.select_atoms('moltype Protein_chain_'+c+' and name CA') for c in 'ABCD'};nz=u.select_atoms('moltype Protein_chain_D and resname LYS and name NZ');ce=u.atoms[[int(next(a.index for a in r.atoms if a.name=='CE')) for r in nz.residues]];require(len(lig)==79 and len(nz)>0,'Ligand/lysine selection failed');require(np.array_equal(nz.resids,previous['lys_resids']),'Lysine scope drift');atoms=pdb_atoms(pdb);ref={c:np.array([a['xyz'] for a in atoms if a['record']=='ATOM' and a['chain']==c and a['name']=='CA']) for c in 'ABCD'}
 for c in 'ABCD':require(ref[c].shape==ca[c].positions.shape,'Reference CA shape mismatch')
 specs=[]
 for name in ['8RX0_relaxed_G76']:
  pack=np.load(ROOT/'references'/f'{name}.npz');xyz=pack['xyz'];radii=np.array([RAD[e] for e in pack['elements']]);tree=cKDTree(xyz)
  for anchor,chain in [('VHL','C'),('ELOC','B')]:
   ii,jj=matched(str(pack[anchor+'_sequence']),''.join(seq1(r) for r in ca[chain].resnames));specs.append({'name':name+'__'+anchor,'scaffold':name,'anchor':anchor,'chain':chain,'ii':ii,'jj':jj,'ref':pack[anchor+'_CA'][ii],'xyz':xyz,'radii':radii,'tree':tree,'donor':pack['donor_C'],'roles':pack['roles']})
 nf=201;nl=len(nz);ne=len(specs[0]['donor']);distances=np.zeros((len(specs),nf,ne,nl),np.float32);cosines=np.zeros_like(distances);exposure=np.zeros((nf,nl),np.float32);clashes=np.zeros((len(specs),nf),np.float32);fits=np.zeros_like(clashes);times=[];contacterrors=[];occluderids=np.concatenate([heavy[c].indices for c in 'ABCD']+[lig.indices]);occluderrad=np.array([RAD[a.name.lstrip('0123456789')[0]] for c in 'ABCD' for a in heavy[c]]+[RAD[a.name.lstrip('0123456789')[0]] for a in lig]);targetrad=np.array([RAD[a.name.lstrip('0123456789')[0]] for a in heavy['D']]);contacts=pd.read_csv(CAN/'../prosettac_robust/existing_result_audit/contact_core_v1/poses'/poi/pose/'contacts_A.tsv',sep='\t');grid=np.arange(nf)/10;frame_indices=[]
 for index in range(nf):
  frame=int(np.argmin(abs(u.trajectory.times-grid[index]*1000))) if hasattr(u.trajectory,'times') else index
  # Reviewed trajectories are201 frames at0.1ns; exact time is checked, never interpolated.
  ts=u.trajectory[frame];require(abs(ts.time/1000-grid[index])<5.1e-5,'Frame timestamp mismatch');whole(ts);ligwhole(ts);im=assemble(ca['B'].positions,ca['C'].positions,ca['D'].positions,ref['B'],ref['C'],ref['D'],ts.dimensions);lc=ligand_connected_images(heavy['C'].positions,lig.positions,heavy['D'].positions,im['vhl_shift_A'],ts.dimensions);predA,_,_=fit_reference_points(ref['B'],ca['B'].positions,ref['A']);_,ashift,_=image_to_reference(ca['A'].positions,predA,ts.dimensions);shifts={'A':ashift,'B':np.zeros(3),'C':im['vhl_shift_A'],'D':lc['poi_shift_A']};target=heavy['D'].positions+shifts['D'];nzpos=nz.positions+shifts['D'];cepos=ce.positions+shifts['D'];occluders=np.vstack([heavy[c].positions+shifts[c] for c in 'ABCD']+[lig.positions+lc['ligand_shift_A']]);exposure[index]=exposed_fraction(nzpos,nz.indices,occluders,occluderids,occluderrad);error=max(abs(lc['VHL_ligand_min_A']-contacts.true_VHL_heavy_A.iloc[index]),abs(lc['ligand_POI_min_A']-contacts.POI_heavy_A.iloc[index]));contacterrors.append(error);require(error<.002,'Audited contact replay failed')
  for si,spec in enumerate(specs):
   mobile=ca[spec['chain']].positions[spec['jj']]+shifts[spec['chain']];# fit template toMD, then invert toquery fixed template spatial index
   _,rmsd,rot=fit_reference_points(spec['ref'],mobile,spec['ref']);fits[si,index]=rmsd;localtarget=(target-mobile.mean(0))@rot+spec['ref'].mean(0);localnz=(nzpos-mobile.mean(0))@rot+spec['ref'].mean(0);localce=(cepos-mobile.mean(0))@rot+spec['ref'].mean(0);clashes[si,index]=overlap_positive(localtarget,targetrad,spec['tree'],spec['xyz'],spec['radii'])[0];vec=spec['donor'][:,None,:]-localnz[None,:,:];dd=np.linalg.norm(vec,axis=2);side=localnz-localce;cc=np.einsum('eli,li->el',vec,side)/np.maximum(dd*np.linalg.norm(side,axis=1)[None,:],1e-12);distances[si,index]=dd;cosines[si,index]=cc
  times.append(ts.time/1000);frame_indices.append(frame)
 u.trajectory.close();require(max(contacterrors)<.002 and np.array_equal(np.array(times).round(4),grid),'Contact/frame replay incomplete');np.savez_compressed(out/'packet.npz',times_ns=np.array(times),frame_indices=np.array(frame_indices),lys_resids=nz.resids,bridge=previous['true_vhl_mask'],reach50=(previous['corrected_min_e2_A']<50).any(axis=1),nz_exposed_fraction=exposure,donor_distance_A=distances,sidechain_outward_cosine=cosines,target_scaffold_max_overlap_A=clashes,anchor_RMSD_A=fits,spec_names=np.array([s['name'] for s in specs]));r={'poi':poi,'pose':pose,'frames':nf,'lysines':nl,'specs':[s['name'] for s in specs],'endpoint_count':ne,'input_sha256':sources,'source_manifest_sha256':frozen,'packet_sha256':sha256(out/'packet.npz'),'max_contact_replay_error_A':max(contacterrors),'seconds':time.monotonic()-start,'historical_chemistry':True,'catalysis_validated':False};atomic_json(receipt,r);return r

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=4);ap.add_argument('--limit',type=int);args=ap.parse_args();labels=json.loads((ROOT/'config/labels.json').read_text())['labels'];cases=sorted((target,p.parent.name) for target in labels for p in (PRIOR/'analysis/poses'/target).glob('*/receipt.json'));require(len(cases)==141,'Expected complete141 retained core cases');cases=cases[:args.limit] if args.limit else cases;rows=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as ex:
  for r in ex.map(work,cases):rows.append(r);print(len(rows),len(cases),r['poi'],r['pose'],round(r['seconds'],1),flush=True)
 atomic_json(ROOT/('probe_completion.json' if args.limit else 'completion.json'),{'cases':len(rows),'frames':sum(r['frames'] for r in rows),'records':rows,'source_manifest_sha256':verify_source()})
if __name__=='__main__':main()
