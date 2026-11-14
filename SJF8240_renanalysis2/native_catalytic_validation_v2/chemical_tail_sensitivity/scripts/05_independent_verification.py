"""Independent coordinate/geometry replay and target-level metric arithmetic."""
from pathlib import Path
import sys,json,importlib.util
import numpy as np,pandas as pd
from Bio.SVDSuperimposer import SVDSuperimposer
from scipy.spatial.distance import cdist
import MDAnalysis as mda
from MDAnalysis.transformations import unwrap
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'))
spec=importlib.util.spec_from_file_location('extract',Path(__file__).with_name('02_extract_chemical_tail_sensitivity.py'));e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)
def main():
 receipt=json.loads((ROOT/'validation/comparison_receipt.json').read_text())
 for name,h in receipt['outputs_sha256'].items():assert e.sha256(ROOT/'results'/name)==h,name
 checks=[]
 for poi,pose in [('ABL1','pose001'),('DDR2','pose001')]:
  packet=np.load(ROOT/'results/poses'/poi/pose/'packet.npz');u=mda.Universe(str(e.OLD/'random_screen'/poi/pose/'prod/prod.tpr'),str(e.OLD/'analysis/random_screen'/poi/pose/'prod_nojump.xtc'),format=e.ReadOnlyXTCReader);protein=u.select_atoms('protein');uw=unwrap(protein);ligall=u.select_atoms('resname UNL');lw=unwrap(ligall);lig=ligall[ligall.masses>1.5];groups={c:u.select_atoms('moltype Protein_chain_'+c) for c in 'ABCD'};heavy={c:g[g.masses>1.5] for c,g in groups.items()};ca={c:u.select_atoms('moltype Protein_chain_'+c+' and name CA') for c in 'ABCD'};nz=u.select_atoms('moltype Protein_chain_D and resname LYS and name NZ');aa=e.pdb_atoms(e.OLD/'random_screen'/poi/pose/'structure/complex_protein.pdb');ref={c:np.array([a['xyz'] for a in aa if a['record']=='ATOM' and a['chain']==c and a['name']=='CA']) for c in 'ABCD'};sc=np.load(ROOT/'references/8RX0_relaxed_G76.npz');ii,jj=e.matched(str(sc['VHL_sequence']),''.join(e.seq1(r) for r in ca['C'].resnames));anchor=sc['VHL_CA'][ii];ids=np.concatenate([heavy[c].indices for c in 'ABCD']+[lig.indices]);radii=np.array([e.RAD[a.name[0]] for c in 'ABCD' for a in heavy[c]]+[e.RAD[a.name[0]] for a in lig]);tar_radii=np.array([e.RAD[a.name[0]] for a in heavy['D']]);sr=np.array([e.RAD[x] for x in sc['elements']])
  for index in [0,100,200]:
   ts=u.trajectory[int(packet['frame_indices'][index])];uw(ts);lw(ts);im=e.assemble(ca['B'].positions,ca['C'].positions,ca['D'].positions,ref['B'],ref['C'],ref['D'],ts.dimensions);lc=e.ligand_connected_images(heavy['C'].positions,lig.positions,heavy['D'].positions,im['vhl_shift_A'],ts.dimensions);predA,_,_=e.fit_reference_points(ref['B'],ca['B'].positions,ref['A']);_,ashift,_=e.image_to_reference(ca['A'].positions,predA,ts.dimensions);shift={'A':ashift,'B':np.zeros(3),'C':im['vhl_shift_A'],'D':lc['poi_shift_A']};mobile=ca['C'].positions[jj]+shift['C'];svd=SVDSuperimposer();svd.set(mobile,anchor);svd.run();rot,t=svd.get_rotran();donors=sc['donor_C']@rot+t;nzpos=nz.positions+shift['D'];direct=cdist(donors,nzpos);distance_error=float(np.max(abs(direct-packet['donor_distance_A'][0,index])));assert distance_error<.002
   target=heavy['D'].positions+shift['D'];scaffold=sc['xyz']@rot+t;over=tar_radii[:,None]+sr[None,:]-cdist(target,scaffold);overlap=max(0.,float(over.max()));overlap_error=abs(overlap-packet['target_scaffold_max_overlap_A'][0,index]);assert overlap_error<.002
   occluders=np.vstack([heavy[c].positions+shift[c] for c in 'ABCD']+[lig.positions+lc['ligand_shift_A']]);ex=[]
   for point,own in zip(nzpos,nz.indices):
    surface=point+e.SPHERE*2.95;keep=ids!=own;near=keep&(np.linalg.norm(occluders-point,axis=1)<=6.5);hidden=(cdist(surface,occluders[near])<radii[near][None,:]+1.4).any(axis=1);ex.append(float((~hidden).mean()))
   exposure_error=float(np.max(abs(np.array(ex)-packet['nz_exposed_fraction'][index])));assert exposure_error<=1/180+1e-6;checks.append({'poi':poi,'pose':pose,'frame':int(index),'distance_max_error_A':distance_error,'overlap_error_A':float(overlap_error),'exposure_max_error':exposure_error,'protein_reconstruction':'independent MDAnalysis unwrap','rigid_fit':'independent Bio.SVDSuperimposer','distance_engine':'SciPy cdist,not KDTree'})
  u.trajectory.close()
 targets=pd.read_csv(ROOT/'results/primary_target_scores.tsv',sep='\t');metrics=pd.read_csv(ROOT/'results/native_metrics.tsv',sep='\t');auc_checks=[]
 for score in ['bridge_fraction','productive_fraction','joint_bridge_productive_fraction','productive_conditional_on_bridge']:
  yes=targets.loc[targets.label==1,score].dropna().to_numpy();no=targets.loc[targets.label==0,score].dropna().to_numpy();auc=float(np.mean([(a>b)+.5*(a==b) for a in yes for b in no]));expected=metrics[(metrics.reference=='8RX0_relaxed_G76')&(metrics.anchor=='VHL')&(metrics.score==score)].AUROC.iloc[0];assert abs(auc-expected)<1e-12;auc_checks.append({'score':score,'independent_pairwise_AUROC':auc,'independent_positive_negative_pairs':len(yes)*len(no)})
 result={'status':'PASS','frame_replays':checks,'AUROC_replays':auc_checks,'all_primary_outputs_hash_checked':True,'independent_labelled_targets':8,'biological_validation':False,'verification_script_sha256':e.sha256(__file__)};(ROOT/'validation/independent_verification.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS: six raw-frame replays; SASA,distances,clashes and four target AUROCs independently checked.')
if __name__=='__main__':main()
