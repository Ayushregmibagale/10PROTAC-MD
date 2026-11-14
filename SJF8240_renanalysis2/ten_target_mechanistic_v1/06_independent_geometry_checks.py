"""Independent SVD, brute-force distances, overlap and exposure on preselected frames."""
from pathlib import Path
import sys,json,importlib.util
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('extract',ROOT/'01_extract_matching_inputs.py');e=importlib.util.module_from_spec(s);s.loader.exec_module(e)
import numpy as np
from Bio.SVDSuperimposer import SVDSuperimposer
from scipy.spatial.distance import cdist
from MDAnalysis.transformations import unwrap

def main():
 checked=[];rows=json.loads((ROOT/'extraction_completion.json').read_text());sc=np.load(e.NATIVE/'references/8R5H_observed_SY8.npz');sr=np.array([e.native.RAD[x] for x in sc['elements']])
 for target in ['ABL1','AXL','SIK3_full','RPS6KB2_full']:
  r=sorted([r for r in rows if r['target']==target],key=lambda x:x['pose'])[0];pose=r['pose'];packet=np.load(ROOT/'extracted'/target/pose/'native_packet.npz');sources=r['source_sha256'];tpr=next(p for p in sources if p.endswith('.tpr'));xtc=next(p for p in sources if p.endswith('.xtc'));pdb=next(p for p in sources if p.endswith('.pdb'));u=e.mda.Universe(tpr,xtc,format=e.ReadOnlyXTCReader);uw=unwrap(u.select_atoms('protein'));ligall=u.select_atoms('resname UNL');lw=unwrap(ligall);lig=ligall[ligall.masses>1.5];groups={c:u.select_atoms('moltype Protein_chain_'+c) for c in 'ABCD'};heavy={c:g[g.masses>1.5] for c,g in groups.items()};ca={c:u.select_atoms('moltype Protein_chain_'+c+' and name CA') for c in 'ABCD'};nz=u.select_atoms('moltype Protein_chain_D and resname LYS and name NZ');atoms=e.pdb_atoms(pdb);ref={c:np.array([a['xyz'] for a in atoms if a['record']=='ATOM' and a['chain']==c and a['name']=='CA']) for c in 'ABCD'};ii,jj=e.native.matched(str(sc['VHL_sequence']),''.join(e.seq1(r) for r in ca['C'].resnames));anchor=sc['VHL_CA'][ii];ids=np.concatenate([heavy[c].indices for c in 'ABCD']+[lig.indices]);rad=np.array([e.native.RAD[a.name.lstrip('0123456789')[0]] for c in 'ABCD' for a in heavy[c]]+[e.native.RAD[a.name.lstrip('0123456789')[0]] for a in lig]);targetrad=np.array([e.native.RAD[a.name.lstrip('0123456789')[0]] for a in heavy['D']])
  for frame in [0,100,200]:
   ts=u.trajectory[int(packet['frame_indices'][frame])];uw(ts);lw(ts);im=e.assemble(ca['B'].positions,ca['C'].positions,ca['D'].positions,ref['B'],ref['C'],ref['D'],ts.dimensions);lc=e.ligand_connected_images(heavy['C'].positions,lig.positions,heavy['D'].positions,im['vhl_shift_A'],ts.dimensions);predA,_,_=e.fit_reference_points(ref['B'],ca['B'].positions,ref['A']);_,ashift,_=e.image_to_reference(ca['A'].positions,predA,ts.dimensions);shifts={'A':ashift,'B':np.zeros(3),'C':im['vhl_shift_A'],'D':lc['poi_shift_A']};mobile=ca['C'].positions[jj]+shifts['C'];svd=SVDSuperimposer();svd.set(mobile,anchor);svd.run();rot,translation=svd.get_rotran();donors=sc['donor_C']@rot+translation;nzxyz=nz.positions+shifts['D'];direct=cdist(donors,nzxyz);distance_error=float(np.max(abs(direct-packet['donor_distance_A'][0,frame])));assert distance_error<.002;targetxyz=heavy['D'].positions+shifts['D'];scaffold=sc['xyz']@rot+translation;overlap=max(0.,float((targetrad[:,None]+sr[None,:]-cdist(targetxyz,scaffold)).max()));overlap_error=abs(overlap-packet['target_scaffold_max_overlap_A'][0,frame]);assert overlap_error<.002
   occluders=np.vstack([heavy[c].positions+shifts[c] for c in 'ABCD']+[lig.positions+lc['ligand_shift_A']]);values=[]
   for point,own in zip(nzxyz,nz.indices):
    keep=(ids!=own)&(np.linalg.norm(occluders-point,axis=1)<=6.5);surface=point+e.native.SPHERE*2.95;values.append(float((~(cdist(surface,occluders[keep])<rad[keep][None,:]+1.4).any(1)).mean()))
   exposure_error=float(np.max(abs(np.array(values)-packet['nz_exposed_fraction'][frame])));assert exposure_error<=1/180+1e-6;checked.append({'target':target,'pose':pose,'frame':frame,'NZ_distance_max_error_A':distance_error,'scaffold_overlap_error_A':float(overlap_error),'exposure_max_error':exposure_error,'SVD_fit_RMSD_A':float(svd.get_rms()),'protein_reconstruction':'MDAnalysis unwrap','fit':'Bio.SVDSuperimposer','distance':'scipy.cdist'})
  u.trajectory.close()
 graphchecks=[]
 for target in ['SIK3_full','RPS6KB2_full']:
  native=np.load(ROOT/'inputs/native_features'/(target+'.npz'));zero=np.load(ROOT/'inputs/zero_padded_control'/(target+'.npz'))
  for w in range(4):
   a=native[f'x{w}'];b=zero[f'x{w}'];assert a.shape==b.shape and a.shape[1]==34;assert np.array_equal(a[:,:29],b[:,:29]) and np.all(b[:,29:]==0);assert np.all(a[a[:,23]==0,29:]==0) and np.isfinite(a).all()
   for k in ['edge_index','edge_attr']:assert np.array_equal(native[f'{k}{w}'],zero[f'{k}{w}'])
  graphchecks.append({'target':target,'windows':4,'matching_29_features_and_edges':True,'native_addition_only_on_POI_lysines':True})
 (ROOT/'validation/independent_geometry.json').write_text(json.dumps({'status':'PASS','preselected_targets':['ABL1','AXL','SIK3_full','RPS6KB2_full'],'selection':'first lexicographic retained pose,frames0,100,200;no score-based selection','frame_checks':checked,'graph_checks':graphchecks,'shared_imaging_convention':True,'physical_catalysis_validated':False},indent=2));print('PASS:12 independent geometry frame checks;8 graph-window contracts')
if __name__=='__main__':main()
