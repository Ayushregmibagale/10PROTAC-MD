"""Recover per-frame fraction of E2 ensemble points within50A, not a min-distance proxy."""
from pathlib import Path
import sys,os,json,concurrent.futures,importlib.util
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[k]='1'
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('extract',ROOT/'01_extract_matching_inputs.py');e=importlib.util.module_from_spec(s);s.loader.exec_module(e)
import numpy as np

def work(r):
 t,pose=r['poi'],r['pose'];out=ROOT/'reach50'/t/pose;out.mkdir(parents=True,exist_ok=True);receipt=out/'receipt.json'
 if receipt.exists():
  rr=json.loads(receipt.read_text());assert e.sha(out/'packet.npz')==rr['packet_sha256'] and rr['script_sha256']==e.sha(__file__);return rr
 sources=r['input_sha256'];tpr=next(p for p in sources if p.endswith('.tpr'));xtc=next(p for p in sources if p.endswith('.xtc'));pdb=next(p for p in sources if p.endswith('.pdb'))
 for p in [tpr,xtc,pdb]:assert e.sha(p)==sources[p]
 oldpath=e.CAN/'full_frame_e2_delta50_v1/analysis/poses'/t/pose/'packet.npz';old=np.load(oldpath);u=e.mda.Universe(tpr,xtc,format=e.ReadOnlyXTCReader);whole=e.WholeProteinTree(u.select_atoms('protein'));ligall=u.select_atoms('resname UNL');lw=e.WholeProteinTree(ligall);lig=ligall[ligall.masses>1.5];ca={c:u.select_atoms('moltype Protein_chain_'+c+' and name CA') for c in 'BCD'};heavy={c:u.select_atoms('moltype Protein_chain_'+c+' and mass 1.5:1000') for c in 'CD'};nz=u.select_atoms('moltype Protein_chain_D and resname LYS and name NZ');atoms=e.pdb_atoms(pdb);ref={c:np.array([a['xyz'] for a in atoms if a['record']=='ATOM' and a['chain']==c and a['name']=='CA']) for c in 'BCD'};ensemble_path=Path('/home/ayushregmibagale/protac-md/reference_structures/e2_c85_ensemble.json');ensemble=np.array(next(iter(json.loads(ensemble_path.read_text()).values())),float);values=[];maximum=0.;bridge=[];time=[]
 for i in range(201):
  ts=u.trajectory[i];assert abs(ts.time/1000-i/10)<5.1e-5;whole(ts);lw(ts);im=e.assemble(ca['B'].positions,ca['C'].positions,ca['D'].positions,ref['B'],ref['C'],ref['D'],ts.dimensions);lc=e.ligand_connected_images(heavy['C'].positions,lig.positions,heavy['D'].positions,im['vhl_shift_A'],ts.dimensions);points,_,_=e.fit_reference_points(ref['B'],ca['B'].positions,ensemble);d=e.distance_array(nz.positions+lc['poi_shift_A'],points);maximum=max(maximum,float(np.max(abs(d.min(1)-old['corrected_min_e2_A'][i]))));values.append((d<50).mean(1));bridge.append(lc['VHL_ligand_min_A']<5 and lc['ligand_POI_min_A']<5);time.append(ts.time/1000)
 u.trajectory.close();assert maximum<.002 and np.array_equal(bridge,old['true_vhl_mask']);np.savez_compressed(out/'packet.npz',times_ns=time,lys_resids=nz.resids,reach50_fraction=values,bridge=bridge);rr={'target':t,'pose':pose,'packet_sha256':e.sha(out/'packet.npz'),'script_sha256':e.sha(__file__),'source_sha256':{**{p:sources[p] for p in [tpr,xtc,pdb]},str(oldpath):e.sha(oldpath),str(ensemble_path):e.sha(ensemble_path)},'max_min_distance_replay_error_A':maximum,'ensemble_points':len(ensemble)};e.save(receipt,rr);return rr

def main():
 todo=json.loads((e.NATIVE/'completion.json').read_text())['records'];rows=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=4) as ex:
  for rr in ex.map(work,todo):rows.append(rr);print(len(rows),rr['target'],rr['pose'],flush=True)
 e.save(ROOT/'reach50_completion.json',rows)
if __name__=='__main__':main()
