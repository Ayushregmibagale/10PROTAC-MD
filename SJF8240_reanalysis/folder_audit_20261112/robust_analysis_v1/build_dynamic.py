"""Generic periodic joint-contact extraction with sequence-verified role selections.
No docking/MD. Explicit per-compound topology/ligand/time-grid configuration.
"""
from pathlib import Path
import json,argparse,os,tempfile
import numpy as np
from contracts import *
def write_npz(path,**arrays):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
 with tempfile.NamedTemporaryFile(dir=path.parent,delete=False) as f:name=f.name;np.savez_compressed(f,**arrays)
 os.replace(name,path)
def sequence(group):
 from Bio.SeqUtils import seq1
 aliases={'HID':'H','HIE':'H','HIP':'H','HSD':'H','HSE':'H','HSP':'H','CYX':'C','MSE':'M'}
 result=''.join(seq1(str(n),custom_map=aliases,undef_code='X') for n in group.resnames);require('X' not in result,'Unknown residue identity');return result

def select_roles(universe,case,references):
 groups={}
 for role,entry in case['protein_roles'].items():
  group=universe.select_atoms(entry['selection']);require(len(group)>0,'Empty protein selection '+role);require(np.all(group.names=='CA') and len(group)==len(group.residues),'Select exactly one CA per residue')
  require(sequence(group)==entry['model_sequence'],'Role sequence mismatch '+role)
  verify_sequence_identity(sequence(group),entry['uniprot'],references)
  if role in ['vhl','elob']:require(entry['uniprot']=={'vhl':'P40337','elob':'Q15370'}[role],'Wrong protein accession for role')
  groups[role]=group
 require(set(groups)=={'vhl','elob','poi'},'Explicit VHL/ELOB/POI roles required')
 for a in groups:
  for b in groups:
   if a!=b:require(not set(groups[a].indices)&set(groups[b].indices),'Overlapping protein selections')
 return groups

def run(manifest,output):
 import MDAnalysis as mda
 from MDAnalysis.coordinates.XTC import XTCReader
 from MDAnalysis.lib.distances import distance_array
 class ReadOnlyXTCReader(XTCReader):
  def _load_offsets(self):self._read_offsets(store=False)
  def _read_offsets(self,store=False):return super()._read_offsets(store=False)
 manifest=Path(manifest).resolve();base=manifest.parent;c=json.loads(manifest.read_text());output=Path(output).resolve();require(not output.exists(),'Fresh output required; old extraction is preserved')
 require(c['units']=={'coordinates':'angstrom','time':'ns'} and c['raw_trajectory_time_unit']=='ps','Explicit supported units required');require(bool(c['compound_id']),'Compound identity missing')
 grid=np.asarray(c['expected_times_ns'],float);require(grid.ndim==1 and len(grid)>1 and np.isfinite(grid).all() and np.all(np.diff(grid)>0),'Strictly increasing expected time grid required');threshold=c['contact_A'];edge_min=c['edge_frequency_min'];require(threshold>0 and 0<edge_min<=1,'Invalid thresholds')
 reference_path=bound_file(base,c['reference_sequences']);references=json.loads(reference_path.read_text())['sequences']
 cases=c['cases'];require(cases and len({(x['target_id'],x['pose_id']) for x in cases})==len(cases),'Unique nonempty cases required');bound={}
 for case in cases:
  require(case['protein_roles']['poi']['uniprot']==c['target_uniprot'][case['target_id']],'POI accession mismatch')
  for k in ['tpr','xtc']:bound[(case['target_id'],case['pose_id'],k)]=bound_file(base,case[k])
  require(all(not p.is_relative_to(output) for p in bound.values()),'Output overlaps source')
 output.mkdir(parents=True);frozen={'manifest':sha(manifest),'script':sha(__file__),'contracts':sha(Path(__file__).with_name('contracts.py'))};atomic_json(output/'intent.json',frozen);records=[];aggregation={}
 for case in cases:
  key=(case['target_id'],case['pose_id']);u=mda.Universe(str(bound[(*key,'tpr')]),str(bound[(*key,'xtc')]),format=ReadOnlyXTCReader);groups=select_roles(u,case,references);lig=u.select_atoms(case['ligand_selection']);lig=lig[lig.masses>1.5];require(len(lig)==case['ligand_heavy_atoms'] and len(lig)>0,'Ligand atom count mismatch');require(not any(set(lig.indices)&set(g.indices) for g in groups.values()),'Ligand overlaps protein')
  sums={r:{'joint':np.zeros((len(groups[r]),len(groups['poi']))),'left':np.zeros(len(groups[r])),'right':np.zeros(len(groups['poi']))} for r in ['elob','vhl']};times=[];last=None;duplicates=0
  for ts in u.trajectory:
   t=float(ts.time)/1000
   if t<grid[0]-1e-6:continue
   if t>grid[-1]+1e-6:break
   if times and abs(t-times[-1])<1e-8:
    require(c['duplicate_policy']=='identical_only' and np.array_equal(last,ts.positions),'Conflicting or disallowed duplicate frame');duplicates+=1;continue
   require(not times or t>times[-1],'Trajectory time order')
   box=np.asarray(ts.dimensions,float);require(box.shape==(6,) and np.isfinite(box).all() and np.all(box[:3]>0) and np.all((box[3:]>0)&(box[3:]<180)),'Invalid periodic box')
   from MDAnalysis.lib.mdamath import triclinic_vectors
   require(np.linalg.det(triclinic_vectors(box,dtype=np.float64))>1e-6,'Degenerate box')
   right=distance_array(lig.positions,groups['poi'].positions,box=box).min(0)
   for role,a in sums.items():
    left=distance_array(lig.positions,groups[role].positions,box=box).min(0);a['joint']+=joint_contact_observation(left,right,threshold);a['left']+=left;a['right']+=right
   times.append(t);last=ts.positions.copy()
  u.trajectory.close();require(np.asarray(times).shape==grid.shape and np.allclose(times,grid,atol=5.1e-5,rtol=0),'Incomplete or unexpected time grid')
  for role,a in sums.items():
   for k in a:a[k]/=len(times)
   packet={**a,'left_resids':groups[role].resids,'right_resids':groups['poi'].resids,'left_resnames':groups[role].resnames.astype('U4'),'right_resnames':groups['poi'].resnames.astype('U4')};p=output/'poses'/case['target_id']/case['pose_id']/(role+'.npz');write_npz(p,**packet);aggregation.setdefault((case['target_id'],role),[]).append(packet)
  records.append({'target':case['target_id'],'pose':case['pose_id'],'frames':len(times),'identical_duplicate_frames':duplicates,'input_sha256':{k:sha(bound[(*key,k)]) for k in ['tpr','xtc']},'role_sequences':{r:sequence(g) for r,g in groups.items()}})
  atomic_json(output/'progress.json',{'completed_cases':len(records),'expected_cases':len(cases),'last_case':list(key)});print(json.dumps({'completed':len(records),'expected':len(cases),'target':key[0],'pose':key[1]}),flush=True)
 for (target,role),packets in aggregation.items():
  ref=packets[0]
  for p in packets:
   require(all(np.array_equal(ref[k],p[k]) for k in ['left_resids','right_resids','left_resnames','right_resnames']),'Pose residue correspondence differs')
  joint=np.mean([p['joint'] for p in packets],0);left=np.mean([p['left'] for p in packets],0);right=np.mean([p['right'] for p in packets],0);nl=len(left);x=np.zeros((nl+len(right),3),np.float32)
  for offset,ids,ri in [(0,ref['left_resids'],0),(nl,ref['right_resids'],1)]:x[offset:offset+len(ids),ri]=1;x[offset:offset+len(ids),2]=(ids-ids.min())/max(int(ids.max()-ids.min()),1)
  i,j=np.nonzero(joint>=edge_min);ei=np.array([np.r_[i,nl+j],np.r_[nl+j,i]],np.int64);ea=np.column_stack([joint[i,j],left[i]/10,right[j]/10]).astype(np.float32);path=output/'targets'/target/role/'graph.npz';write_npz(path,x=x,edge_index=ei,edge_attr=np.concatenate([ea,ea]),compound_id=np.array([c['compound_id']]),poi=np.array([target]),role=np.array([role]));graph_arrays(path,3);atomic_json(path.with_suffix('.json'),{'status':'COMPLETE','poses':len(packets),'graph_sha256':sha(path),'role':role,'compound_id':c['compound_id'],'target_id':target,'equal_pose_weight':True,'joint_not_product_of_marginals':True})
 require(sha(reference_path)==c['reference_sequences']['sha256'],'Reference sequence drift')
 require(sha(manifest)==frozen['manifest'] and sha(__file__)==frozen['script'] and sha(Path(__file__).with_name('contracts.py'))==frozen['contracts'],'Source changed')
 for case in cases:
  for k in ['tpr','xtc']:require(sha(bound[(case['target_id'],case['pose_id'],k)])==case[k]['sha256'],'Raw input drift')
 atomic_json(output/'completion.json',{'status':'COMPLETE','cases':len(records),'records':records,'source':frozen,'new_MD':False})
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--manifest',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.manifest,a.output)
