"""Bounded retry-budget diagnostic on first three failed, in-range GSK geometries.

Selection uses sampled linker range and rejection receipts, never native RMSD or
experimental labels. Does not change the running qualification protocol.
"""
import argparse,csv,hashlib,json,os,sys,time
from pathlib import Path
import numpy as np
from rdkit import Chem

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for k in ('adapter','raw','drift','output'):p.add_argument('--'+k,type=Path,required=True)
 a=p.parse_args()
 for k in ('adapter','raw','drift','output'):setattr(a,k,getattr(a,k).resolve())
 a.output.mkdir(parents=True,exist_ok=False);os.sched_setaffinity(0,[min(os.sched_getaffinity(0))]);sys.path.insert(0,str(a.adapter))
 from constrained_linker import constrained_linker,stereo_key
 from official_parameters import docked_heads
 heads_paths=[a.raw/'vhl_head_H.sdf',a.raw/'poi_head_H.sdf'];heads=[Chem.SDMolSupplier(str(x),removeHs=True)[0] for x in heads_paths];full=Chem.MolFromSmiles((a.raw/'full.smi').read_text().strip())
 eligible=[r for r in csv.DictReader(a.drift.open(),delimiter='\t') if r['reference']=='7PI4_GSK215' and r['outside_sampled_range']=='False' and r['linker_status']=='GEOMETRY_SAMPLING_REJECTED']
 selected=sorted(eligible,key=lambda r:(int(r['global_index']),int(r['local_index'])))[:3]
 if len(selected)!=3:raise ValueError('Need three eligible geometries')
 rows=[];start=time.monotonic()
 for row in selected:
  name=f"pd.{int(row['global_index'])}_docking_{int(row['local_index']):04d}.pdb";pose=a.raw/'Patchdock_Results'/name
  combined=docked_heads(pose,heads_paths,root=a.raw);positions=combined.GetConformer().GetPositions();cut=heads[0].GetNumAtoms();targets=[positions[:cut],positions[cut:]]
  for budget in (20,100):
   if time.monotonic()-start>150:raise TimeoutError('Bounded diagnostic elapsed budget reached; partial report retained')
   before=time.monotonic();mol,result=constrained_linker(full,heads,targets,42,attempts=budget,embed_iterations=20)
   if mol is not None:
    if stereo_key(mol)!=Chem.MolToInchiKey(full):raise ValueError('Recovered linker stereo mismatch')
    with Chem.SDWriter(str(a.output/(name+f'.retries{budget}.sdf'))) as w:w.write(mol)
   rows.append({'pose':name,'pose_sha256':hashlib.sha256(pose.read_bytes()).hexdigest(),'outer_retry_budget':budget,'internal_iteration_limit':20,'seed':42,'accepted':mol is not None,'seconds':time.monotonic()-before,'diagnostic':result})
   report={'status':'IN_PROGRESS','selection':'First3 by global/local index among failed GSK geometries within preexisting sampled anchor range; no native/label selection','rows':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'engine_sha256':hashlib.sha256((a.adapter/'constrained_linker.py').read_bytes()).hexdigest()};(a.output/'EMBEDDING_RETRY_PROBE.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:rows[-1][k] for k in ['pose','outer_retry_budget','accepted','seconds']}),flush=True)
 report['status']='BOUNDED_RETRY_PROBE_COMPLETE';report['limitations']=['Three selected failed geometries and one seed only; not a convergence study','In-range anchor distance does not guarantee feasible head orientations','No protein-clash or docking recovery qualification'];(a.output/'EMBEDDING_RETRY_PROBE.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
