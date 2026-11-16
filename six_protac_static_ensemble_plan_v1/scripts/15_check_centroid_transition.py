"""Separate candidate smoke check for explicit centroid conversion and restoration.

Not a replacement docking pipeline or recovery benchmark. Uses global index1.
"""
import argparse,hashlib,json,os,subprocess,sys
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser(description=__doc__)
 for name in ('raw','generator','output'):p.add_argument('--'+name,type=Path,required=True)
 a=p.parse_args()
 for k in ('raw','generator','output'):setattr(a,k,getattr(a,k).resolve())
 a.output.mkdir(parents=True,exist_ok=False);os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:4]);receipts=[]
 for head in ('PT0','PT1'):
  cmd=[sys.executable,str(a.generator),str(a.raw/(head+'_explicit_H.sdf')),'-n',head,'-p',str(a.output/head),'--centroid','--chain','X']
  result=subprocess.run(cmd,text=True,capture_output=True,timeout=120);(a.output/(head+'_generation.log')).write_text(result.stdout+result.stderr)
  if result.returncode:raise RuntimeError('Centroid parameter generation failed')
  receipts.append({'head':head,'source_sdf_sha256':sha(a.raw/(head+'_explicit_H.sdf')),'centroid_params_sha256':sha(a.output/(head+'.cen.params'))})
 import pyrosetta
 from pyrosetta.rosetta.protocols.simple_moves import SwitchResidueTypeSetMover,ReturnSidechainMover
 from pyrosetta.rosetta.protocols.docking import DockingLowRes,setup_foldtree
 from pyrosetta.rosetta.utility import vector1_int
 pyrosetta.init(f'-extra_res_fa {a.raw / "PT0.params"} -extra_res_fa {a.raw / "PT1.params"} -extra_res_cen {a.output / "PT0.cen.params"} -extra_res_cen {a.output / "PT1.cen.params"} -ignore_unrecognized_res false -constant_seed -jran 42 -out:level 100')
 source=a.raw/'Patchdock_Results/pd.1.pdb';pose=pyrosetta.pose_from_file(str(source));jumps=vector1_int();jumps.append(1);setup_foldtree(pose,'ABCX_DY',jumps)
 def state():
  rs=[r for r in pose.residues if r.is_protein()]
  return {'is_fullatom':pose.is_fullatom(),'protein_residues':len(rs),'protein_residues_with_CEN':sum(r.has('CEN') for r in rs),'nonprotein_residues':[r.name() for r in pose.residues if not r.is_protein()]}
 original=pyrosetta.Pose(pose);restore=ReturnSidechainMover(original);before=state();SwitchResidueTypeSetMover('centroid').apply(pose);centroid=state()
 if pose.is_fullatom() or centroid['protein_residues_with_CEN']<1:raise ValueError('Failed centroid conversion')
 DockingLowRes(pyrosetta.create_score_function('interchain_cen'),1).apply(pose);after_low=state();restore.apply(pose);after_restore=state()
 if not pose.is_fullatom() or after_restore['protein_residues_with_CEN']!=0:raise ValueError('Failed full-atom restoration')
 pose.dump_pdb(str(a.output/'restored_after_lowres.pdb'))
 r={'status':'CENTROID_TRANSITION_SMOKE_PASS','before':before,'centroid':centroid,'after_lowres':after_low,'after_restore':after_restore,'parameter_receipts':receipts,'source_sha256':sha(source),'generator_sha256':sha(a.generator),'script_sha256':sha(Path(__file__)),'limitations':['Single first-index pose, seed42; no full refinement, linker reconstruction or recovery qualification','Head geometry/sterics and phosphorylated MS143 compatibility remain unqualified','Frozen baseline has not been changed']}
 (a.output/'CENTROID_TRANSITION_CHECK.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
if __name__=='__main__':main()
