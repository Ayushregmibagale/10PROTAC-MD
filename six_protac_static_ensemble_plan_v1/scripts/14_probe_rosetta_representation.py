"""Bounded diagnostic: record residue representation around inherited low-res docking.

Uses first global pose by index, never selected from native recovery. No pipeline
source is changed. Does not establish the cause of failed structural recovery.
"""
import argparse,json,hashlib,os
from pathlib import Path

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 a.raw=a.raw.resolve();a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
 # Restrict the diagnostic process to at most four available cores.
 os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:4])
 import pyrosetta
 from pyrosetta.rosetta.protocols.docking import DockingLowRes,setup_foldtree
 from pyrosetta.rosetta.utility import vector1_int
 pyrosetta.init(f'-extra_res_fa {a.raw / "PT0.params"} -extra_res_fa {a.raw / "PT1.params"} -ignore_unrecognized_res false -constant_seed -jran 42 -out:level 100')
 source=a.raw/'Patchdock_Results/pd.1.pdb';pose=pyrosetta.pose_from_file(str(source));jumps=vector1_int();jumps.append(1);setup_foldtree(pose,'ABCX_DY',jumps)
 def state():
  residues=[r for r in pose.residues if r.is_protein()]
  return {'pose_is_fullatom':pose.is_fullatom(),'protein_residues':len(residues),'protein_residues_with_CEN_atom':sum(r.has('CEN') for r in residues),'protein_atom_count':sum(r.natoms() for r in residues)}
 before=state();score=pyrosetta.create_score_function('interchain_cen');mover=DockingLowRes(score,1)
 doc=DockingLowRes.__doc__
 mover.apply(pose);after=state()
 report={'status':'REPRESENTATION_DIAGNOSTIC_COMPLETE','before_DockingLowRes':before,'after_DockingLowRes':after,'source_pose_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'source_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'runtime_version':pyrosetta.version(),'DockingLowRes_doc':doc,'limitations':['Single first-index pose, seed42; diagnostic only','No changes to frozen baseline or claim that this alone explains failed recovery']}
 (a.output/'REPRESENTATION_PROBE.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
