"""Check MS143 AKT phosphorylated-Thr308 survives explicit centroid roundtrip.

No docking, minimization or MD. This is not binding/production qualification.
"""
import argparse,hashlib,json,math,os,subprocess,sys,shutil
from pathlib import Path

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for k in ('inputs','generator','output'):p.add_argument('--'+k,type=Path,required=True)
 p.add_argument('--restore-declared-tpo',action='store_true',help='Separate candidate: restore declared original full-atom TPO aligned to its current backbone')
 p.add_argument('--rigid-perturb-centroid',action='store_true',help='Apply a fixed90degree rotation and translation before restoration; static geometry control')
 a=p.parse_args()
 for k in ('inputs','generator','output'):setattr(a,k,getattr(a,k).resolve())
 a.output.mkdir(parents=True,exist_ok=False);os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:4])
 if any(c.isspace() for c in str(a.output)):
  raise ValueError('PyRosetta parameter directory must have no whitespace')
 for head in ('PT0','PT1'):
  cmd=[sys.executable,str(a.generator),str(a.inputs/(head+'_explicit_H.sdf')),'-n',head,'-p',str(a.output/head),'--centroid','--chain','X']
  r=subprocess.run(cmd,text=True,capture_output=True,timeout=120);(a.output/(head+'_generation.log')).write_text(r.stdout+r.stderr)
  if r.returncode:raise RuntimeError('Centroid head typing failed')
  # Preserve the original full-atom parameter definition exactly.
  shutil.copy2(a.inputs/(head+'.params'),a.output/(head+'.original_fa.params'))
 import pyrosetta
 from pyrosetta.rosetta.protocols.simple_moves import SwitchResidueTypeSetMover,ReturnSidechainMover
 pyrosetta.init(f'-extra_res_fa {a.output / "PT0.original_fa.params"} -extra_res_fa {a.output / "PT1.original_fa.params"} -extra_res_cen {a.output / "PT0.cen.params"} -extra_res_cen {a.output / "PT1.cen.params"} -ignore_unrecognized_res false -constant_seed -jran 42 -out:level 100')
 source=a.inputs/'AKT1_PT1_input.pdb';local_source=a.output/'AKT1_PT1_input.pdb';shutil.copy2(source,local_source)
 pose=pyrosetta.pose_from_file(str(local_source));original=pyrosetta.Pose(pose)
 idx=[i for i in range(1,pose.total_residue()+1) if pose.pdb_info().chain(i)=='D' and pose.pdb_info().number(i)==308]
 if len(idx)!=1:raise ValueError('Ambiguous/missing Thr308')
 i=idx[0]
 def state():
  r=pose.residue(i);return {'is_fullatom':pose.is_fullatom(),'residue_name':r.name(),'has_phosphate_P':r.has('P'),'residue_atoms':[r.atom_name(j).strip() for j in range(1,r.natoms()+1)]}
 before=state();restore=ReturnSidechainMover(original);SwitchResidueTypeSetMover('centroid').apply(pose);centroid=state()
 def expected(xyz):
  if not a.rigid_perturb_centroid:return xyz
  return pyrosetta.rosetta.numeric.xyzVector_double_t(-xyz.y+13,xyz.x-7,xyz.z+5)
 if a.rigid_perturb_centroid:
  for ri in range(1,pose.total_residue()+1):
   for ai in range(1,pose.residue(ri).natoms()+1):
    pose.set_xyz(pyrosetta.rosetta.core.id.AtomID(ai,ri),expected(pose.residue(ri).xyz(ai)))
 restore.apply(pose);uncorrected=state()
 if a.restore_declared_tpo:
  pose.replace_residue(i,original.residue(i),True)
 after=state()
 check={'status':'MODIFICATION_ROUNDTRIP_MEASURED','before':before,'centroid':centroid,'uncorrected_restore':uncorrected,'restore_declared_tpo':a.restore_declared_tpo,'rigid_perturb_centroid':a.rigid_perturb_centroid,'after':after,'source_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest()}
 if not pose.is_fullatom() or not after['has_phosphate_P'] or set(before['residue_atoms'])!=set(after['residue_atoms']):
  check['status']='FAILED_PHOSPHOTHREONINE_NOT_RESTORED'
  (a.output/'MS143_CENTROID_MODIFICATION_CHECK.json').write_text(json.dumps(check,indent=2)+'\n')
  raise ValueError('Phosphorylated Thr308 not restored intact')
 shifts=[(pose.residue(i).xyz(name)-expected(original.residue(i).xyz(name))).norm() for name in before['residue_atoms']]
 score=float(pyrosetta.create_score_function('ref2015')(pose))
 if not math.isfinite(score) or max(shifts)>.003:raise ValueError('Unexpected coordinate/score roundtrip failure')
 pose.dump_pdb(str(a.output/'AKT1_PT1_roundtrip.pdb'))
 report={'status':'MS143_PHOSPHOTHREONINE_CENTROID_ROUNDTRIP_PASS','before':before,'centroid':centroid,'uncorrected_restore':uncorrected,'restore_declared_tpo':a.restore_declared_tpo,'rigid_perturb_centroid':a.rigid_perturb_centroid,'after':after,'max_TPO_atom_displacement_from_expected_frame_A':max(shifts),'fullatom_score_finite':True,'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'source_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'limitations':['No docking or minimization; static residue-type roundtrip and optional fixed rigid transform only','Does not qualify AKT binding pose, complete E3 assembly or reference recovery']}
 (a.output/'MS143_CENTROID_MODIFICATION_CHECK.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
