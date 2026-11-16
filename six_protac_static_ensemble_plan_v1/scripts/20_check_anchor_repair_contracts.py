"""Meaningful positive/negative controls for the frozen anchor repair candidate.

No docking, MD, pipeline source mutation or biological prediction.
"""
import argparse,hashlib,importlib.util,json,sys,shutil,os
from pathlib import Path

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for k in ('adapter','raw','ms143_probe','centroid_params','output'):p.add_argument('--'+k.replace('_','-'),dest=k,type=Path,required=True)
 a=p.parse_args()
 for k in ('adapter','raw','ms143_probe','centroid_params','output'):setattr(a,k,getattr(a,k).resolve())
 a.output.mkdir(parents=True,exist_ok=False)
 if any(c.isspace() for c in str(a.output)):raise ValueError('Use a whitespace-free output directory')
 os.environ['PROTAC_PROSETTAC']=str(a.adapter.parent/'PRosettaC')
 sys.path.insert(0,str(a.adapter));spec=importlib.util.spec_from_file_location('anchor_candidate',a.adapter/'run_docking.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 names,lo,hi=m.read_anchor_contract(a.raw)
 f=m.FlatHarmonicFunc((lo+hi)/2,1.,(hi-lo)/2)
 if any(abs(f.func(x))>1e-9 for x in [lo,(lo+hi)/2,hi]):raise ValueError('Soft restraint penalizes its flat interval')
 if f.func(hi+1)<=0 or f.func(max(0,lo-1))<=0:raise ValueError('Soft restraint fails outside interval')
 bad=a.output/'bad_contract';bad.mkdir()
 for name in ['Protac_params.txt','Init0.pdb','Init1.pdb','anchor_distance_sampling.json']:shutil.copy2(a.raw/name,bad/name)
 pth=bad/'Init0.pdb';pth.write_text('\n'.join(x for x in pth.read_text().splitlines() if not x.startswith('HETATM'))+'\n')
 try:m.read_anchor_contract(bad)
 except ValueError:pass
 else:raise ValueError('Missing named anchor was accepted')
 for name in ['PT0.params','PT1.params','AKT1_PT1_input.pdb']:shutil.copy2(a.ms143_probe/name,a.output/name)
 for head in ('PT0','PT1'):shutil.copy2(a.centroid_params/(head+'.cen.params'),a.output/(head+'.cen.params'))
 import pyrosetta
 from pyrosetta.rosetta.protocols.simple_moves import SwitchResidueTypeSetMover,ReturnSidechainMover
 pyrosetta.init(f'-extra_res_fa {a.output / "PT0.params"} -extra_res_fa {a.output / "PT1.params"} -extra_res_cen {a.output / "PT0.cen.params"} -extra_res_cen {a.output / "PT1.cen.params"} -ignore_unrecognized_res false -constant_seed -jran 42 -out:level 100')
 original=pyrosetta.pose_from_file(str(a.output/'AKT1_PT1_input.pdb'));pose=pyrosetta.Pose(original);restore=ReturnSidechainMover(original);SwitchResidueTypeSetMover('centroid').apply(pose);restore.apply(pose)
 fixed=m.restore_modified_residues(pose,original)
 if len(fixed)!=1 or not pose.residue(fixed[0]).has('P'):raise ValueError('Candidate did not restore exactly the declared phosphothreonine')
 idx=next(i for i in range(1,original.total_residue()+1) if original.residue(i).is_protein() and original.residue(i).name3() not in ('ALA','GLY','THR'))
 ala=next(i for i in range(1,original.total_residue()+1) if original.residue(i).name3()=='ALA')
 broken=pyrosetta.Pose(original);broken.replace_residue(idx,original.residue(ala),True)
 try:m.restore_modified_residues(broken,original)
 except ValueError:pass
 else:raise ValueError('Unexpected protein residue substitution was silently accepted')
 report={'status':'REPAIR_CONTRACT_CONTROLS_PASS','checks':['zero penalty throughout sampled flat interval','positive restraint penalty outside interval','missing anchor rejected','actual candidate restores phosphorylated Thr after centroid roundtrip','unexpected protein residue substitution rejected'],'candidate_source_sha256':hashlib.sha256((a.adapter/'run_docking.py').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Static execution contracts only; not reference recovery qualification'}
 (a.output/'REPAIR_CONTRACT_CHECK.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
