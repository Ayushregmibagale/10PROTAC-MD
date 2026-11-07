"""Audit retained E2 ensemble multiplicity and static frame correspondence."""
from pathlib import Path
import json,sys,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import pdb_atoms,sha256,atomic_json,ValidationError
project=Path('/home/ayushregmibagale/protac-md');ensemble=project/'reference_structures/e2_c85_ensemble.json';data=json.loads(ensemble.read_text());arrays={k:np.asarray(v,dtype=float) for k,v in data.items()}
if any(a.ndim!=2 or a.shape[1]!=3 or not np.isfinite(a).all() for a in arrays.values()):raise ValidationError('Invalid retained E2 coordinates')
array_hashes={k:hashlib.sha256(a.tobytes()).hexdigest() for k,a in arrays.items()}
base=project/'simulations/DDR2/structure/rank1/complex_protein.pdb'
def ca(path):
 atoms=[a for a in pdb_atoms(path) if a['record']=='ATOM' and a['chain']=='B' and a['name']=='CA']
 if len(atoms)!=85:raise ValidationError('Expected 85 Elongin-C Cα atoms: '+str(path))
 return np.asarray([a['xyz'] for a in atoms]),[(a['resnum'],a['icode'],a['resname']) for a in atoms]
reference,keys=ca(base);roster=json.loads((ROOT/'existing_result_audit/bridge_cutoff_v2_core_20ns/pose_provenance.json').read_text());rows=[];sources={str(base):sha256(base)}
for poi in sorted({r['poi'] for r in roster}):
 main=project/'simulations'/poi/'structure/rank1/complex_protein.pdb'
 if poi!='EPHB2' and main.exists():selected=main
 else:
  selected=next((p for p in sorted((project/'random_screen'/poi).glob('pose*/structure/complex_protein.pdb')) if len([a for a in pdb_atoms(p) if a['name']=='CA' and a['chain']=='B'])==85),None)
 if selected is None:raise ValidationError('No retained E2 reference PDB for '+poi)
 xyz,identities=ca(selected)
 if identities!=keys:raise ValidationError('Elongin-C residue correspondence differs: '+poi)
 raw=float(np.sqrt(np.mean(np.sum((xyz-reference)**2,axis=1))));rows.append({'poi':poi,'producer_reference':str(selected),'raw_CA_RMSD_to_common_reference_A':raw});sources[str(selected)]=sha256(selected)
for row in roster:
 p=project/'random_screen'/row['poi']/row['pose_id']/'structure/complex_protein.pdb';xyz,identities=ca(p)
 if identities!=keys:raise ValidationError('Pose reference sequence correspondence differs')
 row['ELOC_CA_RMSD_A']=float(np.sqrt(np.mean(np.sum((xyz-reference)**2,axis=1))));sources[str(p)]=sha256(p)
out=ROOT/'existing_result_audit/e2_reference_frame_v1';out.mkdir(exist_ok=False)
report={'status':'STATIC_E2_REFERENCE_AUDIT_NOT_DYNAMIC_IMAGING_VALIDATION','ensemble_keys':len(arrays),'points_per_key':{k:len(a) for k,a in arrays.items()},'flattened_point_records':sum(len(a) for a in arrays.values()),'unique_ensemble_array_hashes':len(set(array_hashes.values())),'flattened_point_multiplicity':len(arrays) if len(set(array_hashes.values()))==1 else None,'unique_point_coordinates':len(np.unique(np.vstack(list(arrays.values())),axis=0)),'per_key_array_sha256':array_hashes,'producer_references':rows,'random_pose_structures_checked':len(roster),'max_random_pose_ELOC_CA_RMSD_to_common_reference_A':max(r['ELOC_CA_RMSD_A'] for r in roster),'ensemble_sha256':sha256(ensemble),'source_PDB_sha256':sources,'audit_source_sha256':sha256(Path(__file__)),'original_files_modified':False,'new_simulations_executed':False,'limitations':['Equal duplication of identical coordinate arrays leaves minima and equally weighted fractions unchanged, but is not independent ensemble diversity.','This checks static correspondence of the retained reference frames; dynamic unwrapping/coherent imaging and sampling-prior validity are separate questions.','The source ensemble uses a Gaussian positional perturbation model; its point fraction is not an experimentally calibrated reaction probability.']}
atomic_json(out/'summary.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ['source_PDB_sha256','per_key_array_sha256','producer_references','limitations']},indent=2))
