"""Bounded, read-only imaging probe on three preselected retained MD cases.

Uses the old E2 rotation function without executing its original module. Virtual
E2 coordinates are not reinterpreted as periodic simulated particles.
"""
from pathlib import Path
import ast,json,sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis import transformations
from MDAnalysis.lib.distances import distance_array
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import pdb_atoms,sha256,atomic_json,ValidationError
from read_only_trajectory import ReadOnlyXTCReader
project=Path('/home/ayushregmibagale/protac-md');old=project/'scripts/02_analysis/16_all_poi_e2_reachability.py';tree=ast.parse(old.read_text());function=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='kabsch_rotation');namespace={'np':np};exec(compile(ast.Module(body=[function],type_ignores=[]),str(old),'exec'),namespace);rotation=namespace['kabsch_rotation']
ref=project/'simulations/DDR2/structure/rank1/complex_protein.pdb';reference=np.asarray([a['xyz'] for a in pdb_atoms(ref) if a['record']=='ATOM' and a['chain']=='B' and a['name']=='CA'],dtype=np.float32)
if reference.shape!=(85,3):raise ValidationError('Expected common Elongin-C reference')
static=json.loads((ROOT/'existing_result_audit/e2_reference_frame_v1/summary.json').read_text())
if any(r['raw_CA_RMSD_to_common_reference_A']>1e-6 for r in static['producer_references']):raise ValidationError('Producer reference frames differ; review before applying common reference')
ensemble_file=project/'reference_structures/e2_c85_ensemble.json';ensemble=np.asarray(json.loads(ensemble_file.read_text())['DDR2_rank1'],dtype=float)
plan=[('DDR2','pose001'),('MET','pose002'),('AXL','pose011')];frames=[0,50,100,150,200];out=ROOT/'existing_result_audit/e2_dynamic_imaging_probe_v1'
if out.exists():raise ValidationError('Fresh probe output required')
source_hash=sha256(Path(__file__));inputs={str(old):sha256(old),str(ref):sha256(ref),str(ensemble_file):sha256(ensemble_file)};rows=[]
for poi,pose in plan:
 tpr=project/'random_screen'/poi/pose/'prod/prod.tpr';xtc=project/'analysis/random_screen'/poi/pose/'prod_nojump.xtc'
 inputs.update({str(tpr):sha256(tpr),str(xtc):sha256(xtc)})
 u=mda.Universe(str(tpr),str(xtc),format=ReadOnlyXTCReader);protein=u.select_atoms('protein');u.trajectory.add_transformations(transformations.unwrap(protein))
 anchor=u.select_atoms('moltype Protein_chain_B and name CA');vhl=u.select_atoms('moltype Protein_chain_C');target=u.select_atoms('moltype Protein_chain_D');eloc=u.select_atoms('moltype Protein_chain_B');lys=target.select_atoms('resname LYS and name NZ')
 if len(anchor)!=85 or not len(lys) or len(u.trajectory)<201:raise ValidationError('Probe topology/frame correspondence failed')
 heavy=lambda g:g[g.masses>1.5]
 for frame in frames:
  ts=u.trajectory[frame];time=float(ts.time)/1000
  if abs(time-frame*.1)>5.1e-6:raise ValidationError('Unexpected frame time; no truncation/interpolation')
  mobile=anchor.positions;pdb_c=reference.mean(0);mob_c=mobile.mean(0);R=rotation(reference-pdb_c,mobile-mob_c);fitted=(reference-pdb_c)@R.T+mob_c;e2=(ensemble-pdb_c)@R.T+mob_c
  row={'poi':poi,'pose_id':pose,'frame':frame,'time_ns':time,'rotation_determinant':float(np.linalg.det(R)),'ELOC_anchor_fit_RMSD_A':float(np.sqrt(np.mean(np.sum((fitted-mobile)**2,axis=1)))),'direct_min_modeled_E2_to_any_lysine_NZ_A':float(distance_array(e2,lys.positions).min())}
  for name,a,b in [('ELOC_VHL',heavy(eloc),heavy(vhl)),('VHL_POI',heavy(vhl),heavy(target)),('ELOC_POI',heavy(eloc),heavy(target))]:
   direct=float(distance_array(a.positions,b.positions).min());periodic=float(distance_array(a.positions,b.positions,box=ts.dimensions).min());row[name+'_direct_min_A']=direct;row[name+'_periodic_min_A']=periodic;row[name+'_direct_minus_periodic_A']=direct-periodic
  rows.append(row)
if any(sha256(Path(p))!=digest for p,digest in inputs.items()) or sha256(Path(__file__))!=source_hash:raise ValidationError('Probe input/source changed')
out.mkdir(parents=True);atomic_json(out/'rows.json',rows)
gaps={name:max(r[name+'_direct_minus_periodic_A'] for r in rows) for name in ['ELOC_VHL','VHL_POI','ELOC_POI']}
report={'status':'BOUNDED_DYNAMIC_IMAGING_PROBE_NOT_EXHAUSTIVE_VALIDATION','cases':[{'poi':p,'pose_id':s} for p,s in plan],'frames_per_case':frames,'sampled_frames':len(rows),'maximum_direct_minus_periodic_distance_A':gaps,'max_ELOC_anchor_fit_RMSD_A':max(r['ELOC_anchor_fit_RMSD_A'] for r in rows),'minimum_rotation_determinant':min(r['rotation_determinant'] for r in rows),'maximum_rotation_determinant':max(r['rotation_determinant'] for r in rows),'source_input_sha256':inputs,'audit_source_sha256':source_hash,'rows_sha256':sha256(out/'rows.json'),'original_files_modified':False,'new_MD_or_GNN_training_executed':False,'limitations':['Only three preselected cases at five times each; rare imaging events in other frames/cases remain untested.','Direct/periodic protein-distance differences are diagnostics, not proof that a virtual E2 minimum-image distance is the correct biological replacement.','No reachability formula or original descriptor was changed.']};atomic_json(out/'summary.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ['source_input_sha256','limitations']},indent=2))
