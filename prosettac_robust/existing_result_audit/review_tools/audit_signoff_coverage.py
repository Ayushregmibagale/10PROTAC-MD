"""Read-only 19-target coverage and retained contact-window sensitivity audit.

Produces descriptive sensitivity ranges, never confidence intervals or labels.
"""
from pathlib import Path
import csv,json,sys,collections
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import sha256,atomic_json,ValidationError
project=Path('/home/ayushregmibagale/protac-md')
out=ROOT/'existing_result_audit/scientific_signoff_v1'
if out.exists():raise ValidationError('Fresh signoff audit output required')
inputs={}
def read(path):
 inputs[str(path)]=sha256(path)
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
pred=read(project/'analysis/gnn/e2gnn_lopo_results.tsv')
future=read(project/'analysis/gnn/kinome_predictions.tsv')
manifest=read(project/'random_screen/manifest.tsv')
md=read(project/'analysis/gnn/md_contact_per_pose.tsv')
cohort=json.loads((ROOT/'workflow/config/sjf8240_publication_cohort.json').read_text())
assay_path=ROOT/'sjf8240_review/primary_assay_label_audit.json';inputs[str(assay_path)]=sha256(assay_path)
assay={r['poi']:r for r in json.loads(assay_path.read_text())['rows']}
if {r['poi'] for r in pred}!=set(cohort['historical_labels']) or len(pred)!=17:raise ValidationError('Historical cohort differs')
keys=set(cohort['historical_labels'])|{'SIK3','RPS6KB2'}
if len(keys)!=19:raise ValidationError('Expected 19 confirmed targets')
rows=[];sensitivity=[]
for poi in sorted(keys):
 poses=[r for r in md if r['poi']==poi];stages=collections.Counter(r['stage'] for r in manifest if r['poi']==poi)
 pp=[r for r in pred if r['poi']==poi];fp=[r for r in future if r['poi']==poi+'_full']
 row={'poi':poi,'role':'historical_evaluation' if pp else 'prediction_only','historical_label':pp[0]['label'] if pp else '', 'saved_prediction_present':bool(pp or fp),'manifest_stages':dict(stages),'MD_graph_pose_count':len(poses),'MD_graph_sampled_frame_counts':sorted(set(int(r['n_frames']) for r in poses)),'assay_evidence_status':assay.get(poi,{}).get('status','NO_REVIEWED_EXPERIMENTAL_LABEL'),'true_VHL_contact_reanalysis_cases':0,'unconditional_scientific_signoff':'WITHHELD'}
 for p in sorted((ROOT/'existing_result_audit/contact_core_v1/poses'/poi).glob('*/contacts_A.tsv')):
  inputs[str(p)]=sha256(p);a=np.genfromtxt(p,names=True,delimiter='\t');t=a['time_ns']
  if len(t)!=201 or not np.allclose(t,np.arange(201)/10,atol=5.1e-6,rtol=0):raise ValidationError('Common frame window changed')
  distances=np.maximum(a['true_VHL_heavy_A'],a['POI_heavy_A'])
  if not np.isfinite(distances).all() or (distances<0).any():raise ValidationError('Invalid contact distance')
  windows=[('all_0_20',np.ones(len(t),dtype=bool))]+[(f'({lo},{lo+5}]',(t>lo)&(t<=lo+5)) for lo in [0,5,10,15]]
  for cutoff in [4.,5.,6.]:
   for name,mask in windows:
    if name!='all_0_20' and mask.sum()!=50:raise ValidationError('Expected equal 50-frame blocks')
    sensitivity.append({'poi':poi,'pose_id':p.parent.name,'cutoff_A':cutoff,'window_ns':name,'frames':int(mask.sum()),'occupancy':float(np.mean(distances[mask]<cutoff))})
  row['true_VHL_contact_reanalysis_cases']+=1
 rows.append(row)
summary=[]
for poi in sorted({r['poi'] for r in sensitivity}):
 selected=[r for r in sensitivity if r['poi']==poi]
 means={f'{cutoff:g}A/{window}':float(np.mean([r['occupancy'] for r in selected if r['cutoff_A']==cutoff and r['window_ns']==window])) for cutoff in [4.,5.,6.] for window in ['all_0_20','(0,5]','(5,10]','(10,15]','(15,20]']}
 v=np.array([r['occupancy'] for r in selected if r['cutoff_A']==5. and r['window_ns']=='all_0_20'])
 blocks=[means[f'5A/{window}'] for window in ['(0,5]','(5,10]','(10,15]','(15,20]']]
 leave=(v.sum()-v)/(len(v)-1)
 summary.append({'poi':poi,'poses':len(v),'equal_pose_mean_occupancies':means,'5A_block_mean_range':max(blocks)-min(blocks),'5A_pose_occupancy_min':float(v.min()),'5A_pose_occupancy_max':float(v.max()),'5A_leave_one_pose_out_mean_range':[float(leave.min()),float(leave.max())]})
if any(sha256(Path(p))!=h for p,h in inputs.items()):raise ValidationError('Source drift')
out.mkdir(parents=True)
atomic_json(out/'coverage_19_targets.json',rows);atomic_json(out/'contact_window_sensitivity.json',summary)
with (out/'contact_window_per_pose.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(sensitivity[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(sensitivity)
report={'status':'SCIENTIFIC_SIGNOFF_WITHHELD','confirmed_scope_targets':19,'historical_evaluation_targets':17,'prediction_only_targets':2,'contact_reanalysed_targets':len(summary),'contact_reanalysed_poses':sum(r['true_VHL_contact_reanalysis_cases'] for r in rows),'MD_graph_poses':sum(r['MD_graph_pose_count'] for r in rows),'MD_graph_sampled_frame_counts':sorted({n for r in rows for n in r['MD_graph_sampled_frame_counts']}),'input_sha256':inputs,'audit_source_sha256':sha256(Path(__file__)),'output_sha256':{p.name:sha256(p) for p in out.iterdir()},'original_files_modified':False,'new_training_or_simulations_executed':False,'limitations':['Descriptive window/cutoff/leave-one-pose-out sensitivity only; no confidence intervals or independent-replicate claims.','Contact frames and selected poses are correlated; 141 poses are not 141 independent biological targets.','Only eight core targets have this corrected contact-window audit; two prediction targets have chemistry lineage but not this completed role/contact audit.','Historical reference-label table is not an assay-reviewed training release.','No chemical equivalence, convergence, catalytic mechanism, external validation or calibrated probability is inferred.']}
atomic_json(out/'summary.json',report)
print(json.dumps({k:v for k,v in report.items() if 'sha256' not in k},indent=2))
