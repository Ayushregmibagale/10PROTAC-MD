"""Audit saved GNN reporting without training or changing original predictions."""
from pathlib import Path
import json,sys,shutil
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score,roc_auc_score
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import sha256,ValidationError
from gnn_contract import subgroup_aurocs
OLD=Path('/home/ayushregmibagale/protac-md')
out=ROOT/'existing_result_audit/gnn_reporting_v1'
if out.exists():raise ValidationError('Immutable output already exists')
out.mkdir();(out/'sources').mkdir()
paths=[OLD/'analysis/gnn'/name for name in ['e2gnn_lopo_results.tsv','e2gnn_lopo_summary.tsv','dual_lopo_results.tsv','dual_lopo_summary.tsv']]
paths += [OLD/'docs'/name for name in ['gnn_robustness_requirements.md','paper_ml_results.md']]
hashes={str(p):sha256(p) for p in paths}
for p in paths:shutil.copy2(p,out/'sources'/p.name)
e2=pd.read_csv(paths[0],sep='\t');dual=pd.read_csv(paths[2],sep='\t')
if e2.poi.duplicated().any() or not set(e2.label)=={0,1}:raise ValidationError('Unexpected retained E2 cohort')
if not np.array_equal(e2.pred.to_numpy(),(e2.prob_degrader>=.5).astype(int).to_numpy()):raise ValidationError('Saved prediction/score mismatch')
if not np.array_equal(e2.correct.to_numpy(),(e2.pred==e2.label).astype(int).to_numpy()):raise ValidationError('Saved correctness mismatch')
groups=subgroup_aurocs(e2.label,e2.prob_degrader,{'all':np.ones(len(e2),dtype=bool),'md':(e2.has_md==1).to_numpy(),'nonmd':(e2.has_md==0).to_numpy()})
for name,mask in [('all',np.ones(len(e2),dtype=bool)),('md',(e2.has_md==1).to_numpy()),('nonmd',(e2.has_md==0).to_numpy())]:
 d=e2.loc[mask];groups[name]['correct']=int(d.correct.sum());groups[name]['accuracy']=float(d.correct.mean()) if len(d) else None
 groups[name]['targets']=list(d.poi)
current_dual=dual[dual.poi.isin(e2.poi)].copy()
if current_dual.poi.duplicated().any() or set(current_dual.poi)!=set(e2.poi):raise ValidationError('Dual cohort mismatch')
current_dual=current_dual.set_index('poi').loc[e2.poi].reset_index()
if not np.array_equal(current_dual.label,e2.label):raise ValidationError('Historical table labels disagree')
metadata=[]
for _,row in e2.iterrows():
 metadata.append({'poi':row.poi,'label':int(row.label),'has_md':int(row.has_md),'has_e2_saved_flag':int(row.has_e2),'saved_mean_E2_is_finite':bool(pd.notna(row.e2_mean_min_dist_A)),'interpretation':'has_e2 is not a per-field finite-value certificate'})
pd.DataFrame(metadata).to_csv(out/'availability_semantics.tsv',sep='\t',index=False)
report={'status':'RETAINED_SCORE_REPORTING_RECONCILED_NOT_NEW_VALIDATION','E2_saved_score_metrics':groups,
 'current_dual_table':{'rows':len(dual),'reviewed_cohort_rows':len(current_dual),'unknown_label_rows':int((dual.label==-1).sum()),'reviewed_cohort_unique_scores':sorted(map(float,current_dual.prob_degrader.unique())),'reviewed_cohort_has_md_counts':{str(k):int(v) for k,v in current_dual.has_md.value_counts().items()},'reviewed_cohort_AUROC_from_current_scores':float(roc_auc_score(current_dual.label,current_dual.prob_degrader)),'reviewed_cohort_accuracy_from_current_predictions':float(accuracy_score(current_dual.label,current_dual.pred))},
 'not_supported':['A no-MD AUROC of 1.000: this subgroup has no positive labels, so AUROC is undefined.','A paired reconstruction of the older 04d AUROC=0.850 versus 04g AUROC=0.800 using the current dual TSV: its scores and cohort have changed.','A generalization guarantee from nine correct negative calls or a fixed sample-size threshold.','Claiming that adding E2 improves biology or is causal solely from the selected 17-target score table.','Calling a finite docking search failure an experimental negative, or attributing assay disagreement to an experimental false negative from model output.'],
 'source_sha256':hashes,'audit_sha256':sha256(Path(__file__)),'new_training_or_MD':False,'original_outputs_modified':False,'limits':['Scores are retained rounded values, not reconstructed unrounded fold predictions.','Labels are historical and unresolved assay/target correspondence is unchanged.','Availability confounding was audited separately; this report fits no metadata model.','No author-selected assay hierarchy, experimental relabeling or performance-improving case deletion.']}
if any(sha256(Path(p))!=v for p,v in hashes.items()):raise ValidationError('Original file drift')
(out/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print(json.dumps({'groups':groups,'current_dual_table':report['current_dual_table']},indent=2))
