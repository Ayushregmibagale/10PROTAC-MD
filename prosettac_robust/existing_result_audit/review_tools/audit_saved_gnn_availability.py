"""Descriptive audit of saved predictions and availability; fit no model."""
from pathlib import Path
import json,sys
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score,accuracy_score,balanced_accuracy_score,matthews_corrcoef
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import sha256,atomic_json,ValidationError
source=Path('/home/ayushregmibagale/protac-md/analysis/gnn/e2gnn_lopo_results.tsv');digest=sha256(source);df=pd.read_csv(source,sep='\t')
if df.poi.duplicated().any() or set(df.label)!={0,1}:raise ValidationError('Unique binary-labelled saved cohort required')
for col in ['has_md','has_e2','pred']:
 if not set(df[col])<={0,1}:raise ValidationError('Nonbinary metadata field')
if not np.isfinite(df.prob_degrader).all() or ((df.prob_degrader<0)|(df.prob_degrader>1)).any():raise ValidationError('Invalid saved scores')
y=df.label.to_numpy();rows=[]
for name,score,pred in [('saved_GNN',df.prob_degrader.to_numpy(),df.pred.to_numpy()),('MD_available_flag',df.has_md.to_numpy(),df.has_md.to_numpy()),('E2_available_flag',df.has_e2.to_numpy(),df.has_e2.to_numpy())]:
 rows.append({'descriptor':name,'cases':len(y),'correct':int((pred==y).sum()),'accuracy':float(accuracy_score(y,pred)),'balanced_accuracy':float(balanced_accuracy_score(y,pred)),'AUROC':float(roc_auc_score(y,score)),'MCC':float(matthews_corrcoef(y,pred))})
counts=[]
for col in ['has_md','has_e2']:
 for flag in [0,1]:
  subset=df[df[col]==flag];counts.append({'field':col,'available':flag,'cases':len(subset),'recorded_positive_labels':int(subset.label.sum()),'recorded_negative_labels':int((subset.label==0).sum())})
if sha256(source)!=digest:raise ValidationError('Saved predictions changed during audit')
out=ROOT/'existing_result_audit/gnn_availability_audit_v1'
if (out/'summary.json').exists():raise ValidationError('Audit output already used')
pd.DataFrame(rows).to_csv(out/'descriptive_comparison.tsv',sep='\t',index=False,lineterminator='\n')
report={'status':'DESCRIPTIVE_SAVED_COHORT_AUDIT_NOT_NEW_VALIDATION','cases':len(df),'source_sha256':digest,'source':str(source),'audit_source_sha256':sha256(Path(__file__)),'availability_counts':counts,'comparison':rows,'model_fitted':False,'GNN_retrained':False,'new_MD_or_docking_executed':False,'original_results_modified':False,'threshold_for_binary_flags':'.5, without tuning','limitations':['Uses recorded historical labels; does not resolve experimental label provenance.','These availability checks were chosen as a retrospective audit, not preregistered independent validation.','Saved GNN AUROC is recomputed from its retained rounded scores; original unrounded logits/weights are not inferred.','Metadata flags are descriptive rank scores, not calibrated probabilities.','Small selected cohort and availability/selection confounding prevent a causal or generalization conclusion.']}
atomic_json(out/'summary.json',report);print(json.dumps({'comparison':rows,'availability_counts':counts},indent=2))
