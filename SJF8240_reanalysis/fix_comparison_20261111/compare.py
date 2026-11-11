"""Read-only controlled comparison of completed MD GNN fixes; no fitting."""
from pathlib import Path
import json,hashlib
import pandas as pd
ROOT=Path(__file__).resolve().parent
WORK=ROOT.parent
VERSIONS={'unweighted':'corrected_md_only_v1','weighted':'md_weighted_v1','weighted_E2':'md_weighted_e2_v1'}
def resolve(name):
 direct=WORK/name
 return direct if direct.exists() else WORK/('SJF8240_'+name)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 rows=[];models={};labels=None;graphs={};folds={};reference_settings=None
 for method,directory in VERSIONS.items():
  root=resolve(directory);protocol=json.loads((root/'config/protocol.json').read_text());inputs=json.loads((root/'config/input_manifest.json').read_text())['entries'];maps={(r['poi'],r['role']):r['sha256'] for r in inputs if not r['prediction_only']}
  settings={k:protocol[k] for k in ['epochs','seed','hidden','heads','dropout','learning_rate','weight_decay','optimizer','scheduler','decision_threshold']}
  if reference_settings is None:reference_settings=settings
  else:assert settings==reference_settings,'Training schedule/settings changed across methods'
  for entry in inputs:assert sha(root/entry['path'])==entry['sha256'],'Graph file drift'
  if labels is None:labels=protocol['labels'];graphs=maps
  else:assert labels==protocol['labels'] and graphs==maps,'Cohort/labels/graph arrays differ across methods'
  for role in ['elob','vhl']:
   folder=root/'analysis'/role;completion=json.loads((folder/'completion.json').read_text());review=json.loads((root/'validation'/(role+'_verification.json')).read_text());assert sha(folder/'completion.json')==review['completion_sha256']
   assert review['status']=='INDEPENDENT_CPU_REPLAY_AND_METRICS_PASS'
   for name,digest in completion['outputs_sha256'].items():assert sha(folder/name)==digest
   for poi in protocol['cohort']:
    receipt=json.loads((folder/f'fold_{poi}.json').read_text());assert receipt['training_targets']==[p for p in protocol['cohort'] if p!=poi]
    if method!='unweighted':
     expected=folds.get((poi,role))
     if expected is None:folds[(poi,role)]=receipt['class_weights']
     else:assert expected==receipt['class_weights'],'Loss weights differ between Option1/2'
   row={'method':method,'role':role,**completion['metrics'],'model_completion_sha256':sha(folder/'completion.json'),'verification_sha256':sha(root/'validation'/(role+'_verification.json'))};row['confusion_matrix_labels_0_1']=json.dumps(row['confusion_matrix_labels_0_1']);rows.append(row);models[(method,role)]=completion['metrics']
 pd.DataFrame(rows).to_csv(ROOT/'metrics.tsv',sep='\t',index=False)
 records=[]
 for method,directory in VERSIONS.items():
  for role in ['elob','vhl']:
   frame=pd.read_csv(resolve(directory)/'analysis'/role/'predictions.tsv',sep='\t');frame['method']=method;frame['role']=role;records.append(frame)
 pd.concat(records,ignore_index=True).to_csv(ROOT/'per_target_predictions.tsv',sep='\t',index=False)
 table='| Method (true VHL) | AUROC | Accuracy | Balanced accuracy | Sensitivity | Specificity |\n|---|---:|---:|---:|---:|---:|\n'
 for method,title in [('unweighted','Original corrected-label MD-only'),('weighted','Option 1: balanced loss'),('weighted_E2','Option 2: balanced loss + E2')]:
  m=models[(method,'vhl')];table+=f"| {title} | {m['AUROC']:.4f} | {m['accuracy']:.1%} | {m['balanced_accuracy']:.1%} | {m['sensitivity']:.1%} | {m['specificity']:.1%} |\n"
 text='''# Class weighting and E2 GNN comparison

Both requested interventions are implemented and evaluated using the same eight MD targets, user-confirmed six-versus-two labels, raw graph inputs, 300 epochs and seed42. Only ABL1/AXL are non-degraded; EPHA2 is degraded. Every target has MD. Original versions, data, weights and results are preserved.

'''+table+'''
Option 1 uses balanced per-target BCE weights N_train/(2*N_class) from the seven training targets only. With degrader=1, pos_weight=3 would emphasize degraders, so it is not used. Negative/positive weight ratios are6 or2.5 depending on the held-out class. Tests verify equal total class weighting, the gradient direction and exclusion of held-out labels. This balances optimization; it does not supply more negative examples.

Option 2 appends audited mean minimum E2/lysine distance and delta reach80 to the MD GNN embedding before classification. Means/std are fitted on the seven training targets only, recorded per fold and independently checked. The GAT backbone remains unchanged; the classifier first layer adds two inputs (128 additional weights). Both role controls receive identical scalar values. The feature hypothesis was motivated by already-known target behavior, so this retrospective LOPO is not an independent validation of that feature-selection hypothesis.

The corrected E2 scalars reproduce all141 existing core cases and28,341 frames. They use a ligand-connected coherent image and true-VHL heavy-contact bridge mask. ABL1 delta80 is approximately-0.036848 under this definition; the legacy-0.0466 value uses older definitions. Mean minimum distance is unconditional over frames; delta80 is an average of paired within-pose contrasts (4–10paired poses/target). These are modeled geometric reachability proxies, not proof that all18 lysines are buried, a full E2 steric model or catalytic probability. Source candidate/imaging limitations remain explicit.

The six-versus-two cohort is small and imbalanced; the majority-class accuracy reference is75%. Each negative-held-out fold has just one negative in training. Read accuracy alongside AUROC, balanced accuracy and specificity. No threshold inversion, new seed search, label changes or outcome-based tuning was performed. A score improvement on this cohort would not establish independent biological validity.

Both fixes ran with paired ELOB controls; metrics.tsv contains all six model/method rows and confusion matrices. per_target_predictions.tsv retains full precision held-out scores. All completed model arms pass saved-weight replay, source/input hash checks, exact held-out membership and independent CPU metric recomputation. E2 scaler lineage and fold-local class weights are checked explicitly. E2-augmented extension predictions are omitted because audited E2 descriptors for SIK3/RPS6KB2 are unavailable; no fabricated zeros are used. No new docking, MD or conventional predictor was run.

PyTorch loss definition: https://docs.pytorch.org/docs/2.10/generated/torch.nn.BCEWithLogitsLoss.html
'''
 (ROOT/'COMPARISON.md').write_text(text)
 (ROOT/'verification.json').write_text(json.dumps({'status':'THREE_METHOD_COMPARISON_VERIFIED','labels':labels,'same16_evaluation_graph_hashes':graphs and len(graphs)==16,'same_class_weights_for_weighted_and_E2':True,'same_cohort_seed_epoch_schedule':True,'metrics_sha256':sha(ROOT/'metrics.tsv'),'per_target_predictions_sha256':sha(ROOT/'per_target_predictions.tsv'),'report_sha256':sha(ROOT/'COMPARISON.md'),'comparison_script_sha256':sha(__file__),'new_fitting_executed':False},indent=2)+'\n')
 print(table)
if __name__=='__main__':main()
