"""Report completed paired raw-graph reanalysis without promoting it to validation."""
from pathlib import Path
import json,hashlib
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 outputs={}
 protocol=json.loads((ROOT/'config/protocol.json').read_text())
 for role in ['elob','vhl']:
  run=ROOT/'analysis/gnn'/role/'seed_42';completion=json.loads((run/'completion.json').read_text())
  if completion['status']!='NEW_EXPLORATORY_LOPO_COMPLETE':raise ValueError('Missing completed model arm')
  if completion.get('role')!=role or sorted(completion['cohort'])!=sorted(protocol['primary_cohort']) or completion.get('cases')!=16:raise ValueError('Model cohort/role mismatch')
  for name,digest in completion['output_sha256'].items():
   if Path(name).name!=name or sha(run/name)!=digest:raise ValueError('Model output drift: '+name)
  predictions=pd.read_csv(run/'predictions.tsv',sep='\t')
  if len(predictions)!=16 or predictions.poi.duplicated().any() or set(predictions.poi)!=set(protocol['primary_cohort']):raise ValueError('Incomplete/duplicate held-out predictions')
  if any(protocol['historical_labels'][p]!=label for p,label in zip(predictions.poi,predictions.label)):raise ValueError('Prediction labels drift')
  outputs[role]=(predictions,completion)
 a,sa=outputs['elob'];b,sb=outputs['vhl']
 if set(a.poi)!=set(b.poi):raise ValueError('Paired cohort differs')
 paired=a.merge(b,on=['poi','label'],suffixes=('_elob','_vhl'),validate='one_to_one');
 if len(paired)!=16:raise ValueError('Paired label/cohort mismatch')
 paired['score_change_vhl_minus_elob']=paired.score_vhl-paired.score_elob
 out=ROOT/'analysis/comparison';out.mkdir(exist_ok=False);paired.to_csv(out/'paired_predictions.tsv',sep='\t',index=False)
 experimental=pd.read_csv(ROOT/'experimental_data/target_exports/user_experimental_summary_v2/named_saved_tables_SJF_vs_DMSO.tsv',sep='\t')
 # Preserve all three endpoints; make no new binary label from significance.
 experimental.merge(paired,on='poi',how='left').to_csv(out/'experimental_endpoint_comparison.tsv',sep='\t',index=False)
 unknown=[]
 for role in ['elob','vhl']:
  p=ROOT/'analysis/gnn'/role/'seed_42/prediction_only_fold_scores.tsv';df=pd.read_csv(p,sep='\t')
  for poi,g in df.groupby('poi'):unknown.append({'role':role,'poi':poi,'mean_fold_score':float(g.score.mean()),'min_fold_score':float(g.score.min()),'max_fold_score':float(g.score.max()),'folds':len(g),'experimental_negative_inferred':False})
 pd.DataFrame(unknown).to_csv(out/'prediction_only_ensemble.tsv',sep='\t',index=False)
 report={'status':'PAIRED_NEW_REANALYSIS_COMPLETE','models':{role:{'AUROC':data[1]['AUROC'],'accuracy':data[1]['accuracy'],'cases':data[1]['cases']} for role,data in outputs.items()},'cohort':sa['cohort'],'CDK17_identity_exclusion_prespecified':True,'not_historical_run_replay':True,'experimental_labels_not_changed':True,'new_docking_or_MD':False,'independent_validation_certified':False,'protocol_sha256':sha(ROOT/'config/protocol.json'),'source_model_receipts_sha256':{role:sha(ROOT/'analysis/gnn'/role/'seed_42/completion.json') for role in outputs},'report_script_sha256':sha(__file__),'outputs_sha256':{p.name:sha(p) for p in out.iterdir() if p.is_file()}}
 (out/'completion.json').write_text(json.dumps(report,indent=2)+'\n')
 (ROOT/'RESULTS.md').write_text('# New paired raw-graph reanalysis results\n\n'+f"Completed16-target exploratory LOPO. ELOB control AUROC {sa['AUROC']:.4f}, accuracy {sa['accuracy']:.4f}; true-VHL AUROC {sb['AUROC']:.4f}, accuracy {sb['accuracy']:.4f}.\n\n"+'Both arms use the same repaired pipeline, retained structures, historical labels and fixed400-epoch protocol. They are new runs, not recovery of the missing historical0.850 model. CDK17/CDK14 was excluded before scores. RPE1 endpoints and prediction-only fold ensembles are separate exports. No independent validation, calibrated degradation probability or publication sign-off is inferred. Source hashes, per-fold saved weights and replay checks are retained. No new docking or MD.\n')
 print(json.dumps(report),flush=True)
if __name__=='__main__':main()
