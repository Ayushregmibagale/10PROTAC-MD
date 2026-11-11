"""Independent CPU replay of both completed MD-only arms and final comparison."""
from pathlib import Path
import sys,json,datetime
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'))
from common import *
import pandas as pd

def review(role):
 verify_frozen();p,graphs=load_inputs(role,'cpu');torch.set_num_threads(1)
 run=ROOT/'analysis'/role;c=json.loads((run/'completion.json').read_text());settings=c['settings']
 require(c['status']=='CORRECTED_MD_ONLY_LOPO_COMPLETE' and c['cases']==8 and c['role']==role,'Model arm incomplete')
 require(settings['labels']==p['labels'] and settings['cohort']==p['cohort'] and settings['protocol_sha256']==sha(ROOT/'config/protocol.json') and settings['execution_manifest_sha256']==verify_frozen(),'Completed model source/label mismatch')
 for name,digest in c['outputs_sha256'].items():require(Path(name).name==name and sha(run/name)==digest,'Model artifact drift: '+name)
 frame=pd.read_csv(run/'predictions.tsv',sep='\t');extra=pd.read_csv(run/'prediction_only_fold_scores.tsv',sep='\t')
 require(len(frame)==8 and not frame.poi.duplicated().any() and set(frame.poi)==set(p['cohort']),'Primary fold coverage mismatch')
 require(len(extra)==16 and set(extra.poi)==set(p['prediction_only']) and not extra.duplicated(['held_out','poi']).any(),'Prediction-only coverage mismatch')
 rows=[];replays=[];maximum=0.
 for held in p['cohort']:
  receipt=json.loads((run/f'fold_{held}.json').read_text());checkpoint=run/f'fold_{held}.pt';state=torch.load(checkpoint,map_location='cpu',weights_only=True)
  training=[x for x in p['cohort'] if x!=held];row=frame.loc[frame.poi==held].iloc[0]
  require(state['settings']==settings and state['held_out']==held and state['training_targets']==training==receipt['training_targets'] and state['seed']==fold_seed(p['seed'],held),'Held-out leakage, seed or settings mismatch')
  require(receipt['training_class_counts']=={'positive':sum(p['labels'][x] for x in training),'negative':sum(1-p['labels'][x] for x in training)},'Training class counts mismatch')
  require(sha(checkpoint)==receipt['checkpoint_sha256'] and receipt['prediction']==state['prediction'],'Saved weight/prediction mismatch')
  require(row.label==p['labels'][held] and row.has_MD and row.pred==int(row.score>=.5) and abs(row.score-receipt['prediction']['score'])<=1e-12,'Prediction table label/score mismatch')
  require(abs(receipt['replayed_score']-row.score)<=1e-7,'Original saved GPU replay did not pass')
  model=MDEncoder();model.load_state_dict(state['model_state_dict']);score=predict(model,graphs[held]);difference=abs(score-row.score);maximum=max(maximum,difference);require(difference<=1e-5,'CPU held-out replay difference exceeds tolerance')
  for poi in p['prediction_only']:
   erow=extra.loc[(extra.held_out==held)&(extra.poi==poi)].iloc[0];diff=abs(predict(model,graphs[poi])-erow.score);maximum=max(maximum,diff);require(diff<=1e-5 and pd.isna(erow.experimental_label),'CPU prediction-only replay/label mismatch')
  rows.append({'poi':held,'label':int(row.label),'score':float(row.score)});replays.append({'poi':held,'GPU_score':float(row.score),'CPU_score':score,'absolute_difference':difference,'checkpoint_sha256':sha(checkpoint),'training_targets':training})
 require(metrics(rows)==c['metrics'],'Independent metrics mismatch')
 result={'status':'INDEPENDENT_CPU_REPLAY_AND_METRICS_PASS','role':role,'folds':8,'prediction_only_inferences':16,'maximum_absolute_CPU_score_difference':maximum,'CPU_tolerance':1e-5,'original_GPU_tolerance':1e-7,'CPU_torch':str(torch.__version__),'completion_sha256':sha(run/'completion.json'),'execution_manifest_sha256':verify_frozen(),'metrics':c['metrics'],'folds_replayed':replays,'new_training_in_review':False,'validator_sha256':sha(__file__)}
 atomic_json(ROOT/'validation'/f'{role}_verification.json',result);return frame,c,extra,result

def main():
 outputs={role:review(role) for role in ['elob','vhl']};p=config();out=ROOT/'results';out.mkdir(exist_ok=False)
 paired=outputs['elob'][0].merge(outputs['vhl'][0],on=['poi','label'],suffixes=('_elob','_vhl'),validate='one_to_one');require(len(paired)==8,'Paired cohort mismatch');paired.to_csv(out/'paired_predictions.tsv',sep='\t',index=False)
 external=[]
 for role,(_,c,df,reviewed) in outputs.items():
  for poi,g in df.groupby('poi'):external.append({'role':role,'poi':poi,'folds':len(g),'mean_uncalibrated_score':float(g.score.mean()),'min_fold_score':float(g.score.min()),'max_fold_score':float(g.score.max()),'experimental_label':None})
 pd.DataFrame(external).to_csv(out/'prediction_only_ensemble.tsv',sep='\t',index=False)
 elob=outputs['elob'][1]['metrics'];vhl=outputs['vhl'][1]['metrics'];summary={'status':'CORRECTED_MD_ONLY_PAIRED_ANALYSIS_COMPLETE','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'cohort':p['cohort'],'labels':p['labels'],'cases':8,'positive_count':6,'negative_count':2,'EPHA2_degraded':True,'all_evaluation_targets_have_MD':True,'static_graphs_used':False,'training_protocol_sha256':sha(ROOT/'config/protocol.json'),'execution_manifest_sha256':verify_frozen(),'metrics':{'elob':elob,'vhl':vhl},'model_completion_sha256':{role:sha(ROOT/'analysis'/role/'completion.json') for role in outputs},'independent_verification_sha256':{role:sha(ROOT/'validation'/f'{role}_verification.json') for role in outputs},'outputs_sha256':{path.name:sha(path) for path in out.iterdir() if path.is_file()},'new_docking_or_MD':False,'independent_biological_validation':False}
 atomic_json(out/'completion.json',summary)
 text=f'''# Corrected-label MD-only GNN results

Completed a new eight-target MD-only LOPO comparison with the user-confirmed EPHA2 label. Six degraded binders: DDR2, EPHA2, EPHB2, MAPK14_full, MET and RIPK2_full. Two non-degraded binders: ABL1 and AXL. Every evaluated target has an MD graph; static docking graphs and E2 scalar features are absent from the model.

| Model arm | AUROC | Accuracy | Balanced accuracy | Sensitivity | Specificity |
|---|---:|---:|---:|---:|---:|
| ELOB control | {elob['AUROC']:.4f} | {elob['accuracy']:.1%} | {elob['balanced_accuracy']:.1%} | {elob['sensitivity']:.1%} | {elob['specificity']:.1%} |
| True VHL (primary) | {vhl['AUROC']:.4f} | {vhl['accuracy']:.1%} | {vhl['balanced_accuracy']:.1%} | {vhl['sensitivity']:.1%} | {vhl['specificity']:.1%} |

The cohort has only two negative targets. An always-degraded call would achieve 75% accuracy; interpret accuracy together with AUROC, balanced accuracy, sensitivity and specificity. These are exploratory internal LOPO results, not independent biological validation or calibrated degradation probabilities. Each negative-held-out fold trains on only one other negative target. No optimization was performed after seeing scores.

The frozen original 04h MD-only model uses a two-layer GAT with mean/max pooling and its original classifier head, hidden size64, four heads, dropout0.3, 300 epochs, Adam lr0.001/weight decay0.0001, cosine schedule and seed42 with target-bound fold seeds. The original model class AST was verified exactly. Training holds out the whole target; saved weights and inference replay are recorded for every fold. Independent CPU replay checks all16 fold models and32 prediction-only fold inferences. Maximum CPU score differences: ELOB {outputs['elob'][3]['maximum_absolute_CPU_score_difference']:.3g}, VHL {outputs['vhl'][3]['maximum_absolute_CPU_score_difference']:.3g}.

This version reuses the existing MD contact graphs, correcting EPHA2 label metadata without changing graph features. No new docking or MD was run. Old16-target dual-encoder results remain frozen under their old labels. This new comparison changes cohort, labels and architecture together and cannot isolate any one change relative to that historical run. Within this version the ELOB/VHL arms use the same model, target labels and training settings.

The retained MD sources cover141 systems across the eight evaluation targets. Prediction-only SIK3 and RPS6KB2 reuse33 additional systems; their fold-ensemble scores are exported separately and are not validation labels. The retained structure is MAPK14_full; no MAPK15 structure was substituted. Historical experimental endpoint correspondence, retained chemistry/simulation heterogeneity, external generalization and calibration remain limitations. LF/PS assay measurements and all prior models, metrics and source data are unchanged.

`results/paired_predictions.tsv` gives all held-out scores; `results/completion.json` includes confusion matrices and average precision. Protocol/source/input hashes, per-fold weights, prediction-only scores and independent review receipts are retained in this version.
'''
 (ROOT/'RESULTS.md').write_text(text);atomic_json(ROOT/'validation/final_review.json',{'status':'CORRECTED_MD_ONLY_RESULTS_REVIEWED','results_sha256':sha(ROOT/'RESULTS.md'),'completion_sha256':sha(out/'completion.json'),'validator_sha256':sha(__file__),'no_new_training_during_review':True});print(text)
if __name__=='__main__':main()
