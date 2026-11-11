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
 require(len(extra)==8*len(p['prediction_only']) and set(extra.poi)==set(p['prediction_only']) and not extra.duplicated(['held_out','poi']).any(),'Prediction-only coverage mismatch')
 rows=[];replays=[];maximum=0.
 for held in p['cohort']:
  receipt=json.loads((run/f'fold_{held}.json').read_text());checkpoint=run/f'fold_{held}.pt';state=torch.load(checkpoint,map_location='cpu',weights_only=True)
  training=[x for x in p['cohort'] if x!=held];row=frame.loc[frame.poi==held].iloc[0]
  require(state.get('class_weights')=={str(k):v for k,v in balanced_class_weights(p['labels'],training).items()} and receipt.get('class_weights')==state['class_weights'] and state['settings']==settings and state['held_out']==held and state['training_targets']==training==receipt['training_targets'] and state['seed']==fold_seed(p['seed'],held),'Held-out leakage, seed or settings mismatch')
  require(receipt['training_class_counts']=={'positive':sum(p['labels'][x] for x in training),'negative':sum(1-p['labels'][x] for x in training)},'Training class counts mismatch')
  require(sha(checkpoint)==receipt['checkpoint_sha256'] and receipt['prediction']==state['prediction'],'Saved weight/prediction mismatch')
  require(row.label==p['labels'][held] and row.has_MD and row.pred==int(row.score>=.5) and abs(row.score-receipt['prediction']['score'])<=1e-12,'Prediction table label/score mismatch')
  require(abs(receipt['replayed_score']-row.score)<=1e-7,'Original saved GPU replay did not pass')
  require(state['e2_normalization']==receipt['e2_normalization'],'E2 scaler receipt mismatch');apply_e2_features(graphs,p,training,saved_scaler=state['e2_normalization'])
  model=MDEncoder();model.load_state_dict(state['model_state_dict']);score=predict(model,graphs[held]);difference=abs(score-row.score);maximum=max(maximum,difference);require(difference<=1e-5,'CPU held-out replay difference exceeds tolerance')
  for poi in p['prediction_only']:
   erow=extra.loc[(extra.held_out==held)&(extra.poi==poi)].iloc[0];diff=abs(predict(model,graphs[poi])-erow.score);maximum=max(maximum,diff);require(diff<=1e-5 and pd.isna(erow.experimental_label),'CPU prediction-only replay/label mismatch')
  rows.append({'poi':held,'label':int(row.label),'score':float(row.score)});replays.append({'poi':held,'GPU_score':float(row.score),'CPU_score':score,'absolute_difference':difference,'checkpoint_sha256':sha(checkpoint),'training_targets':training})
 require(metrics(rows)==c['metrics'],'Independent metrics mismatch')
 result={'status':'INDEPENDENT_CPU_REPLAY_AND_METRICS_PASS','role':role,'folds':8,'prediction_only_inferences':8*len(p['prediction_only']),'maximum_absolute_CPU_score_difference':maximum,'CPU_tolerance':1e-5,'original_GPU_tolerance':1e-7,'CPU_torch':str(torch.__version__),'completion_sha256':sha(run/'completion.json'),'execution_manifest_sha256':verify_frozen(),'metrics':c['metrics'],'folds_replayed':replays,'new_training_in_review':False,'validator_sha256':sha(__file__)}
 atomic_json(ROOT/'validation'/f'{role}_verification.json',result);return frame,c,extra,result

def main():
 outputs={role:review(role) for role in ['elob','vhl']};p=config();out=ROOT/'results';out.mkdir(exist_ok=False)
 paired=outputs['elob'][0].merge(outputs['vhl'][0],on=['poi','label'],suffixes=('_elob','_vhl'),validate='one_to_one');require(len(paired)==8,'Paired cohort mismatch');paired.to_csv(out/'paired_predictions.tsv',sep='\t',index=False)
 external=[]
 for role,(_,c,df,reviewed) in outputs.items():
  for poi,g in df.groupby('poi'):external.append({'role':role,'poi':poi,'folds':len(g),'mean_uncalibrated_score':float(g.score.mean()),'min_fold_score':float(g.score.min()),'max_fold_score':float(g.score.max()),'experimental_label':None})
 pd.DataFrame(external,columns=['role','poi','folds','mean_uncalibrated_score','min_fold_score','max_fold_score','experimental_label']).to_csv(out/'prediction_only_ensemble.tsv',sep='\t',index=False)
 elob=outputs['elob'][1]['metrics'];vhl=outputs['vhl'][1]['metrics'];summary={'status':'CORRECTED_MD_ONLY_PAIRED_ANALYSIS_COMPLETE','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'cohort':p['cohort'],'labels':p['labels'],'cases':8,'positive_count':6,'negative_count':2,'EPHA2_degraded':True,'all_evaluation_targets_have_MD':True,'E2_features_used':True,'E2_feature_columns':p['E2_columns'],'static_graphs_used':False,'training_protocol_sha256':sha(ROOT/'config/protocol.json'),'execution_manifest_sha256':verify_frozen(),'metrics':{'elob':elob,'vhl':vhl},'model_completion_sha256':{role:sha(ROOT/'analysis'/role/'completion.json') for role in outputs},'independent_verification_sha256':{role:sha(ROOT/'validation'/f'{role}_verification.json') for role in outputs},'outputs_sha256':{path.name:sha(path) for path in out.iterdir() if path.is_file()},'new_docking_or_MD':False,'independent_biological_validation':False}
 atomic_json(out/'completion.json',summary)
 text=f'''# Class-weighted MD GNN with E2 reachability features

Eight-target corrected-label LOPO: six degraded binders (DDR2, EPHA2, EPHB2, MAPK14_full, MET, RIPK2_full); two non-degraded binders (ABL1, AXL). All targets have MD. The original two-layer MD GAT backbone receives two graph-level E2 geometry descriptors at its classifier head. No static encoder or conventional predictive model is used.

| Arm | AUROC | Accuracy | Balanced accuracy | Sensitivity | Specificity |
|---|---:|---:|---:|---:|---:|
| ELOB control | {elob['AUROC']:.4f} | {elob['accuracy']:.1%} | {elob['balanced_accuracy']:.1%} | {elob['sensitivity']:.1%} | {elob['specificity']:.1%} |
| True VHL | {vhl['AUROC']:.4f} | {vhl['accuracy']:.1%} | {vhl['balanced_accuracy']:.1%} | {vhl['sensitivity']:.1%} | {vhl['specificity']:.1%} |

Loss uses balanced per-target weights N_train/(2*N_class), derived from the seven training targets only. E2 means and standard deviations are also fitted only to those seven targets, saved in each fold checkpoint, and applied to the held-out target without refitting. No missing target scalar is zero-filled. The epoch budget (300), seed (42), optimizer, dropout and threshold remain fixed before fitting. No outcome-based tuning was performed. An always-degraded call gives 75% accuracy in this six-versus-two cohort; only two negatives limit the interpretation of specificity and generalization.

The E2 inputs are: (1) equal-pose mean of the time-averaged minimum E2-point/lysine-NZ distance in Angstrom, across all sampled frames; (2) the equal mean of within-pose bridge-minus-nonbridge 80-Angstrom geometric reachability fractions, using only poses containing both states. They use the audited ligand-connected image convention and true-VHL heavy-contact state mask. All 141 core cases/28,341 frames and the target feature table were independently reproduced from existing packet/contact files. Paired-state coverage varies from 4 to 10 poses per target and is retained in provenance. Both graph roles receive the same E2 values, so role comparison changes the MD graph selection alone.

These quantities are modeled geometric proxies, not experimental E2 occupancy, steric accessibility, catalytic probability, or evidence that all lysines are buried. Corrected ABL1 delta80 is approximately -0.03685, not the legacy -0.0466 value. Frames with alternative image sensitivity remain documented; no pose is removed to improve a model score. E2-candidate data are used at the user's explicit request for exploratory fitting, not upgraded to scientific certification.

All 16 saved model folds passed independent CPU replay. Maximum differences: ELOB {outputs['elob'][3]['maximum_absolute_CPU_score_difference']:.3g}; VHL {outputs['vhl'][3]['maximum_absolute_CPU_score_difference']:.3g}. Source/input hashes, class weights, exact held-out membership, training-only E2 scalers and metrics were independently checked. No new docking or MD was run. The unweighted and weighted MD-only versions remain unchanged. Adding E2 expands the head by two inputs, so any difference includes both descriptor information and the small added head capacity.

No E2-augmented SIK3/RPS6KB2 predictions are produced: audited E2 descriptors for those extension targets are absent, and no fabricated zeros are used. This retrospective, eight-protein LOPO evaluation does not establish independent biological validation or calibrated degradation probabilities. Original structures, trajectories, assay tables and earlier results are preserved.
'''
 (ROOT/'RESULTS.md').write_text(text);atomic_json(ROOT/'validation/final_review.json',{'status':'CORRECTED_MD_ONLY_RESULTS_REVIEWED','results_sha256':sha(ROOT/'RESULTS.md'),'completion_sha256':sha(out/'completion.json'),'validator_sha256':sha(__file__),'no_new_training_during_review':True});print(text)
if __name__=='__main__':main()
