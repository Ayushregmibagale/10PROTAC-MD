"""Independent CPU replay of both completed MD-only arms and final comparison."""
from pathlib import Path
import sys,json,datetime
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'))
from common import *
import pandas as pd

def review(role):
 verify_frozen();p,graphs=load_inputs(role,'cpu');torch.set_num_threads(1)
 run=ROOT/'analysis'/role;c=json.loads((run/'completion.json').read_text());settings=c['settings']
 require(c['status']=='CORRECTED_MD_ONLY_LOPO_COMPLETE' and c['cases']==10 and c['role']==role,'Model arm incomplete')
 require(settings['labels']==p['labels'] and settings['cohort']==p['cohort'] and settings['protocol_sha256']==sha(ROOT/'config/protocol.json') and settings['execution_manifest_sha256']==verify_frozen(),'Completed model source/label mismatch')
 for name,digest in c['outputs_sha256'].items():require(Path(name).name==name and sha(run/name)==digest,'Model artifact drift: '+name)
 frame=pd.read_csv(run/'predictions.tsv',sep='\t');extra=pd.read_csv(run/'prediction_only_fold_scores.tsv',sep='\t')
 require(len(frame)==10 and not frame.poi.duplicated().any() and set(frame.poi)==set(p['cohort']),'Primary fold coverage mismatch')
 require(len(extra)==0 and set(extra.poi)==set(p['prediction_only']) and not extra.duplicated(['held_out','poi']).any(),'Prediction-only coverage mismatch')
 rows=[];replays=[];maximum=0.
 for held in p['cohort']:
  receipt=json.loads((run/f'fold_{held}.json').read_text());checkpoint=run/f'fold_{held}.pt';state=torch.load(checkpoint,map_location='cpu',weights_only=True)
  training=[x for x in p['cohort'] if x!=held];row=frame.loc[frame.poi==held].iloc[0]
  require(state['settings']==settings and state['held_out']==held and state['training_targets']==training==receipt['training_targets'] and state['seed']==fold_seed(p['seed'],held),'Held-out leakage, seed or settings mismatch')
  expected_hyp=[x for x in training if x in p['hypothetical_targets']]
  require(state['training_hypothetical_targets']==expected_hyp==receipt['training_hypothetical_targets'],'Assumed-label training membership mismatch')
  require(bool(row.label_is_hypothetical)==(held in p['hypothetical_targets']),'Hypothetical-label flag lost')
  require(pd.isna(row.experimental_label) if held in p['hypothetical_targets'] else row.experimental_label==p['labels'][held],'Experimental label provenance misrepresented')
  require(receipt['training_class_counts']=={'positive':sum(p['labels'][x] for x in training),'negative':sum(1-p['labels'][x] for x in training)},'Training class counts mismatch')
  require(sha(checkpoint)==receipt['checkpoint_sha256'] and receipt['prediction']==state['prediction'],'Saved weight/prediction mismatch')
  require(row.label==p['labels'][held] and row.has_MD and row.pred==int(row.score>=.5) and abs(row.score-receipt['prediction']['score'])<=1e-12,'Prediction table label/score mismatch')
  require(abs(receipt['replayed_score']-row.score)<=1e-7,'Original saved GPU replay did not pass')
  model=MDEncoder();model.load_state_dict(state['model_state_dict']);score=predict(model,graphs[held]);difference=abs(score-row.score);maximum=max(maximum,difference);require(difference<=1e-5,'CPU held-out replay difference exceeds tolerance')
  for poi in p['prediction_only']:
   erow=extra.loc[(extra.held_out==held)&(extra.poi==poi)].iloc[0];diff=abs(predict(model,graphs[poi])-erow.score);maximum=max(maximum,diff);require(diff<=1e-5 and pd.isna(erow.experimental_label),'CPU prediction-only replay/label mismatch')
  rows.append({'poi':held,'label':int(row.label),'score':float(row.score)});replays.append({'poi':held,'GPU_score':float(row.score),'CPU_score':score,'absolute_difference':difference,'checkpoint_sha256':sha(checkpoint),'training_targets':training})
 require(metrics(rows)==c['metrics'],'Independent metrics mismatch')
 result={'status':'INDEPENDENT_CPU_REPLAY_AND_METRICS_PASS','role':role,'folds':10,'prediction_only_inferences':0,'maximum_absolute_CPU_score_difference':maximum,'CPU_tolerance':1e-5,'original_GPU_tolerance':1e-7,'CPU_torch':str(torch.__version__),'completion_sha256':sha(run/'completion.json'),'execution_manifest_sha256':verify_frozen(),'metrics':c['metrics'],'folds_replayed':replays,'hypothetical_targets_explicitly_held_out':p['hypothetical_targets'],'no_hypothetical_label_is_experimental_validation':True,'new_training_in_review':False,'validator_sha256':sha(__file__)}
 atomic_json(ROOT/'validation'/f'{role}_verification.json',result);return frame,c,extra,result

def main():
 outputs={role:review(role) for role in ['elob','vhl']};p=config();out=ROOT/'results';out.mkdir(exist_ok=False)
 paired=outputs['elob'][0].merge(outputs['vhl'][0],on=['poi','label'],suffixes=('_elob','_vhl'),validate='one_to_one');require(len(paired)==10,'Ten-target paired cohort mismatch');paired.to_csv(out/'paired_predictions.tsv',sep='\t',index=False)
 comparison=[];hyp_rows=[];decomposition=[]
 for role,(frame,c,extra,verification) in outputs.items():
  baseline=json.loads((ROOT/'provenance'/('v1_'+role+'_completion.json')).read_text())
  old=pd.read_csv(ROOT/'provenance'/('v1_'+role+'_predictions.tsv'),sep='\t');old_rows=old.to_dict('records');require(metrics(old_rows,expected_count=8)==baseline['metrics'],'V1 reference metrics mismatch')
  confirmed=frame[~frame.poi.isin(p['hypothetical_targets'])];confirmed_metrics=metrics(confirmed.to_dict('records'),expected_count=8)
  for scope,m in [('v1_original_eight',baseline['metrics']),('v2_all_ten_assumed_labels',c['metrics']),('v2_original_eight_subset',confirmed_metrics)]:comparison.append({'role':role,'scope':scope,**{k:v for k,v in m.items() if not isinstance(v,list)}})
  for _,row in frame[frame.poi.isin(p['hypothetical_targets'])].iterrows():
   receipt=json.loads((ROOT/'analysis'/role/f'fold_{row.poi}.json').read_text())
   require(row.poi not in receipt['training_targets'] and len(receipt['training_targets'])==9,'Hypothetical target leaked into its own fold')
   hyp_rows.append({'role':role,'poi':row.poi,'assumed_label':0,'confirmed_experimental_label':None,'held_out_score':float(row.score),'predicted_class':int(row.pred),'predicted_as_hypothesized_non_degrader':bool(row.pred==0),'fully_excluded_from_own_training':True,'training_targets':';'.join(receipt['training_targets']),'experimental_validation':False})
  positives=frame[frame.label==1]
  for scope,negative in [('confirmed_negatives_ABL1_AXL',frame[frame.poi.isin(['ABL1','AXL'])]),('hypothetical_negatives',frame[frame.poi.isin(p['hypothetical_targets'])])]:
   decomposition.append({'role':role,'scope':scope,'positive_negative_pairs':len(positives)*len(negative),'correctly_ordered_pairs':sum(float(a)>float(b) for a in positives.score for b in negative.score),'ties':sum(float(a)==float(b) for a in positives.score for b in negative.score)})
 pd.DataFrame(comparison).to_csv(out/'v1_v2_comparison.tsv',sep='\t',index=False);pd.DataFrame(hyp_rows).to_csv(out/'hypothetical_target_held_out_scores.tsv',sep='\t',index=False);pd.DataFrame(decomposition).to_csv(out/'rank_pair_decomposition.tsv',sep='\t',index=False)
 summary={'status':'V2_EXPLORATORY_TEN_TARGET_LOPO_COMPLETE','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'cohort':p['cohort'],'labels':p['labels'],'cases':10,'positive_count':6,'negative_count':4,'hypothetical_negative_targets':p['hypothetical_targets'],'hypothetical_labels_are_experimental_validation':False,'metrics':{role:item[1]['metrics'] for role,item in outputs.items()},'original_eight_subset_metrics':{role:metrics(item[0][~item[0].poi.isin(p['hypothetical_targets'])].to_dict('records'),expected_count=8) for role,item in outputs.items()},'training_protocol_sha256':sha(ROOT/'config/protocol.json'),'execution_manifest_sha256':verify_frozen(),'model_completion_sha256':{role:sha(ROOT/'analysis'/role/'completion.json') for role in outputs},'independent_verification_sha256':{role:sha(ROOT/'validation'/f'{role}_verification.json') for role in outputs},'outputs_sha256':{path.name:sha(path) for path in out.iterdir() if path.is_file()},'new_docking_or_MD':False,'independent_biological_validation':False}
 atomic_json(out/'completion.json',summary)
 text='# corrected_md_only_v2: exploratory hypothetical-negative cohort\n\nSame v1 MD-only model, original architecture, 300 epochs, base seed 42, target-bound fold seeds, unweighted BCE loss, optimizer and threshold. All ten targets have retained MD graphs. Six degraded targets: DDR2, EPHA2, EPHB2, MAPK14_full, MET and RIPK2_full. Confirmed/user-labelled non-degraders: ABL1 and AXL. SIK3_full and RPS6KB2_full are assumed non-degraders for this version only; their experimental degradation labels remain unknown.\n\n| Arm/scope | AUROC | Balanced accuracy | Sensitivity | Specificity | Accuracy |\n|---|---:|---:|---:|---:|---:|\n'
 for role in ['elob','vhl']:
  for row in [x for x in comparison if x['role']==role]:text+=f"| {role}: {row['scope']} | {row['AUROC']:.4f} | {row['balanced_accuracy']:.1%} | {row['sensitivity']:.1%} | {row['specificity']:.1%} | {row['accuracy']:.1%} |\n"
 text+='\n## Hypothetical targets: genuine held-out predictions\n\n| Arm | Target | Held-out score | Predicted class | Agrees with assumed non-degrader label | Excluded from own training |\n|---|---|---:|---:|---|---|\n'
 for row in hyp_rows:text+=f"| {row['role']} | {row['poi']} | {row['held_out_score']:.8f} | {row['predicted_class']} | {row['predicted_as_hypothesized_non_degrader']} | True |\n"
 text+='\nEach model trained on the other nine proteins. When SIK3 is held out, its own assumed label does not enter its fitting; RPS6KB2 remains in that training fold under its assumed label, and vice versa. Independent verification checks actual saved training membership, labels/provenance, source/input/weight hashes, held-out score replay and metrics. All 20 saved models were replayed independently on CPU. Agreement against a hypothetical label is not biological correctness or experimental validation.\n\nThe overall negative count increases from two to four and the majority-class accuracy reference falls from 75% to 60%. The comparison above reports whether aggregate results change and separately evaluates v2 scores on the original eight targets, so changes driven solely by the two assumed negatives remain visible. Rank-pair decomposition separates the 12 positive-versus-confirmed-negative pairs from the 12 positive-versus-hypothetical-negative pairs. Adding targets changes geometry, cohort composition and optimizer updates per epoch as well as class balance; it does not isolate a causal effect of balance alone.\n\nA higher AUROC would show consistency with this assumed partition in this retrospective test, not confirm the degradation labels of SIK3/RPS6KB2 or prove general geometric discrimination. Failure would show that this particular frozen GNN/cohort failed; it would not prove MD contact graphs incapable of separation under every model or class balance. No new seed search, score inversion, threshold tuning or relabelling was performed after outcomes. No new docking/MD or graph geometry rebuild occurred: four private NPZ copies change only unknown label metadata from -1 to the hypothesis value 0. Original graph files and the v1 framework/results remain unchanged.\n\nFull-precision scores, per-fold weights/receipts, hypothetical-target flags and both independent CPU audits are retained. No class-weight or E2 intervention is used in v2.\n'
 (ROOT/'RESULTS.md').write_text(text);atomic_json(ROOT/'validation/final_review.json',{'status':'V2_TEN_TARGET_HYPOTHESIS_RESULTS_REVIEWED','results_sha256':sha(ROOT/'RESULTS.md'),'completion_sha256':sha(out/'completion.json'),'validator_sha256':sha(__file__),'hypothetical_labels_are_experimental_validation':False,'no_new_training_during_review':True});print(text)
if __name__=='__main__':main()
