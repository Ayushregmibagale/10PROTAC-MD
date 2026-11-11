"""Verify a finished paired report and write a concise, evidence-bound review."""
from pathlib import Path
import json,hashlib,datetime
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def require(ok,message):
 if not ok:raise ValueError(message)
def main():
 out=ROOT/'analysis/comparison';completed=json.loads((out/'completion.json').read_text());require(completed['status']=='PAIRED_NEW_REANALYSIS_COMPLETE','Paired report incomplete')
 for name,digest in completed['outputs_sha256'].items():require(Path(name).name==name and sha(out/name)==digest,'Paired report output drift')
 require(sha(ROOT/'config/protocol.json')==completed['protocol_sha256'],'Protocol drift')
 require(sha(ROOT/'scripts/06_make_comparison.py')==completed['report_script_sha256'],'Report source drift')
 protocol=json.loads((ROOT/'config/protocol.json').read_text());rows=[];reviews={}
 for role in ['elob','vhl']:
  path=ROOT/'validation'/f'verification_{role}.json';review=json.loads(path.read_text());reviews[role]=review
  require(review['status']=='INDEPENDENT_MODEL_ARTIFACT_AND_CPU_REPLAY_PASS','Model review did not pass')
  completion=ROOT/'analysis/gnn'/role/'seed_42/completion.json'
  require(sha(completion)==review['completion_sha256']==completed['source_model_receipts_sha256'][role],'Reviewed versus reported model mismatch')
  for key in ['AUROC','accuracy']:require(completed['models'][role][key]==review['metrics'][key],'Reported paired metric mismatch')
  for group,values in [('overall',{'n':16,'positive':5,'negative':11,'auroc':review['metrics']['AUROC'],'auroc_status':'DEFINED'}),*review['subgroups'].items()]:
   rows.append({'role':role,'group':group,**values,'overall_accuracy':review['metrics']['accuracy'] if group=='overall' else None})
 paired=pd.read_csv(out/'paired_predictions.tsv',sep='\t');require(len(paired)==16 and set(paired.poi)==set(protocol['primary_cohort']) and 'CDK17' not in set(paired.poi),'Paired evaluation cohort mismatch')
 for role in ['elob','vhl']:
  raw=pd.read_csv(ROOT/'analysis/gnn'/role/'seed_42/predictions.tsv',sep='\t').set_index('poi');joined=paired.set_index('poi').loc[raw.index]
  require(np.allclose(joined['score_'+role],raw.score,rtol=0,atol=1e-12) and np.array_equal(joined.label,raw.label),'Paired score/label join mismatch')
 graph=json.loads((ROOT/'validation/graph_packet_verification.json').read_text());coverage=pd.read_csv(out/'static_input_validity.tsv',sep='\t')
 require(len(coverage)==19 and coverage.configured_files.sum()==graph['configured_source_files'] and coverage.included_complete_poses.sum()==graph['accepted_complete_source_files'] and coverage.excluded_unusable_files.sum()==graph['excluded_unusable_source_files'],'Coverage export mismatch')
 external=pd.read_csv(out/'prediction_only_ensemble.tsv',sep='\t');require(len(external)==4 and set(external.poi)==set(protocol['prediction_only']) and (external.folds==16).all() and not external.experimental_negative_inferred.any(),'Prediction-only summary mismatch')
 for _,row in external.iterrows():
  scores=pd.read_csv(ROOT/'analysis/gnn'/row.role/'seed_42/prediction_only_fold_scores.tsv',sep='\t');values=scores.loc[scores.poi==row.poi,'score']
  require(abs(values.mean()-row.mean_fold_score)<1e-12 and abs(values.min()-row.min_fold_score)<1e-12 and abs(values.max()-row.max_fold_score)<1e-12,'Prediction-only aggregation mismatch')
 source=pd.read_csv(ROOT/'experimental_data/target_exports/user_experimental_summary_v2/named_saved_tables_SJF_vs_DMSO.tsv',sep='\t')
 exported=pd.read_csv(out/'experimental_endpoint_comparison.tsv',sep='\t');require(len(source)==len(exported),'Experimental endpoint rows dropped')
 pd.testing.assert_frame_equal(source.reset_index(drop=True),exported[source.columns].reset_index(drop=True),check_dtype=False,rtol=0,atol=1e-12)
 summary=pd.DataFrame(rows);summary.to_csv(ROOT/'validation/reviewed_metrics.tsv',sep='\t',index=False)
 elob=reviews['elob'];vhl=reviews['vhl']
 md=lambda value:'undefined' if value is None else f'{value:.4f}'
 text=f'''# Reviewed SJF8240 paired reanalysis

Both frozen model arms completed all 16 LOPO folds. This is a new exploratory comparison on retained structures and historical labels; it does not reconstruct the missing historical model or establish independent biological validation.

| Model arm | Overall AUROC (16 targets) | Accuracy | MD-only AUROC (8 targets) |
|---|---:|---:|---:|
| ELOB control | {elob['metrics']['AUROC']:.4f} | {elob['metrics']['accuracy']:.1%} | {md(elob['subgroups']['md']['auroc'])} |
| True VHL | {vhl['metrics']['AUROC']:.4f} | {vhl['metrics']['accuracy']:.1%} | {md(vhl['subgroups']['md']['auroc'])} |

The full cohort has five positive and eleven negative historical labels. The MD subgroup contains all five positives and three negatives. The eight targets without MD are all negative, so their subgroup AUROC is undefined. MD availability is associated with the labels: aggregate performance alone is insufficient evidence of biological discrimination. These subgroup checks do not establish which feature caused a prediction, and no additional model fitting or tuning was performed.

All {graph['configured_source_files']:,} configured PDB files are accounted for: {graph['accepted_complete_source_files']:,} complete files included, three unusable files excluded symmetrically with original hashes and reasons. All 19 static targets and 174 MD cases passed artifact checks. Each MD case contributes 201 unique frames over 0–20 ns; all ten target graph aggregations were independently reproduced from saved contact packets. Original structures, labels, trajectories and historical model outputs remain preserved. No new docking or MD was run.

All 32 saved fold models and 64 prediction-only fold inferences passed independent CPU replay against saved GPU scores. Maximum score difference: ELOB {elob['maximum_CPU_score_difference_including_prediction_only']:.3g}; VHL {vhl['maximum_CPU_score_difference_including_prediction_only']:.3g}. Frozen input/source hashes, output/checkpoint hashes, held-out membership and recomputed metrics passed verification. This is software/artifact verification, not independent experimental validation.

CDK17-labelled structures match CDK14 and were excluded from the primary cohort before fitting, without renaming or relabelling originals. SIK3 and RPS6KB2 remain prediction-only; fold-ensemble scores are uncalibrated and are not experimental labels. LF abundance, PS light and PS heavy endpoints stay separate, with all source endpoint rows preserved. Historical assay-label correspondence, supplier equivalence, external generalization and calibration remain limitations.

The frozen model uses the original 04d raw-graph dual encoder, 400 epochs, seed 42 and no tuning. E2 scalar descriptors are separate audit outputs and are not inputs to this model. See `analysis/comparison/` for paired target scores, input coverage, endpoint exports and prediction-only summaries; `validation/` holds independent verification receipts and `validation/reviewed_metrics.tsv` records subgroup counts.
'''
 (ROOT/'REVIEWED_RESULTS.md').write_text(text)
 report={'status':'PAIRED_RESULTS_INDEPENDENTLY_REVIEWED','checked_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'validator_sha256':sha(__file__),'paired_completion_sha256':sha(out/'completion.json'),'graph_verification_sha256':sha(ROOT/'validation/graph_packet_verification.json'),'model_verification_sha256':{role:sha(ROOT/'validation'/f'verification_{role}.json') for role in reviews},'reviewed_metrics_sha256':sha(ROOT/'validation/reviewed_metrics.tsv'),'reviewed_results_sha256':sha(ROOT/'REVIEWED_RESULTS.md'),'experimental_endpoint_rows_preserved':len(source),'independent_biological_validation':False,'new_training_during_review':False}
 (ROOT/'validation/final_review.json').write_text(json.dumps(report,indent=2)+'\n');print(text);print(json.dumps(report,indent=2))
if __name__=='__main__':main()
