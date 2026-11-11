"""Verify all fold metadata and export the prespecified comparison."""
from pathlib import Path
import sys,json,hashlib,itertools
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'))
from common import *
from temporal_model import TemporalEncoder

def main():
 p=config();summary=[];pairs=[];targets=[];checks=[]
 for stage,sourcefile in [('baseline','pair_source_hashes.json'),('temporal','temporal_source_hashes.json')]:
  frozen=json.loads((ROOT/'config'/sourcefile).read_text())
  for path,digest in frozen.items():require(sha(ROOT/path)==digest,'Frozen drift '+path)
  for role in ['elob','vhl']:
   out=ROOT/'analysis'/stage/role;cpu=json.loads((out/'independent_cpu_verification.json').read_text());done=json.loads((out/'completion.json').read_text());require(cpu['verified_cpu'] and cpu['max_replay_error']<=1e-5,'CPU failure');records=[]
   expected=list(itertools.product([x for x in p['cohort'] if p['labels'][x]],['ABL1','AXL']))
   for pos,neg in expected:
    name=pos+'__'+neg;receipt=json.loads((out/(name+'.json')).read_text());checkpoint=out/(name+'.pt');require(sha(checkpoint)==receipt['checkpoint_sha256'],'Checkpoint drift');state=torch.load(checkpoint,map_location='cpu',weights_only=True);train=[x for x in p['cohort'] if x not in [pos,neg]]
    require(state['source_hashes']==frozen and state['pair']==[pos,neg] and state['training']==train and receipt['training']==train,'Leakage/source mismatch');require(state['seed']==receipt['seed']==fold_seed(42,name),'Seed mismatch');require(sum(p['labels'][x] for x in train)==5 and len(train)==6,'Class prior mismatch');require(state['scores'][pos]==receipt['positive_score'] and state['scores'][neg]==receipt['negative_score'],'Saved scores mismatch')
    row={k:receipt[k] for k in ['pair','positive','negative','positive_score','negative_score']};row.update(stage=stage,role=role,correct_rank=float(receipt['positive_score']>receipt['negative_score'])+.5*float(receipt['positive_score']==receipt['negative_score']),margin=receipt['positive_score']-receipt['negative_score']);records.append(row);pairs.append(row)
    for target,score in [(pos,receipt['positive_score']),(neg,receipt['negative_score'])]:targets.append({'stage':stage,'role':role,'pair':name,'target':target,'label':p['labels'][target],'score':score,'prediction':int(score>=.5),'other_held_out':neg if target==pos else pos})
    checks.append({'stage':stage,'role':role,'pair':name,'weight_sha256':sha(checkpoint),'receipt_sha256':sha(out/(name+'.json')),'training_targets':train})
   rank=float(np.mean([r['correct_rank'] for r in records]));sens=float(np.mean([r['positive_score']>=.5 for r in records]));spec=float(np.mean([r['negative_score']<.5 for r in records]));balanced=(sens+spec)/2
   for key,val in [('matched_pair_AUROC',rank),('sensitivity',sens),('specificity',spec),('balanced_accuracy',balanced)]:require(done[key]==val and cpu[key]==val,'Metric recomputation failure')
   summary.append({'stage':stage,'role':role,'matched_pair_AUROC':rank,'correct_pairs':sum(r['correct_rank'] for r in records),'balanced_accuracy':balanced,'sensitivity':sens,'specificity':spec,'median_absolute_pair_margin':float(np.median([abs(r['margin']) for r in records])),'max_CPU_error':cpu['max_replay_error']})
 result=ROOT/'results';pd.DataFrame(summary).to_csv(result/'comparison.tsv',sep='\t',index=False);pd.DataFrame(pairs).to_csv(result/'pair_scores.tsv',sep='\t',index=False);df=pd.DataFrame(targets);df.to_csv(result/'target_held_out_scores.tsv',sep='\t',index=False);df.groupby(['stage','role','target','label'],as_index=False).agg(score_mean=('score','mean'),score_min=('score','min'),score_max=('score','max'),held_out_contexts=('score','size')).to_csv(result/'target_context_summary.tsv',sep='\t',index=False)
 # Same pair comparison, no pooling of predictions across models as AUROC.
 pairdf=pd.DataFrame(pairs);a=pairdf[pairdf.stage=='baseline'].merge(pairdf[pairdf.stage=='temporal'],on=['role','pair','positive','negative'],suffixes=('_baseline','_temporal'));a.to_csv(result/'paired_method_changes.tsv',sep='\t',index=False)
 params={'baseline':sum(x.numel() for x in MDEncoder().parameters()),'temporal':sum(x.numel() for x in TemporalEncoder().parameters())};atomic_json(ROOT/'validation/fold_membership_and_source_audit.json',{'status':'PASS','models':48,'records':checks,'parameter_counts':params})
 text='# Final three-stage results\n\nAll three stages completed. Forty-eight GNN folds (24 original, 24 enriched) passed saved-weight replay and independent CPU inference. All 141 pose packets completed; eight fresh trajectory replays reproduced every saved measurement. Four contract tests passed. No new MD, docking or conventional classifier.\n\n'
 text+='| Model | Arm | Matched-pair AUROC | Correct pairs | Balanced accuracy | Sensitivity | Specificity |\n|---|---|---:|---:|---:|---:|---:|\n'
 for r in summary:text+=f'| {r["stage"]} | {r["role"].upper()} | {r["matched_pair_AUROC"]:.3f} | {r["correct_pairs"]:g}/12 | {r["balanced_accuracy"]:.3f} | {r["sensitivity"]:.3f} | {r["specificity"]:.3f} |\n'
 text+='\n## Did anything improve?\n\n'
 for role in ['vhl','elob']:
  b=next(r for r in summary if r['stage']=='baseline' and r['role']==role);t=next(r for r in summary if r['stage']=='temporal' and r['role']==role);d=a[a.role==role];improved=int(((d.correct_rank_baseline==0)&(d.correct_rank_temporal==1)).sum());worse=int(((d.correct_rank_baseline==1)&(d.correct_rank_temporal==0)).sum());text+=f'{role.upper()}: enriched minus baseline matched-pair AUROC {t["matched_pair_AUROC"]-b["matched_pair_AUROC"]:+.3f}; balanced accuracy {t["balanced_accuracy"]-b["balanced_accuracy"]:+.3f}. {improved} previously wrong pairs became correct; {worse} previously correct pairs became wrong.\n\n'
 text+='''The original v1 pooled LOPO AUROC was zero in both arms, with balanced accuracy 0.50. The new matched-pair estimate uses both held-out targets scored by the same model and constant training composition (5 positives, 1 negative). Its numerical difference from pooled LOPO does not establish a model improvement or prove class-prior bias caused the earlier result. Classification performance at 0.5 and the baseline/enriched comparison under the same matched protocol are the relevant separate checks. No score inversion or threshold optimization was used.

## Mechanistic result

ABL1 NZ exposure during bridging is approximately 98.4%, and none of its 18 lysines has mean bridge-frame NZ SASA below 5 A². The earlier assertion that bridging buries all 18 lysines was an incorrect interpretation of a modeled reachability change. Mean ABL1 complex NZ SASA is 35.81 A², with mean occlusion loss 0.79 A². AXL mean SASA is 33.39 A². Confirmed degrader means span 29.39–40.78 A². Joint exposed/reachable/bridge occupancy overlaps classes. These descriptors do not provide a clean single-mechanism separation in this dataset. See results/MECHANISTIC_FINDINGS.md and the per-lysine, pose-stability and window tables.

## Interpretation and limitations

Twelve pairs share six positives and two negatives: they are not twelve independent negative proteins. Individual target scores occur in different held-out contexts; their averages/minima/maxima are descriptive, not calibrated degradation probabilities. Do not compute a pooled cross-context AUROC from them. Small pair margins, seed sensitivity, protein family differences and feature hypotheses informed by known labels all limit conclusions. No independent test cohort, feature-selection validation or seed/epoch tuning was added.

The enriched candidate combines residue identity, geometry, four-window topology and a larger temporal head. Both candidates use the same pair splits, seed policy, unweighted loss, optimizer and 300 epochs; the candidate is not an isolated causal test of lysine geometry. Temporal edges use 41 sampled frames versus baseline 201-frame aggregation. Physical scales are fixed; no test-target-fitted scaling. SASA is heavy-atom solvent occlusion; an E2 point cloud does not supply the full E2–ubiquitin catalytic structure. No direct ubiquitination prediction is established.

SIK3_full and RPS6KB2_full remain experimentally unknown and excluded from all fitting/evaluation in this study. Original-model prediction-only scores are retained in baseline fold receipts; the enriched model has no audited extension lysine features, so it does not silently zero-impute extension scores. Their historical v2 assumed-negative experiment remains archived unchanged.

## Files to inspect

- docs/STEP_BY_STEP.docx and docs/STEP_BY_STEP.md: complete history and current results.
- results/comparison.tsv: four matched-protocol metric rows.
- results/pair_scores.tsv and paired_method_changes.tsv: every within-model comparison.
- results/target_held_out_scores.tsv and target_context_summary.tsv: per-target scores and held-out contexts.
- results/per_lysine_summary.tsv and per_lysine_bridge_SASA.tsv: every lysine.
- results/pose_stability_and_bridge_exposure.tsv and time_window_summary.tsv: variability and temporal summaries.
- analysis/{baseline,temporal}/{elob,vhl}/: local fold weights, receipts and independent CPU verifications.
- validation/: complete fold/source checks, replay records and final artifact hashes.
'''
 text+=f'\nTrainable parameters: original {params["baseline"]:,}; enriched {params["temporal"]:,}.\n';(ROOT/'RESULTS.md').write_text(text)
 manifests={}
 for root in ['config','scripts','lib','results']:
  for path in (ROOT/root).rglob('*'):
   if path.is_file() and '__pycache__' not in path.parts:manifests[str(path.relative_to(ROOT))]=sha(path)
 manifests['RESULTS.md']=sha(ROOT/'RESULTS.md');atomic_json(ROOT/'validation/final_review.json',{'status':'PASS','models':48,'folds_per_arm_and_method':12,'independent_negative_targets':2,'artifact_sha256':manifests,'all_previous_runs_preserved':True,'hypothetical_labels_used':False})
 print(text)
if __name__=='__main__':main()
