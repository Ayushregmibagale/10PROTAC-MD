"""Full matched-ablation findings;rank and class discrimination kept separate."""
from pathlib import Path
import sys,json,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/sjf-matched-ablation-mpl')
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'));from contracts import *
def main():
 require(json.loads((ROOT/'validation/independent_cpu_verification.json').read_text())['status']=='PASS','Replaymissing');metric=pd.read_csv(ROOT/'results/metrics.tsv',sep='\t');byneg=pd.read_csv(ROOT/'results/per_negative_target.tsv',sep='\t');targets=pd.read_csv(ROOT/'results/per_target_scores.tsv',sep='\t');pairs=pd.read_csv(ROOT/'results/paired_changes.tsv',sep='\t');receipt=json.loads((ROOT/'validation/comparison_receipt.json').read_text());ctrl=metric[metric.variant=='zero_padded_control'].iloc[0];aug=metric[metric.variant=='native_features'].iloc[0];text='''# Controlled catalytic-feature GNN ablation — 14 November 2026

This is a new exploratory training comparison on the archived historical molecular model. It is not supplier-SJF8240 validation. Both input variants use the same34-feature TemporalEncoder,training budget,folds,seeds,initial weights and thresholds. All old models,labels and analyses remain unchanged.

## Primary result

| Variant | Matched-pair AUROC | Ranking wins/12 | Balanced accuracy at0.5 | Sensitivity | Specificity |
|---|---:|---:|---:|---:|---:|
'''
 for r in metric.itertuples():text+=f'| {r.variant} | {r.matched_pair_AUROC:.3f} | {r.ranking_wins_equivalent:g} | {r.balanced_accuracy_at05:.3f} | {r.sensitivity_at05:.3f} | {r.specificity_at05:.3f} |\n'
 text+=f'\nAugmented-minus-control AUROC={receipt["primary_AUROC_delta"]:+.3f};{receipt["improved_pairs"]} pair rankings improve,{receipt["worsened_pairs"]} worsen and{receipt["unchanged_pairs"]} are unchanged. Balanced-accuracy change={receipt["balanced_accuracy_delta"]:+.3f};specificity change={receipt["specificity_delta"]:+.3f}.\n\n'
 if receipt['primary_AUROC_delta']>0:text+='The added features improve within-pair ranking in this fixed run. '
 elif receipt['primary_AUROC_delta']<0:text+='The added features reduce within-pair ranking in this fixed run. '
 else:text+='The added features do not change the aggregate within-pair ranking in this fixed run. '
 if aug.specificity_at05==0:text+='The augmented classifier still calls every held-out non-degrader a degrader at the frozen0.5 threshold. This does not solve the negative-class problem,even if the ranking AUROC is high.\n\n'
 else:text+='Non-degrader classification at the frozen0.5 threshold is shown above;it remains an eight-target retrospective result.\n\n'
 text+='## ABL1 andAXL discrimination\n\n| Variant | Held-out non-degrader | Matched ranking AUROC | Correct non-degrader calls/6 | Mean negative score |\n|---|---|---:|---:|---:|\n'
 for r in byneg.itertuples():text+=f'| {r.variant} | {r.negative_target} | {r.matched_pair_AUROC:.3f} | {r.nondegrader_predictions_at05} | {r.mean_negative_score:.4f} |\n'
 text+='''
The primary AUROC compares the held-out positive and negative evaluated by the same fitted model. It is not a pooled ROC from probabilities produced by different models. Per-target mean/min/max scores are saved for inspection but remain descriptive because their fitted-model calibrations differ. Each positive has two held-out scores;each negative has six. These are eight independent labelled targets,not twelve independent biological pairs or48 independent predictions. No naive binomial/McNemar significance claim is made for the repeated pairs.

## Controlled protocol

- Six user-confirmed degraders and ABL1/AXL as the two non-degraders. SIK3/RPS6KB2 are unknown and excluded.
-12 held-out positive–negative pairs per variant,24 fitted checkpoints. Training in every fold uses the other six proteins:five positives and one negative. Every pose/window for a held-out protein stays out of training.
-300epochs,Adam learning rate0.001,weight decay0.0001,CosineAnnealingLR300,unweighted binary cross-entropy with logits,hidden64,4heads,dropout0.3,seed42 with deterministic SHA256 fold seed. Same target update order,no early stopping,no tuning or best-checkpoint selection.
- Existing TemporalEncoder with firstGATConv input width29→34 only. Both compared variants have the same34-input architecture. The control's newfive columns arezero;the augmented columns are the previously frozen native-reference descriptors. All original29 features and every graph edge remain identical. The earlier29-input historical0.750 run is not the correct ablation control because wider-layer initialization changes.
- Added lysine-node fields:mean NZ–SY8C3 distance/100Å,reachable20Å fraction,exposure fraction,exposed/reachable/scaffold-compatible fraction,and joint bridge/productive fraction. This is feature engineering on existing graphs,not a new full-scaffold graph architecture.
- No cohort-fitted scaling. Fixed physical unit scaling only. No held-out labels influence feature generation,normalization,optimizer,model selection or classification threshold.

## Weights,receipts and independent CPU verification

Every fold saves a checkpoint with training target IDs,held-out pair,seed,initial-state digest and frozen source manifest hash. Immediate new-instance replay and a separate CPU-verification process reproduce all48 held-out scores. That process also reconstructs the initial weights from each seed,checks semantic final-state digests,confirms paired initial weights match,and audits whole-target split exclusion. Every fold uses exactly1800 optimizer updates.

All weights remain local in results/<variant>/*.pt;compact receipts,predictions and code are checkpointed separately. No old weight file is overwritten. Training was CPU-only with process-local guards against incompatible optional PyG binaries;installed packages were not modified.

## Experimental evidence and chemistry

EXPERIMENTAL_EVIDENCE_REVIEW.md and the version-explicit tables review the supplied RPE1 LF/PS evidence without changing labels. In the named LF table,five of the six positive labels decrease significantly,whileEPHA2 is inconclusive. ABL1 is near-zero/inconclusive;AXL increases versusDMSO. Nonsignificance is not proof of no degradation,and contrasts/cell contexts must not be conflated. PS heavy/newly synthesized and light signals remain separate. Source hashes and all saved result versions are retained. This parallel review did not select the model or labels.

The historical molecular connectivity/actual-coordinate stereochemistry differs from supplier-matched SJF8240. The ablation preserves the actual archived model and cannot by itself validate a supplier-SJF8240 degradation predictor. The added catalytic descriptor is itself reference-dependent and uses a trapped reaction-center surrogate with incomplete atoms. Improved neural ranking,if observed,does not remove those limitations or create independent biological replication.

## Files and decision

results/all_fold_predictions.tsv and paired_changes.tsv show every comparison;per_target_scores.tsv and per_negative_target.tsv expose both non-degraders. validation/independent_cpu_verification.json records all24 model replays. config/protocol.json,input_manifest.json andexecution_manifest.json bind the comparison. The experimental review is a separate label-evidence audit.

Keep all outcomes,including unchanged/worsened folds and false-positive non-degraders. Do not tune thresholds,loss weights,seeds or choose a different catalytic endpoint from these eight outcomes and call the result validation. Broader kinome predictions remain exploratory until chemical/input provenance and independent target-level performance are established.
'''
 (ROOT/'FINDINGS.md').write_text(text)
 from docx import Document
 doc=Document()
 for line in text.splitlines():
  if line.startswith('# '):doc.add_heading(line[2:],0)
  elif line.startswith('## '):doc.add_heading(line[3:],1)
  elif line:doc.add_paragraph(line)
 doc.save(ROOT/'FINDINGS.docx')
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,axes=plt.subplots(1,2,figsize=(11,5),layout='constrained');xx=np.arange(len(pairs));axes[0].bar(xx-.18,pairs.control_positive_score-pairs.control_negative_score,width=.36,label='Zero-padded control');axes[0].bar(xx+.18,pairs.augmented_positive_score-pairs.augmented_negative_score,width=.36,label='Native features');axes[0].axhline(0,color='black',lw=1);axes[0].set_xticks(xx,pairs.pair,rotation=90,fontsize=7);axes[0].set_ylabel('Held-out positive minus negative score');axes[0].legend()
 for variant,color in [('zero_padded_control','#437c9c'),('native_features','#bd7851')]:
  g=targets[(targets.variant==variant)&(targets.label==0)];offset=-.18 if variant=='zero_padded_control' else .18;axes[1].bar(np.arange(2)+offset,g.mean_score,width=.36,label=variant,color=color)
 axes[1].set_xticks(np.arange(2),['ABL1','AXL']);axes[1].axhline(.5,color='black',ls='--',label='Fixed0.5threshold');axes[1].set_ylim(0,1);axes[1].set_ylabel('Mean held-out non-degrader score');axes[1].legend(fontsize=8);fig.suptitle('Matched GNN ablation:ranking versus non-degrader classification');fig.savefig(ROOT/'results/ablation_comparison.png',dpi=180);fig.savefig(ROOT/'results/ablation_comparison.pdf');plt.close(fig)
 atomic_json(ROOT/'validation/scientific_status.json',{'training_complete':True,'models':24,'heldout_scores':48,'CPU_replay':'PASS','AUROC_delta':receipt['primary_AUROC_delta'],'specificity_control':float(ctrl.specificity_at05),'specificity_augmented':float(aug.specificity_at05),'experimental_labels_changed':False,'historical_chemistry':True,'supplier_equivalence':False,'independent_biological_validation':False,'production_degrader_classifier_release':False});print('Saved complete report andplot. AUROC delta:',receipt['primary_AUROC_delta'])
if __name__=='__main__':main()
