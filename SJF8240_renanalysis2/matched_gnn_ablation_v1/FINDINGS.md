# Controlled catalytic-feature GNN ablation — 14 November 2026

This is a new exploratory training comparison on the archived historical molecular model. It is not supplier-SJF8240 validation. Both input variants use the same34-feature TemporalEncoder,training budget,folds,seeds,initial weights and thresholds. All old models,labels and analyses remain unchanged.

## Primary result

| Variant | Matched-pair AUROC | Ranking wins/12 | Balanced accuracy at0.5 | Sensitivity | Specificity |
|---|---:|---:|---:|---:|---:|
| native_features | 0.833 | 10 | 0.500 | 1.000 | 0.000 |
| zero_padded_control | 0.750 | 9 | 0.500 | 1.000 | 0.000 |

Augmented-minus-control AUROC=+0.083;2 pair rankings improve,1 worsen and9 are unchanged. Balanced-accuracy change=+0.000;specificity change=+0.000.

The added features improve within-pair ranking in this fixed run. The augmented classifier still calls every held-out non-degrader a degrader at the frozen0.5 threshold. This does not solve the negative-class problem,even if the ranking AUROC is high.

## ABL1 andAXL discrimination

| Variant | Held-out non-degrader | Matched ranking AUROC | Correct non-degrader calls/6 | Mean negative score |
|---|---|---:|---:|---:|
| native_features | ABL1 | 0.833 | 0 | 0.7725 |
| native_features | AXL | 0.833 | 0 | 0.7766 |
| zero_padded_control | ABL1 | 0.667 | 0 | 0.8073 |
| zero_padded_control | AXL | 0.833 | 0 | 0.8350 |

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
