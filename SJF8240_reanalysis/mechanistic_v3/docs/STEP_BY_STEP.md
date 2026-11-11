# SJF8240 analysis: step-by-step record

Last updated: 2026-11-11. New three-stage investigation is COMPLETE. Original source trajectories and earlier results are preserved.

## 1. Data and chain audit
Reused existing simulations: 141 cases for eight labelled targets; 33 extension cases. Each case contains 201 frames spanning 20 ns. Corrected the legacy ELOB/true-VHL naming: chain A is ELOB, chain C is VHL. Audited input damage, chemistry/provenance and periodic images. No new docking or MD.

## 2. Historical paired model
The 16-target dual model gave AUROC 0.727 in both arms, but its MD subgroup AUROC was zero. Cross-group comparisons explained aggregate performance. This used the historical EPHA2 label and is archived (d5c34f1).

## 3. Corrected labels and MD-only baseline
Six degraded targets: DDR2, EPHA2, EPHB2, MAPK14_full, MET, RIPK2_full. Two non-degraded targets: ABL1, AXL. MAPK14 is the actual available system. SIK3_full and RPS6KB2_full remain experimentally unknown.
The corrected_md_only_v1 GNN used 300 epochs, seed 42 and unweighted BCE. Both ELOB/VHL: AUROC 0, balanced accuracy 0.50, sensitivity 1, specificity 0. All 16 weights independently replayed (b5ed245).

## 4. Class-weight correction
A positive weight of 3 would emphasize degraders, contrary to the intended negative emphasis. Instead used training-fold balanced per-target weights. md_weighted_v1: ELOB AUROC 0.333, balanced accuracy 0.333; VHL AUROC 0.083, balanced accuracy 0.25. No useful improvement.

## 5. Two E2 graph features
Added audited mean E2-to-lysine distance and paired bridge-associated reach80 contrast, scaled on training targets only. md_weighted_e2_v1: both AUROC 0; balanced accuracy ELOB 0.083, VHL 0.167. All 32 weighted/E2 weights passed independent CPU replay (c9f17af). Distance/reach are geometric proxies, not evidence that lysines are buried or ubiquitination occurs.

## 6. Hypothetical ten-target cohort
corrected_md_only_v2 imposed label zero on the experimentally unknown SIK3/RPS6KB2, solely as a hypothesis test. Full cohort ELOB AUROC 0.0417, balanced accuracy 0.4167; VHL AUROC 0.1667, balanced accuracy 0.125. SIK3 scores: ELOB 0.52507, VHL 0.55878. RPS6KB2: ELOB 0.55559, VHL 0.52446. All four exceed 0.5. These are not experimental validation. Twenty models passed CPU replay (b6f095f).

## 7. New evaluation: leave-pair-out — completed
Hold out one degraded and one non-degraded target together: 12 pairs per arm. Every training fold has five positives and one negative. Compare held-out scores within the same model; do not pool scores from differently trained models as an AUC. Report fraction of correctly ordered pairs, sensitivity and specificity at fixed 0.5. Pairs share targets and are not independent biological replicates. Same original model, 300 epochs, seed 42 and unweighted loss.

## 8. Per-lysine mechanistic analysis — completed
Use 41 existing frames per pose (0, 0.5, ..., 20 ns), equal pose weighting. Measure lysine NZ solvent-accessible area in the coherent complex and isolated target, E2-point distance and exposed/reachable bridge occupancy. Four fixed time windows retain coarse dynamics. No label-selected thresholds. These are geometric proxies; full E2–ubiquitin catalytic orientation is not available from a point cloud.

## 9. Lysine-aware temporal GNN — completed
Compare with the same 12 held-out pairs, fixed training budget and seed policy. Add amino-acid identity and measured lysine geometry, retaining four time windows. Freeze the feature/model specification before fitting. Unknown extension targets are never used as experimental negatives.

## Reading the files
Canonical analyses: /home/ayushregmibagale/10PROTAC-MD/SJF8240_reanalysis/. Earlier folders remain unchanged. New work: mechanistic_v3; workspace mirror: SJF8240_mechanistic_v3. Per-fold weights, source hashes and replay receipts will accompany completed stages. No claim of general predictive validity is possible from only two confirmed negatives.

## Recovered failures and safeguards

- Three unusable static PDB files were excluded with ledgers; coordinate-overflow parsing was repaired and tested. Original files were preserved.
- GPU preflight initially rejected missing graph-schema metadata; corrected explicit schema checks before production fitting.
- The first corrected-v1 CPU replay inherited a Python 3.10 package overlay under Python 3.13; removing that overlay fixed replay without retraining.
- A historical report-only Path/string hash bug was fixed after archiving its failure; model outputs were unchanged.
- Failed logs, original source snapshots and all earlier model versions remain in their versioned folders. See the corresponding validation/recovery directories and project memory.

## Detailed audit trail

Earlier complete reports: REVIEWED_RESULTS.md; corrected_md_only_v1/RESULTS.md; fix_comparison_20261111/COMPARISON.md; corrected_md_only_v2/RESULTS.md and INTERPRETATION.md, all beneath the canonical SJF8240_reanalysis directory. Current labels: config/current_md_labels.json. Historical graph methods: docs/ACTUAL_GRAPH_METHODS_20261111.md. Previous checkpoints: d5c34f1, b5ed245, c9f17af, b6f095f.

This document distinguishes completed results from pending stages. Live results below supersede the pending markers in the original new-study plan.

## Live results and progress

Updated 2026-11-11T18:50:58.997479+00:00.

- baseline, elob: 12/12 saved folds. Matched-pair AUROC 0.5000; balanced accuracy 0.5000; sensitivity 1.0000; specificity 0.0000. Independent CPU replay passed.
- baseline, vhl: 12/12 saved folds. Matched-pair AUROC 0.5833; balanced accuracy 0.5000; sensitivity 1.0000; specificity 0.0000. Independent CPU replay passed.
- temporal, elob: 12/12 saved folds. Matched-pair AUROC 0.0833; balanced accuracy 0.5000; sensitivity 1.0000; specificity 0.0000. Independent CPU replay passed.
- temporal, vhl: 12/12 saved folds. Matched-pair AUROC 0.7500; balanced accuracy 0.4583; sensitivity 0.9167; specificity 0.0000. Independent CPU replay passed.

Lysine extraction: 141/141 cases saved. Complete: 5,781 sampled frames.

# Final three-stage results

All three stages completed. Forty-eight GNN folds (24 original, 24 enriched) passed saved-weight replay and independent CPU inference. All 141 pose packets completed; eight fresh trajectory replays reproduced every saved measurement. Four contract tests passed. No new MD, docking or conventional classifier.

| Model | Arm | Matched-pair AUROC | Correct pairs | Balanced accuracy | Sensitivity | Specificity |
|---|---|---:|---:|---:|---:|---:|
| baseline | ELOB | 0.500 | 6/12 | 0.500 | 1.000 | 0.000 |
| baseline | VHL | 0.583 | 7/12 | 0.500 | 1.000 | 0.000 |
| temporal | ELOB | 0.083 | 1/12 | 0.500 | 1.000 | 0.000 |
| temporal | VHL | 0.750 | 9/12 | 0.458 | 0.917 | 0.000 |

## Did anything improve?

VHL: enriched minus baseline matched-pair AUROC +0.167; balanced accuracy -0.042. 4 previously wrong pairs became correct; 2 previously correct pairs became wrong.

ELOB: enriched minus baseline matched-pair AUROC -0.417; balanced accuracy +0.000. 0 previously wrong pairs became correct; 5 previously correct pairs became wrong.

The original v1 pooled LOPO AUROC was zero in both arms, with balanced accuracy 0.50. The new matched-pair estimate uses both held-out targets scored by the same model and constant training composition (5 positives, 1 negative). Its numerical difference from pooled LOPO does not establish a model improvement or prove class-prior bias caused the earlier result. Classification performance at 0.5 and the baseline/enriched comparison under the same matched protocol are the relevant separate checks. No score inversion or threshold optimization was used.

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

Trainable parameters: original 21,825; enriched 39,873.
