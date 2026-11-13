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
