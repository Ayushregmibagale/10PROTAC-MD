# SJF8240 analysis: step-by-step record

Last updated: 2026-11-11. New three-stage investigation is IN PROGRESS. Original source trajectories and earlier results are preserved.

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

## 7. New evaluation: leave-pair-out — pending
Hold out one degraded and one non-degraded target together: 12 pairs per arm. Every training fold has five positives and one negative. Compare held-out scores within the same model; do not pool scores from differently trained models as an AUC. Report fraction of correctly ordered pairs, sensitivity and specificity at fixed 0.5. Pairs share targets and are not independent biological replicates. Same original model, 300 epochs, seed 42 and unweighted loss.

## 8. Per-lysine mechanistic analysis — pending
Use 41 existing frames per pose (0, 0.5, ..., 20 ns), equal pose weighting. Measure lysine NZ solvent-accessible area in the coherent complex and isolated target, E2-point distance and exposed/reachable bridge occupancy. Four fixed time windows retain coarse dynamics. No label-selected thresholds. These are geometric proxies; full E2–ubiquitin catalytic orientation is not available from a point cloud.

## 9. Lysine-aware temporal GNN — pending
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

Updated 2026-11-11T18:08:30.635106+00:00.

- baseline, elob: 12/12 saved folds. Matched-pair AUROC 0.5000; balanced accuracy 0.5000; sensitivity 1.0000; specificity 0.0000. Independent CPU replay passed.
- baseline, vhl: 12/12 saved folds. Matched-pair AUROC 0.5833; balanced accuracy 0.5000; sensitivity 1.0000; specificity 0.0000. Independent CPU replay passed.
- temporal, elob: 3/12 saved folds. CPU verification pending.
- temporal, vhl: 9/12 saved folds. CPU verification pending.

Lysine extraction: 141/141 cases saved. Complete: 5,781 sampled frames.
