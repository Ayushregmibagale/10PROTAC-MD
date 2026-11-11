# Reviewed SJF8240 paired reanalysis

Both frozen model arms completed all 16 LOPO folds. This is a new exploratory comparison on retained structures and historical labels; it does not reconstruct the missing historical model or establish independent biological validation.

| Model arm | Overall AUROC (16 targets) | Accuracy | MD-only AUROC (8 targets) |
|---|---:|---:|---:|
| ELOB control | 0.7273 | 75.0% | 0.0000 |
| True VHL | 0.7273 | 68.8% | 0.0000 |

The full cohort has five positive and eleven negative historical labels. The MD subgroup contains all five positives and three negatives. The eight targets without MD are all negative, so their subgroup AUROC is undefined. MD availability is associated with the labels: aggregate performance alone is insufficient evidence of biological discrimination. These subgroup checks do not establish which feature caused a prediction, and no additional model fitting or tuning was performed.

All 307,770 configured PDB files are accounted for: 307,767 complete files included, three unusable files excluded symmetrically with original hashes and reasons. All 19 static targets and 174 MD cases passed artifact checks. Each MD case contributes 201 unique frames over 0–20 ns; all ten target graph aggregations were independently reproduced from saved contact packets. Original structures, labels, trajectories and historical model outputs remain preserved. No new docking or MD was run.

All 32 saved fold models and 64 prediction-only fold inferences passed independent CPU replay against saved GPU scores. Maximum score difference: ELOB 5.96e-08; VHL 1.19e-07. Frozen input/source hashes, output/checkpoint hashes, held-out membership and recomputed metrics passed verification. This is software/artifact verification, not independent experimental validation.

CDK17-labelled structures match CDK14 and were excluded from the primary cohort before fitting, without renaming or relabelling originals. SIK3 and RPS6KB2 remain prediction-only; fold-ensemble scores are uncalibrated and are not experimental labels. LF abundance, PS light and PS heavy endpoints stay separate, with all source endpoint rows preserved. Historical assay-label correspondence, supplier equivalence, external generalization and calibration remain limitations.

The frozen model uses the original 04d raw-graph dual encoder, 400 epochs, seed 42 and no tuning. E2 scalar descriptors are separate audit outputs and are not inputs to this model. See `analysis/comparison/` for paired target scores, input coverage, endpoint exports and prediction-only summaries; `validation/` holds independent verification receipts and `validation/reviewed_metrics.tsv` records subgroup counts.
