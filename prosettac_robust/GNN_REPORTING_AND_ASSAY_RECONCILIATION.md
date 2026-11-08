# Retained GNN reporting and assay correspondence — 2026-11-08

## What was checked and why

Read-only audits bind saved score tables, draft documents, the retained primary assay workbook, cohort accessions and prior structure-identity evidence to SHA-256 hashes. No docking, MD or biological model fitting occurred. Existing labels and scores remain unchanged. This prevents manuscript claims from drifting away from the files actually available for review.

## Retained E2-GNN scores

| Cohort | Positive / negative labels | Correct calls | AUROC from saved rounded scores |
|---|---|---|---|
| All 17 targets | 5 / 12 | 15 / 17 | 0.800 |
| Eight targets with MD | 5 / 3 | 6 / 8 | 0.600 |
| Nine targets without MD | 0 / 9 | 9 / 9 | Undefined |

A single-class subgroup has no positive-negative pairs for AUROC. Nine correct negative calls do not establish sensitivity, perfect discrimination or a result that generalizes to any sample size. The original robustness draft's no-MD AUROC=1.000 wording is unsupported. No fixed target count alone establishes publication quality. These are descriptive results on historical labels, not a new biological evaluation.

The isolated runner previously evaluated both subgroup AUROCs in one try/except block. On versions that raise for a single-class subgroup, that could erase the valid MD-group result too. New independent subgroup evaluation returns an explicit undefined status and class counts. The TSV keeps NaN for undefined AUROC; it never substitutes zero or one. Invalid labels, nonfinite scores and mismatched masks fail explicitly. Two new regression tests cover a valid group alongside single-class/empty groups and malformed input. Architecture and historical predictions are unchanged.

## Which baseline can be compared

The current retained `dual_lopo_results.tsv` has 45 rows, including 28 unknown-labelled extensions/candidates. All 17 reviewed historical rows have score 0.0 and has_md=0. On those rows its current-score AUROC is 0.500 and accuracy 12/17. That file therefore cannot reconstruct the older draft's 04d AUROC=0.850, accuracy 14/17, or a paired architecture improvement comparison. This does not prove the historical 0.850 run never occurred. Exact original per-fold scores, input snapshots and weights would be needed; do not silently regenerate them or present a fresh run as historical recovery.

The retained E2 table's has_e2 flag is not a certificate that every scalar is finite: MAP4K5 and SLK have flag 1 with missing mean-distance values. Per-field missingness must remain explicit. The earlier availability audit already establishes selection confounding; no new metadata predictor was trained here. Architecture search and rationale development on this same small cohort are exploratory. Molecular/target-family generalization, causal feature attribution and calibrated degradation probabilities remain unestablished.

## Accession-level assay check

A fresh accession join confirms 11/17 historical targets have both compound-1 abundance measurements in the retained workbook. Missing targets remain missing, not negative. The measured CDK17/Q00537 row is Excel row 49. The retained simulated structure matches CDK14/O94921; its separate workbook row 114 is named PFTK1 and says Not Determined for abundance. Binding-column values do not supply a degradation endpoint. The CDK17 measurement cannot be transferred to the CDK14 structure. Neither original structure nor label was changed or dropped.

The previously recorded ABL1 endpoint discordance and MAP4K5 disagreement remain unresolved. An assay hierarchy must be selected independently of scores. Model disagreement does not establish an experimental false negative. The external article browser returned a challenge and the alternate full-text endpoint failed during this pass; no new article claim was inferred from those failures. This audit uses the already retained, hash-verified workbook and prior evidence.

## Evidence and limits

`existing_result_audit/gnn_reporting_v1` contains source snapshots, subgroup metrics and flag semantics. `existing_result_audit/assay_accession_v1` contains accession joins, source Excel row numbers and hashes. Software checkpoint v10 records 141 analysis passes and 137 docking passes / four PyG-only skips. Passing software checks does not clear assay, target-identity, exact historical GNN replay or independent validation gates. No new scientific sign-off is issued.
