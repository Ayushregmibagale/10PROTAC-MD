# Saved GNN: availability control audit

No model was fitted or retrained. This reads the existing 17-target prediction table and compares its retained out-of-fold scores with two fixed metadata flags. Historical labels and all original outputs remain unchanged; assay/chemical-identity concerns are not resolved by this exercise.

| Recorded descriptor | Correct / 17 | AUROC | Balanced accuracy |
|---|---:|---:|---:|
| Saved GNN | 15 | 0.8000 | 0.8583 |
| MD available flag | 14 | 0.8750 | 0.8750 |
| E2 available flag | 12 | 0.7917 | 0.7917 |

All five recorded positives have MD; nine targets without MD are all recorded negatives. Eight targets have MD (five positives/three negatives). E2 availability has a similar imbalance. The saved GNN AUROC is recalculated from its retained rounded scores, without reconstructing unavailable historical unrounded logits or weights.

These are retrospective descriptive controls, not preregistered independent benchmarks or calibrated probabilities. Availability can encode a legitimate stage-selection outcome in a cascaded workflow; it is not automatically leakage. However, here the inherited screening/contact roles and experimental labels require review, so one must distinguish feasibility evidence from unexamined missingness/selection patterns.

The saved GNN makes one more correct classification than the MD flag on these cases, while its ranking AUROC is lower. With only 17 selected targets, that cannot establish general incremental value, mechanistic learning or lack of usefulness. Preserve the result as exploratory. Any claim of added value should include simple availability/scaffold baselines, appropriately grouped evaluation and independently reviewed experimental labels. A corrected-feature evaluation is separate from the historical reported model and need not require new MD.

Source hashes, metadata contingency counts and calculations are retained in `gnn_availability_audit_v1/`; executable audit source is `review_tools/audit_saved_gnn_availability.py`.
