# Corrected-label MD-only GNN results

Completed a new eight-target MD-only LOPO comparison with the user-confirmed EPHA2 label. Six degraded binders: DDR2, EPHA2, EPHB2, MAPK14_full, MET and RIPK2_full. Two non-degraded binders: ABL1 and AXL. Every evaluated target has an MD graph; static docking graphs and E2 scalar features are absent from the model.

| Model arm | AUROC | Accuracy | Balanced accuracy | Sensitivity | Specificity |
|---|---:|---:|---:|---:|---:|
| ELOB control | 0.0000 | 75.0% | 50.0% | 100.0% | 0.0% |
| True VHL (primary) | 0.0000 | 75.0% | 50.0% | 100.0% | 0.0% |

The cohort has only two negative targets. An always-degraded call would achieve 75% accuracy; interpret accuracy together with AUROC, balanced accuracy, sensitivity and specificity. These are exploratory internal LOPO results, not independent biological validation or calibrated degradation probabilities. Each negative-held-out fold trains on only one other negative target. No optimization was performed after seeing scores.

The frozen original 04h MD-only model uses a two-layer GAT with mean/max pooling and its original classifier head, hidden size64, four heads, dropout0.3, 300 epochs, Adam lr0.001/weight decay0.0001, cosine schedule and seed42 with target-bound fold seeds. The original model class AST was verified exactly. Training holds out the whole target; saved weights and inference replay are recorded for every fold. Independent CPU replay checks all16 fold models and32 prediction-only fold inferences. Maximum CPU score differences: ELOB 5.96e-08, VHL 5.96e-08.

This version reuses the existing MD contact graphs, correcting EPHA2 label metadata without changing graph features. No new docking or MD was run. Old16-target dual-encoder results remain frozen under their old labels. This new comparison changes cohort, labels and architecture together and cannot isolate any one change relative to that historical run. Within this version the ELOB/VHL arms use the same model, target labels and training settings.

The retained MD sources cover141 systems across the eight evaluation targets. Prediction-only SIK3 and RPS6KB2 reuse33 additional systems; their fold-ensemble scores are exported separately and are not validation labels. The retained structure is MAPK14_full; no MAPK15 structure was substituted. Historical experimental endpoint correspondence, retained chemistry/simulation heterogeneity, external generalization and calibration remain limitations. LF/PS assay measurements and all prior models, metrics and source data are unchanged.

`results/paired_predictions.tsv` gives all held-out scores; `results/completion.json` includes confusion matrices and average precision. Protocol/source/input hashes, per-fold weights, prediction-only scores and independent review receipts are retained in this version.
