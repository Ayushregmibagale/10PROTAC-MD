# Class weighting and E2 GNN comparison

Both requested interventions are implemented and evaluated using the same eight MD targets, user-confirmed six-versus-two labels, raw graph inputs, 300 epochs and seed42. Only ABL1/AXL are non-degraded; EPHA2 is degraded. Every target has MD. Original versions, data, weights and results are preserved.

| Method (true VHL) | AUROC | Accuracy | Balanced accuracy | Sensitivity | Specificity |
|---|---:|---:|---:|---:|---:|
| Original corrected-label MD-only | 0.0000 | 75.0% | 50.0% | 100.0% | 0.0% |
| Option 1: balanced loss | 0.0833 | 37.5% | 25.0% | 50.0% | 0.0% |
| Option 2: balanced loss + E2 | 0.0000 | 25.0% | 16.7% | 33.3% | 0.0% |

Option 1 uses balanced per-target BCE weights N_train/(2*N_class) from the seven training targets only. With degrader=1, pos_weight=3 would emphasize degraders, so it is not used. Negative/positive weight ratios are6 or2.5 depending on the held-out class. Tests verify equal total class weighting, the gradient direction and exclusion of held-out labels. This balances optimization; it does not supply more negative examples.

Option 2 appends audited mean minimum E2/lysine distance and delta reach80 to the MD GNN embedding before classification. Means/std are fitted on the seven training targets only, recorded per fold and independently checked. The GAT backbone remains unchanged; the classifier first layer adds two inputs (128 additional weights). Both role controls receive identical scalar values. The feature hypothesis was motivated by already-known target behavior, so this retrospective LOPO is not an independent validation of that feature-selection hypothesis.

The corrected E2 scalars reproduce all141 existing core cases and28,341 frames. They use a ligand-connected coherent image and true-VHL heavy-contact bridge mask. ABL1 delta80 is approximately-0.036848 under this definition; the legacy-0.0466 value uses older definitions. Mean minimum distance is unconditional over frames; delta80 is an average of paired within-pose contrasts (4–10paired poses/target). These are modeled geometric reachability proxies, not proof that all18 lysines are buried, a full E2 steric model or catalytic probability. Source candidate/imaging limitations remain explicit.

The six-versus-two cohort is small and imbalanced; the majority-class accuracy reference is75%. Each negative-held-out fold has just one negative in training. Read accuracy alongside AUROC, balanced accuracy and specificity. No threshold inversion, new seed search, label changes or outcome-based tuning was performed. A score improvement on this cohort would not establish independent biological validity.

Both fixes ran with paired ELOB controls; metrics.tsv contains all six model/method rows and confusion matrices. per_target_predictions.tsv retains full precision held-out scores. All completed model arms pass saved-weight replay, source/input hash checks, exact held-out membership and independent CPU metric recomputation. E2 scaler lineage and fold-local class weights are checked explicitly. E2-augmented extension predictions are omitted because audited E2 descriptors for SIK3/RPS6KB2 are unavailable; no fabricated zeros are used. No new docking, MD or conventional predictor was run.

PyTorch loss definition: https://docs.pytorch.org/docs/2.10/generated/torch.nn.BCEWithLogitsLoss.html
