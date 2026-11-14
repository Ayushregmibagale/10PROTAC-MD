# SJF8240 robustness follow-up

Authorized by user: run five-seed robustness, simple baselines, training-only thresholds; also check SIK3/RPS6KB2 prediction-only. No original analysis is overwritten.

Run with the existing miniconda Python runtime, CPU native PyG. Numbered scripts:

1. `01_run.py --workers 4`: train and save 120 matched GNN folds (5 seeds × 12 pairs × 2 feature variants), 300 epochs each. Existing 8-target protocol retained; EPHA2=1. Both variants match initial weights within seed/pair. Freeze checks guard original inputs and training/report code.
2. `02_compare_verify.py`: independent CPU replay of all eight protein scores per checkpoint, matched pair metrics at 0.5 and training-only thresholds, fixed regularized logistic and geometry baselines with identical held-out splits. Training scores are used only for threshold fitting, never counted as held-out performance.
3. `03_unknown_legacy.py`: independent replay of SIK3/RPS6KB2 predictions from the compatible older corrected MD-only v1 checkpoints. Never train on hypothetical labels. New temporal catalytic inputs are unavailable for these two targets; older graph schemas must not be silently substituted. Preserve distinction between model predictions and experimental status.
4. `04_finalize.py`: replay baseline saved coefficients/scaler, summarize seed variability without treating correlated repeats as biological replicates, write step-by-step Markdown/Word findings and validation receipt.

Threshold calibration fits apparent training scores because each outer training set has only one negative. No valid stratified inner negative holdout is available. Thus these are outer-held-out tests of a resubstitution threshold-fitting procedure, not nested cross-validation. No threshold is selected against the held-out pair.

Weights are saved locally; Git keeps scripts, compact inputs/parameters, reports and receipts, not checkpoint binaries. The manifest records original absolute source paths; moving scripts requires preserving/rebinding provenance explicitly. Historical chemistry limitation remains. SIK3/RPS6KB2 are unknown, not validation labels.
