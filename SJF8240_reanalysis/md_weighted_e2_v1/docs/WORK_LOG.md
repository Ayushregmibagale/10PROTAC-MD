# Option2 feature integration

Read all existing E2 audit definitions and retained core packets. The legacy ABL1 delta80=-0.0466 was based on old images/legacy masks and a different aggregation. Reproduced all141 saved ligand-imaged E2 cases (201 frames each) against saved contact hashes, time axes, original topology/trajectory source hashes and per-pose summaries. Recomputed the two target scalars and matched reviewed target tables. Do not call these descriptors direct lysine burial or productive E2 probability.

The chosen features are unconditional equal-pose mean minimum E2/lysine NZ distance and mean within-pose true-VHL-heavy bridge-minus-nonbridge reach80 contrasts for poses with both states. Undefined single-state pose contrasts remain missing and are not set to zero. All eight target values are finite; paired-state coverage4–10poses is recorded. Use the same scalar values for ELOB/VHL graph controls. No pose selection uses labels or scores.

Append two z-normalized graph scalars to the128-dimensional MD GAT embedding before the original64-hidden classifier head (130input channels). GAT backbone AST remains unchanged; classifier head has128 additional weights. Balance class loss from the seven training targets exactly as Option1. E2 normalization mean/population standard deviation are fitted only on training-target features, saved in checkpoint/receipt, recomputed independently and reused at inference. Missing primary-target descriptors are fatal; prediction-only extension descriptors are unavailable, so those extension predictions are omitted explicitly. No conventional predictor is trained.

Eleven software tests include class-weight direction/gradient, held-out-label exclusion, original GAT backbone preservation, held-out-feature exclusion from normalization, zero-variance behavior, corrected raw graph labels and raw-index validation. Freeze all scripts/protocol/config/source feature provenance before fitting. No tuning based on any E2 run score.

## Completed evaluation

Option2 balanced loss plus audited E2 finished both8-fold LOPO arms with frozen settings. Every fold saved weights and passed same-runtime replay. Independent CPU replay, input/source hashes, training membership and class weights pass; Option2 also verifies training-only scaler values and scalar provenance. No run was tuned after seeing scores.

Primary VHL metrics: {"AUROC": 0.0, "accuracy": 0.25, "balanced_accuracy": 0.16666666666666666, "average_precision": 0.594047619047619, "sensitivity": 0.3333333333333333, "specificity": 0.0, "confusion_matrix_labels_0_1": [[0, 2], [4, 2]], "positive_count": 6, "negative_count": 2, "majority_class_accuracy_reference": 0.75}

ELOB control metrics: {"AUROC": 0.0, "accuracy": 0.125, "balanced_accuracy": 0.08333333333333333, "average_precision": 0.594047619047619, "sensitivity": 0.16666666666666666, "specificity": 0.0, "confusion_matrix_labels_0_1": [[0, 2], [5, 1]], "positive_count": 6, "negative_count": 2, "majority_class_accuracy_reference": 0.75}

Neither variant reliably distinguishes degraded and non-degraded binders in this cohort. Refer to the three-method comparison for the unchanged corrected-label unweighted baseline. These adverse results are retained, not re-labelled or inverted.
