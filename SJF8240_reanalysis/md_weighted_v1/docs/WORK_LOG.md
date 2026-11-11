# Option1: class balance

User requested emphasis on the two non-degraded targets. Verified PyTorch2.10 BCEWithLogitsLoss documentation: pos_weight multiplies the positive term; pos_weight3 would emphasize the six degraded targets. Implemented balanced per-target weights N_train/(2*N_class), computed independently inside each seven-target LOPO training fold. The negative/positive ratios are6 or2.5 depending on the held-out class; no held-out label enters the computation. All old graphs/labels/model architecture and300-epoch seed42 protocol are preserved. Only loss weighting changes.

This fix changes optimization contributions, not sample size. It cannot guarantee learning the negative class from two proteins. No conventional predictive model, new docking or MD is run. E2 features are reviewed separately for Option2.

## Completed evaluation

Option1 class balance finished both8-fold LOPO arms with frozen settings. Every fold saved weights and passed same-runtime replay. Independent CPU replay, input/source hashes, training membership and class weights pass; Option2 also verifies training-only scaler values and scalar provenance. No run was tuned after seeing scores.

Primary VHL metrics: {"AUROC": 0.08333333333333333, "accuracy": 0.375, "balanced_accuracy": 0.25, "average_precision": 0.6218253968253968, "sensitivity": 0.5, "specificity": 0.0, "confusion_matrix_labels_0_1": [[0, 2], [3, 3]], "positive_count": 6, "negative_count": 2, "majority_class_accuracy_reference": 0.75}

ELOB control metrics: {"AUROC": 0.3333333333333333, "accuracy": 0.25, "balanced_accuracy": 0.3333333333333333, "average_precision": 0.6968253968253968, "sensitivity": 0.16666666666666666, "specificity": 0.5, "confusion_matrix_labels_0_1": [[1, 1], [5, 1]], "positive_count": 6, "negative_count": 2, "majority_class_accuracy_reference": 0.75}

Neither variant reliably distinguishes degraded and non-degraded binders in this cohort. Refer to the three-method comparison for the unchanged corrected-label unweighted baseline. These adverse results are retained, not re-labelled or inverted.
