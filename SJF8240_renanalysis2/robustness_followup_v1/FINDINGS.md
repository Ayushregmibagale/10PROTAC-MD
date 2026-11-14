# Main findings — 14 November 2026

The original seed-42 AUROC improvement (0.750 to 0.833) is not robust across the five fixed seeds. The control averages 0.733 (range 0.583–0.833); catalytic features average 0.667 (range 0.417–0.833). Features improve one seed, tie one and worsen three. No seed was discarded.

The fixed logistic baseline and geometry-only score each have matched-pair AUROC 0.500. This bounded comparison gives no evidence that the engineered catalytic summaries alone reliably discriminate the eight labelled targets. It does not prove that all geometry-based approaches must fail.

Threshold fitting on training proteins provides little benefit. Control balanced accuracy rises from 0.500 to 0.517 on average, with specificity only 0.033 (2 of 60 repeated negative evaluations). Catalytic balanced accuracy remains 0.492 with specificity 0.000. These repeats are correlated, not 60 distinct negatives. Thresholds were never selected using held-out targets; apparent training calibration is an explicit limitation because only one negative exists in each outer training set.

SIK3 and RPS6KB2: compatible older corrected MD-only v1 checkpoints predict both above 0.5 in all eight folds. Mean true-VHL scores are SIK3 0.579771 and RPS6KB2 0.579751. The model also struggles with known negatives, so this is weak, non-confirmatory evidence. These are not predictions from the new temporal catalytic model: compatible temporal/catalytic inputs for the unknown targets remain unavailable. Their aggregated old graphs cannot be substituted for missing new inputs.

EPHA2 remains a confirmed degrader (label 1). Unknown targets never enter training or performance calculations. Historical molecular identity is retained; no supplier-SJF8240 equivalence or experimental validation is claimed.

120 GNN models, 300 epochs each, all saved locally. Independent CPU replay checked 960 scores (all eight proteins per model, of which 240 scores are held out). Logistic saved-parameter replay differs by at most 3.04e-8, consistent with floating-point feature export. The initial baseline tolerance of 1e-10 was too strict for float32 feature export; this was recorded and changed to 1e-6 without refitting models or changing results. An earlier unknown-replay export failed because an imported ROOT variable shadowed the output path; it attempted no changes to original data, was corrected and rerun successfully.

Next methodological priority: more independently labelled non-degraders and prospective evaluation. More tuning on these eight proteins risks fitting their identities. The current data do not establish a reliable degrader/non-degrader classifier.


# Five-seed robustness, baselines and threshold assessment

EPHA2 remains label 1. SIK3 and RPS6KB2 are unknown and excluded from fitting/evaluation. Historical chemistry retained; no supplier equivalence claimed.

## Results

| model                   |   seed | policy        |   matched_pair_AUROC |   sensitivity |   specificity |   balanced_accuracy |
|:------------------------|-------:|:--------------|---------------------:|--------------:|--------------:|--------------------:|
| geometry_joint_fraction |     42 | fixed_0.5     |             0.5      |      0        |      1        |            0.5      |
| geometry_joint_fraction |     42 | training_only |             0.5      |      0.416667 |      0.583333 |            0.5      |
| logistic_catalytic      |     42 | fixed_0.5     |             0.5      |      0.75     |      0        |            0.375    |
| logistic_catalytic      |     42 | training_only |             0.5      |      0.75     |      0        |            0.375    |
| native_features         |      7 | fixed_0.5     |             0.833333 |      1        |      0        |            0.5      |
| native_features         |      7 | training_only |             0.833333 |      1        |      0        |            0.5      |
| native_features         |     42 | fixed_0.5     |             0.833333 |      1        |      0        |            0.5      |
| native_features         |     42 | training_only |             0.833333 |      1        |      0        |            0.5      |
| native_features         |    123 | fixed_0.5     |             0.583333 |      1        |      0        |            0.5      |
| native_features         |    123 | training_only |             0.583333 |      1        |      0        |            0.5      |
| native_features         |   2026 | fixed_0.5     |             0.416667 |      1        |      0        |            0.5      |
| native_features         |   2026 | training_only |             0.416667 |      1        |      0        |            0.5      |
| native_features         |  31415 | fixed_0.5     |             0.666667 |      0.916667 |      0        |            0.458333 |
| native_features         |  31415 | training_only |             0.666667 |      0.916667 |      0        |            0.458333 |
| zero_padded_control     |      7 | fixed_0.5     |             0.833333 |      1        |      0        |            0.5      |
| zero_padded_control     |      7 | training_only |             0.833333 |      1        |      0.166667 |            0.583333 |
| zero_padded_control     |     42 | fixed_0.5     |             0.75     |      1        |      0        |            0.5      |
| zero_padded_control     |     42 | training_only |             0.75     |      1        |      0        |            0.5      |
| zero_padded_control     |    123 | fixed_0.5     |             0.666667 |      1        |      0        |            0.5      |
| zero_padded_control     |    123 | training_only |             0.666667 |      1        |      0        |            0.5      |
| zero_padded_control     |   2026 | fixed_0.5     |             0.583333 |      1        |      0        |            0.5      |
| zero_padded_control     |   2026 | training_only |             0.583333 |      1        |      0        |            0.5      |
| zero_padded_control     |  31415 | fixed_0.5     |             0.833333 |      1        |      0        |            0.5      |
| zero_padded_control     |  31415 | training_only |             0.833333 |      1        |      0        |            0.5      |

## Interpretation and method

Five predeclared seeds: 42, 7, 123, 2026, 31415. Same 12 outer pair splits, 300 epochs, unweighted BCE, architecture and optimizer as prior ablation. No seed, hyperparameter or feature selection based on outcomes. All checkpoints saved and independently replayed on CPU.

Logistic regression: 15 fixed catalytic summaries (five lysine means across windows, temporal standard deviations and end-minus-start); StandardScaler fitted only on outer training targets, L2 C=1. Geometry score is mean simultaneous exposed/reachable/scaffold-compatible/bridging fraction (feature 5); higher is fixed to degrader-like. The geometry score is not a calibrated probability, so its 0.5 threshold is only a numerical comparator.

Threshold maximizes balanced accuracy on the six training proteins; ties choose closest to 0.5 then smaller threshold. Outer held-out labels never select thresholds. This is training-resubstitution calibration, not nested CV: outer training has only one negative, so a non-degenerate inner negative holdout is unavailable. Training fit can be optimistic and score calibration may shift to unseen proteins.

Eight independent labelled proteins and two negatives; pairs, seeds and repeated predictions are correlated, not new biological replicates. Report seed variability descriptively without a naive significance test.

## Across-seed summary

| model                   | policy        |   AUROC_mean |   AUROC_min |   AUROC_max |   balanced_accuracy_mean |   sensitivity_mean |   specificity_mean |
|:------------------------|:--------------|-------------:|------------:|------------:|-------------------------:|-------------------:|-------------------:|
| geometry_joint_fraction | fixed_0.5     |     0.5      |    0.5      |    0.5      |                 0.5      |           0        |          1         |
| geometry_joint_fraction | training_only |     0.5      |    0.5      |    0.5      |                 0.5      |           0.416667 |          0.583333  |
| logistic_catalytic      | fixed_0.5     |     0.5      |    0.5      |    0.5      |                 0.375    |           0.75     |          0         |
| logistic_catalytic      | training_only |     0.5      |    0.5      |    0.5      |                 0.375    |           0.75     |          0         |
| native_features         | fixed_0.5     |     0.666667 |    0.416667 |    0.833333 |                 0.491667 |           0.983333 |          0         |
| native_features         | training_only |     0.666667 |    0.416667 |    0.833333 |                 0.491667 |           0.983333 |          0         |
| zero_padded_control     | fixed_0.5     |     0.733333 |    0.583333 |    0.833333 |                 0.5      |           1        |          0         |
| zero_padded_control     | training_only |     0.733333 |    0.583333 |    0.833333 |                 0.516667 |           1        |          0.0333333 |

## Paired seed changes

|   seed |   native_features |   zero_padded_control |   native_minus_control |
|-------:|------------------:|----------------------:|-----------------------:|
|      7 |          0.833333 |              0.833333 |              0         |
|     42 |          0.833333 |              0.75     |              0.0833333 |
|    123 |          0.583333 |              0.666667 |             -0.0833333 |
|   2026 |          0.416667 |              0.583333 |             -0.166667  |
|  31415 |          0.666667 |              0.833333 |             -0.166667  |

Catalytic ranking improves in 1/5 seeds, ties in 1, worsens in 3. No best-seed selection is used.

# SIK3 and RPS6KB2 prediction-only assessment

| model                       | role   | unknown_target   |   mean_score |   min_score |   max_score |   degrader_vote_fraction |
|:----------------------------|:-------|:-----------------|-------------:|------------:|------------:|-------------------------:|
| legacy_corrected_md_only_v1 | elob   | RPS6KB2_full     |     0.562268 |    0.52303  |    0.659893 |                        1 |
| legacy_corrected_md_only_v1 | elob   | SIK3_full        |     0.561903 |    0.523786 |    0.659404 |                        1 |
| legacy_corrected_md_only_v1 | vhl    | RPS6KB2_full     |     0.579751 |    0.541105 |    0.634672 |                        1 |
| legacy_corrected_md_only_v1 | vhl    | SIK3_full        |     0.579771 |    0.540586 |    0.634567 |                        1 |

These are predictions from the older corrected MD-only v1 model, independently replayed on CPU. Neither unknown was labelled or used in fitting. This is not experimental degradation evidence. ELOB is an alternative historical arm, not true VHL.

The AUROC 0.8333 temporal catalytic model cannot yet score these targets: their matching four-window 29-feature graphs and five catalytic-geometry annotations have not been built. Existing aggregated 3-feature graphs are not compatible. No missing catalytic features were silently zero-filled. Newest-model unknown classification remains unavailable until matching inputs are prepared and verified from existing trajectories.
