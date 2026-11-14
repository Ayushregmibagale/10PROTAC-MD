# Main findings — 14 November 2026

The original seed-42 AUROC improvement (0.750 to 0.833) is not robust across the five fixed seeds. The control averages 0.733 (range 0.583–0.833); catalytic features average 0.667 (range 0.417–0.833). Features improve one seed, tie one and worsen three. No seed was discarded.

The fixed logistic baseline and geometry-only score each have matched-pair AUROC 0.500. This bounded comparison gives no evidence that the engineered catalytic summaries alone reliably discriminate the eight labelled targets. It does not prove that all geometry-based approaches must fail.

Threshold fitting on training proteins provides little benefit. Control balanced accuracy rises from 0.500 to 0.517 on average, with specificity only 0.033 (2 of 60 repeated negative evaluations). Catalytic balanced accuracy remains 0.492 with specificity 0.000. These repeats are correlated, not 60 distinct negatives. Thresholds were never selected using held-out targets; apparent training calibration is an explicit limitation because only one negative exists in each outer training set.

SIK3 and RPS6KB2: compatible older corrected MD-only v1 checkpoints predict both above 0.5 in all eight folds. Mean true-VHL scores are SIK3 0.579771 and RPS6KB2 0.579751. The model also struggles with known negatives, so this is weak, non-confirmatory evidence. These are not predictions from the new temporal catalytic model: compatible temporal/catalytic inputs for the unknown targets remain unavailable. Their aggregated old graphs cannot be substituted for missing new inputs.

EPHA2 remains a confirmed degrader (label 1). Unknown targets never enter training or performance calculations. Historical molecular identity is retained; no supplier-SJF8240 equivalence or experimental validation is claimed.

120 GNN models, 300 epochs each, all saved locally. Independent CPU replay checked 960 scores (all eight proteins per model, of which 240 scores are held out). Logistic saved-parameter replay differs by at most 3.04e-8, consistent with floating-point feature export. The initial baseline tolerance of 1e-10 was too strict for float32 feature export; this was recorded and changed to 1e-6 without refitting models or changing results. An earlier unknown-replay export failed because an imported ROOT variable shadowed the output path; it attempted no changes to original data, was corrected and rerun successfully.

Next methodological priority: more independently labelled non-degraders and prospective evaluation. More tuning on these eight proteins risks fitting their identities. The current data do not establish a reliable degrader/non-degrader classifier.
