# Interpretation of the v2 cohort test

The unweighted-loss and full target-LOPO procedure are intentionally unchanged from v1. Holding out a positive leaves five positives and four negatives (55.6% positive) in training. Holding out a negative leaves six positives and three negatives (66.7% positive). A predictor dominated by training prevalence could therefore rank negative-held-out targets above positive-held-out targets. This is a conditional explanation, not proof that the fitted GNN uses only prevalence. The saved fold_class_priors.tsv lets this be inspected directly without another fitted model or any score adjustment.

Moving from 6:2 to 6:4 does not remove every imbalance-related issue: positives remain the majority in every training fold. Cohort composition, graph geometry and the number of per-epoch graph updates also change. The hypothetical labels of SIK3/RPS6KB2 are assumptions used in this version, not new experimental truth. Their held-out scores test consistency with those assumptions; they do not validate the assumptions biologically.

Use the full ten-target metrics together with the original-eight subset, hypothetical-target table and rank-pair decomposition. A higher aggregate metric may be driven by the two assumed negatives. Success or failure of this specific frozen GNN cannot prove general geometric separability or impossibility for all MD contact-graph methods. No additional fitting, threshold change, score inversion or model selection occurs in this descriptive audit.

## Observed result of the proposed balance change

V1 had AUROC 0.0, balanced accuracy 50%, sensitivity 100% and specificity 0% in both roles. V2 full-cohort ELOB metrics are AUROC 0.0417, balanced accuracy 41.7%, sensitivity 83.3% and specificity 0%. V2 VHL metrics are AUROC 0.1667, balanced accuracy 12.5%, sensitivity 0% and specificity 25%. The numbers change, but the 6:4 hypothetical cohort does not rescue useful classification in this frozen test.

On the original eight targets, v2 ELOB AUROC remains 0.0; v2 VHL AUROC is 0.3333. These subset results come from v2 models trained on nine other targets, not refits on the original seven, so they describe sensitivity to the cohort change rather than an independent replication of v1.

Neither arm predicts either hypothetical target as class 0: ELOB SIK3/RPS6KB2 scores are 0.52506757/0.55559129; VHL scores are 0.55877954/0.52445984. Every relevant saved training list excludes the target being scored. These outcomes disagree with the imposed negative assumptions; they do not demonstrate actual degradation or non-degradation. The unconfirmed targets must not be added to an experimentally validated negative list.
