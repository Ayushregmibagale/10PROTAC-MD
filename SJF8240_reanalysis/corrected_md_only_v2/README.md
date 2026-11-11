# corrected_md_only_v2: ten-target hypothetical-negative test

Exact corrected_md_only_v1 model bytes, architecture,300epochs,base seed42 with existing target-bound fold seeds, optimizer, unweighted BCE and threshold. Six degraded targets remain DDR2, EPHA2, EPHB2, MAPK14_full, MET and RIPK2_full. ABL1/AXL remain known/user-labelled negatives. SIK3_full/RPS6KB2_full are assumed negatives only for this exploratory version; their experimental degradation labels remain unknown.

The two existing dynamic graph sources are copied without rebuilding geometry; four local NPZ copies change only label metadata from-1 to0. All other inputs reuse v1 bytes. Each of ten LOPO folds trains on the other nine targets. Both ELOB and VHL arms save weights and same-runtime replay receipts, followed by independent CPU verification. No weighted loss or E2 scalar branch is used.

The report compares v1(6:2) and v2(6:4), gives v2 scores on the original eight targets, explicit hypothetical-target held-out scores, training exclusion flags and pair-count decomposition. Agreement with assumed negative labels is not experimental validation. Cohort composition changes along with class balance; improvement/failure does not prove general geometric separability/impossibility. Check STATUS.json for actual progress.
