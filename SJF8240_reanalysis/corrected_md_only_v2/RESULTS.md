# corrected_md_only_v2: exploratory hypothetical-negative cohort

Same v1 MD-only model, original architecture, 300 epochs, base seed 42, target-bound fold seeds, unweighted BCE loss, optimizer and threshold. All ten targets have retained MD graphs. Six degraded targets: DDR2, EPHA2, EPHB2, MAPK14_full, MET and RIPK2_full. Confirmed/user-labelled non-degraders: ABL1 and AXL. SIK3_full and RPS6KB2_full are assumed non-degraders for this version only; their experimental degradation labels remain unknown.

| Arm/scope | AUROC | Balanced accuracy | Sensitivity | Specificity | Accuracy |
|---|---:|---:|---:|---:|---:|
| elob: v1_original_eight | 0.0000 | 50.0% | 100.0% | 0.0% | 75.0% |
| elob: v2_all_ten_assumed_labels | 0.0417 | 41.7% | 83.3% | 0.0% | 50.0% |
| elob: v2_original_eight_subset | 0.0000 | 41.7% | 83.3% | 0.0% | 62.5% |
| vhl: v1_original_eight | 0.0000 | 50.0% | 100.0% | 0.0% | 75.0% |
| vhl: v2_all_ten_assumed_labels | 0.1667 | 12.5% | 0.0% | 25.0% | 10.0% |
| vhl: v2_original_eight_subset | 0.3333 | 25.0% | 0.0% | 50.0% | 12.5% |

## Hypothetical targets: genuine held-out predictions

| Arm | Target | Held-out score | Predicted class | Agrees with assumed non-degrader label | Excluded from own training |
|---|---|---:|---:|---|---|
| elob | RPS6KB2_full | 0.55559129 | 1 | False | True |
| elob | SIK3_full | 0.52506757 | 1 | False | True |
| vhl | RPS6KB2_full | 0.52445984 | 1 | False | True |
| vhl | SIK3_full | 0.55877954 | 1 | False | True |

Each model trained on the other nine proteins. When SIK3 is held out, its own assumed label does not enter its fitting; RPS6KB2 remains in that training fold under its assumed label, and vice versa. Independent verification checks actual saved training membership, labels/provenance, source/input/weight hashes, held-out score replay and metrics. All 20 saved models were replayed independently on CPU. Agreement against a hypothetical label is not biological correctness or experimental validation.

The overall negative count increases from two to four and the majority-class accuracy reference falls from 75% to 60%. The comparison above reports whether aggregate results change and separately evaluates v2 scores on the original eight targets, so changes driven solely by the two assumed negatives remain visible. Rank-pair decomposition separates the 12 positive-versus-confirmed-negative pairs from the 12 positive-versus-hypothetical-negative pairs. Adding targets changes geometry, cohort composition and optimizer updates per epoch as well as class balance; it does not isolate a causal effect of balance alone.

A higher AUROC would show consistency with this assumed partition in this retrospective test, not confirm the degradation labels of SIK3/RPS6KB2 or prove general geometric discrimination. Failure would show that this particular frozen GNN/cohort failed; it would not prove MD contact graphs incapable of separation under every model or class balance. No new seed search, score inversion, threshold tuning or relabelling was performed after outcomes. No new docking/MD or graph geometry rebuild occurred: four private NPZ copies change only unknown label metadata from -1 to the hypothesis value 0. Original graph files and the v1 framework/results remain unchanged.

Full-precision scores, per-fold weights/receipts, hypothetical-target flags and both independent CPU audits are retained. No class-weight or E2 intervention is used in v2.
