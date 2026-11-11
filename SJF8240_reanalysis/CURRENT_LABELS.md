# Current MD labels and evaluated interventions

Latest experiment: corrected_md_only_v2 adds SIK3_full and RPS6KB2_full as hypothetical negatives only. Their experimental labels remain unknown. This does not change the confirmed eight-target label map below. Both ten-target LOPO arms completed; neither arm predicts either hypothetical target as a non-degrader. See corrected_md_only_v2/RESULTS.md for assumption-scoped metrics and held-out scores.

| Class | Targets |
|---|---|
| Degraded / 1 | DDR2, EPHA2, EPHB2, MAPK14_full, MET, RIPK2_full |
| Not degraded / 0 | ABL1, AXL |

These are the user-confirmed labels for eight MD evaluation targets. The retained system is MAPK14_full. Non-MD target labels and original LF/PS measurements are not inferred or changed. config/current_md_labels.json is an immutable provenance snapshot; its creation-time training flag is historical. Current execution state is STATUS.json.

The corrected-label MD-only baseline, fold-balanced loss and fold-balanced loss plus E2 variants have all completed. See fix_comparison_20261111/COMPARISON.md for the controlled comparison. Neither tested intervention gives useful discrimination: the primary VHL specificity remains 0% in both requested fixes. All old runs, weights and metrics remain preserved, with no post-score relabelling or inversion.

E2 inputs use audited ligand-imaged geometry and paired true-VHL-heavy state contrasts, not legacy values or a claim that all lysines are buried. Independent CPU replay verifies all 32 model folds in the two requested interventions. This small retrospective cohort is not independent biological validation.
