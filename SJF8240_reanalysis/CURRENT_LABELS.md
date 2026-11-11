# Current SJF8240 MD labels and completed corrected run

The user confirmed EPHA2 is degraded. In the eight-target MD comparison:

| Class | Targets |
|---|---|
| Degraded /1 | DDR2, EPHA2, EPHB2, MAPK14_full, MET, RIPK2_full |
| Not degraded /0 | ABL1, AXL |

The actual retained system is MAPK14_full. Labels for the other non-MD targets and LF/PS measurements are not inferred or changed. The immutable user-label map is config/current_md_labels.json. Its training flag records the state when that label snapshot was created; current completion is in STATUS.json.

The corrected-label MD-only LOPO run is now complete in corrected_md_only_v1/. Both ELOB and true-VHL arms finished8folds with the original MD-only architecture,300epochs and seed42; all targets have MD graphs. Both arms have AUROC0.0,75% accuracy,50% balanced accuracy and0% specificity. All targets were called degraded, including ABL1/AXL. This run therefore does not distinguish degraded and non-degraded binders. The small cohort contains only2negative targets, and these exploratory results do not establish independent biological validation.

All16 saved models and32 prediction-only scores passed independent CPU replay, with maximum score difference5.96e-8. Current findings: corrected_md_only_v1/RESULTS.md. The earlier16-target dual encoder, EPHA2=0 labels, weights and metrics remain historical artifacts and were not overwritten. No new docking/MD or assay changes occurred.
