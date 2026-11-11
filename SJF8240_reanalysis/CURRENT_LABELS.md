# Current SJF8240 MD labels — user correction

EPHA2 is degraded (label 1), as explicitly confirmed by the user on 11 November 2026 after the completed historical-label run. In the eight-target MD comparison, ABL1 and AXL are the only non-degraded targets.

| Label | Targets |
|---|---|
| Degraded / 1 | DDR2, EPHA2, EPHB2, MAPK14_full, MET, RIPK2_full |
| Not degraded / 0 | ABL1, AXL |

The actual retained system name remains MAPK14_full. This correction applies to the MD comparison; labels for other non-MD targets and original LF/PS assay measurements are not inferred or changed.

The current label map is `config/current_md_labels.json`. Versioned EPHA2 static and dynamic graph copies, for both ELOB and VHL, contain label 1 with all non-label arrays exactly preserved. `label_updates/20261111_epha2/manifest.json` binds the source and derivative hashes and references the unchanged graphs for the other seven targets. Overlay paths in that manifest are relative to its version folder. The builder accepts explicit source, label-map and output paths and refuses to overwrite an existing overlay.

The frozen protocol, original graph files, saved weights, predictions, metrics and reviewed results from checkpoint d5c34f1 remain unchanged as historical artifacts. They used EPHA2 = 0. Their reported aggregate AUROC 0.727 and MD-subgroup AUROC 0.0 do not describe the user-corrected label set. Most LOPO training folds included EPHA2, so merely changing its label in a prediction table would not produce corrected-label LOPO results. A separately versioned fit/evaluation is needed. No new training, docking or MD has been executed for this label update.
