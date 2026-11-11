# SJF8240 retained-structure reanalysis

The latest completed run is the [corrected-label MD-only GNN](corrected_md_only_v1/RESULTS.md). It evaluates eight MD targets: six degraded binders including EPHA2, and two non-degraded binders, ABL1/AXL. All targets have MD. The original 04h MD-only model was run for300epochs with seed42, no tuning, and paired ELOB/true-VHL roles.

Both arms return AUROC0.0,75% accuracy,50% balanced accuracy and0% specificity: every held-out target was called degraded. This equals the all-degraded accuracy reference and does not distinguish the two classes in this run. All16 saved fold models and32 prediction-only inferences passed independent CPU replay; protocol/source/input hashes and held-out membership were verified. No new docking or MD was run.

Current labels and interpretation: [CURRENT_LABELS.md](CURRENT_LABELS.md). Latest paired scores: [corrected_md_only_v1/results/paired_predictions.tsv](corrected_md_only_v1/results/paired_predictions.tsv). The label-map JSON records provenance as of its creation; STATUS.json gives current training completion.

The previous16-target static-plus-MD dual encoder remains frozen in the top-level config/protocol.json, analysis/gnn/, REVIEWED_RESULTS.md and results/paired_20261111/. That historical run used EPHA2=0 and its metrics do not describe the corrected label set. The new run changes cohort, architecture and labels together; comparisons with the old run cannot isolate a label-only effect. All old structures, graphs, models, predictions, assays and reports remain preserved.
