# SJF8240 corrected-feature reanalysis

**Current label correction:** EPHA2 is degraded. For the MD comparison, only ABL1 and AXL are non-degraded. See [CURRENT_LABELS.md](CURRENT_LABELS.md) and [the current label map](config/current_md_labels.json). The completed results below used the previous EPHA2 = 0 label and have not been retrained with this correction.

Completed 11 November 2026: all 19 static targets, 174 MD cases and both 16-fold model arms are finished and independently checked. Read [REVIEWED_RESULTS.md](REVIEWED_RESULTS.md) and [the actual methods](docs/ACTUAL_GRAPH_METHODS_20261111.md). Both overall AUROCs are 0.7273; both MD-subgroup AUROCs are 0.0. The aggregate result does not establish a reliable degradation predictor. Original data and historical results remain preserved. `config/protocol.json` records the protocol frozen before fitting; `STATUS.json` records current completion.

Separate derivative of the original protac-md workflow, created at the user's explicit request. Existing docking structures and MD trajectories are reused; no new docking or MD.

Folders: `config/` frozen protocol and source lists; `experimental_data/` LF/PS input tables and version-explicit target exports; `analysis/static/` rebuilt pose graphs; `analysis/dynamic/` rebuilt MD graphs; `analysis/e2/` existing corrected E2 summaries; `analysis/gnn/` new paired LOPO results; `scripts/` extraction/model tools; `logs/` progress and controls; `docs/` what/how/why notes.

Main comparison: 16 identity-matched historical targets, original 04d raw-graph dual encoder, paired ELOB and true-VHL selections within the repaired graph pipeline, 400 epochs, seed42, no tuning. CDK17-labelled/CDK14 structure is preserved but excluded before scores; SIK3/RPS6KB2 remain prediction-only. Experimental RPE1 conditions: 500 nM, 12 h; PS heavy marks newly synthesized proteins. Historical labels remain exploratory; experimental versions and uncertainty are explicit.

The ELOB control is a NEW same-pipeline rerun, not reconstruction of missing old weights or the 0.850 claim. Main graph-only model does not consume E2 scalars. Corrected E2 summaries are kept separately.

Check `logs/` and completion receipts for actual progress. A launched process is not a completed result. Licensed software, private runtime and bulky inputs are excluded from Git.
