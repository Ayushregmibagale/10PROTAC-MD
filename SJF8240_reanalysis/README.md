# SJF8240 corrected-feature reanalysis

Separate derivative of the original protac-md workflow, created at the user's explicit request. Existing docking structures and MD trajectories are reused; no new docking or MD.

Folders: `config/` frozen protocol and source lists; `experimental_data/` LF/PS input tables and version-explicit target exports; `analysis/static/` rebuilt pose graphs; `analysis/dynamic/` rebuilt MD graphs; `analysis/e2/` existing corrected E2 summaries; `analysis/gnn/` new paired LOPO results; `scripts/` extraction/model tools; `logs/` progress and controls; `docs/` what/how/why notes.

Main comparison: 16 identity-matched historical targets, original 04d raw-graph dual encoder, paired ELOB and true-VHL selections within the repaired graph pipeline, 400 epochs, seed42, no tuning. CDK17-labelled/CDK14 structure is preserved but excluded before scores; SIK3/RPS6KB2 remain prediction-only. Experimental RPE1 conditions: 500 nM, 12 h; PS heavy marks newly synthesized proteins. Historical labels remain exploratory; experimental versions and uncertainty are explicit.

The ELOB control is a NEW same-pipeline rerun, not reconstruction of missing old weights or the 0.850 claim. Main graph-only model does not consume E2 scalars. Corrected E2 summaries are kept separately.

Check `logs/` and completion receipts for actual progress. A launched process is not a completed result. Licensed software, private runtime and bulky inputs are excluded from Git.
