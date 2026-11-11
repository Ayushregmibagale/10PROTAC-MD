# Corrected-label MD-only SJF8240 run

User-authorized new comparison using only retained MD contact graphs on eight targets. EPHA2, DDR2, EPHB2, MAPK14_full, MET and RIPK2_full are degraded; ABL1 and AXL are not degraded. All eight targets have MD. This version uses the exact original 04h MD-only GNN class, its 300-epoch schedule and seed42 with fixed target-bound fold seeds. No static encoder, E2 scalars, tuning or new docking/MD. True VHL is primary; ELOB is a paired control.

The protocol, current user labels, graph checksums and executable sources are frozen before fitting. Each LOPO fold trains on the other seven targets, saves weights and replays inference. Two prediction-only targets remain separate. Independent CPU replay, metrics and per-target score exports follow both completed arms. Check STATUS.json and logs for actual completion; this file is not a result receipt. The old16-target dual encoder remains unchanged and is not directly comparable as a label-only correction.

Local inputs and weights are excluded from Git; input_manifest.json records their source paths and hashes. The original script snapshot under provenance is evidence only; use scripts/run_pipeline.py as the execution entrypoint.
