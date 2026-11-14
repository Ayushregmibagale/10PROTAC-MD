# Ten-target mechanistic follow-up

This study audits existing geometry, creates matching unknown-target GNN inputs and predicts without fitting. EPHA2 stays a degrader. SIK3/RPS6KB2 stay unknown. Original dynamic graphs, trajectories, labels and models are preserved.

Use `/home/ayushregmibagale/miniconda3/bin/python3`. Run numbered scripts in order:

1. `01_extract_matching_inputs.py --workers 4`: all33 unknown poses plus the first lexicographic pose for each of8 known targets. Reconstruct whole components, identify coherent images, check independent periodic contacts and timestamps, compute exact archived base/native descriptor definitions. Source hashes and sequence-role evidence bound to TPR.
2. `02_build_and_validate_graphs.py`: compare8 re-extracted controls with archived inputs, then build4 matching unknown temporal34feature graphs. Preserve29base features, identical edges and zero-padded control.
3. `03_score_unknowns_frozen_models.py`: use120 previously saved checkpoints; predict only, replay each unknown from fresh model, reproduce240 original held-out reference scores. Save score ranges and votes by variant/seed.
4. `04_collect_ensemble_reach50.py`: recover continuous per-lysine E2-ensemble reach50 for141 core trajectories. Reproduce all stored minimum distances and trueVHL bridge masks. Does not equate nearest-point presence with ensemble fraction.
5. `05_mechanistic_audit.py`: audit all174 retained poses,10targets,4time windows,all modeled lysines and6 reference/anchor combinations. Preserve undefined conditional deltas when no bridging/nonbridging frames exist. No per-frame/per-lysine pseudo-replication p-values.
6. `06_independent_geometry_checks.py`: independent unwrap/SVD/brute-force distance,clash,exposure checks at frames0,100,200 for first poses ofABL1/AXL/SIK3/RPS6KB2. Verify all8 unknown graph-window schema contracts.
7. `06b_map_residue_numbering.py`: sequence-align each modeled fragment to the stored canonical reference; flag construct differences and ambiguous mappings, add verified canonical lysine IDs without changing models. Run after05, before07; no report artifact is required.
8. `07_write_report.py`: reproduce archived8target geometry scores, write target overview, figures and findings in Markdown/Word; require all validation receipts PASS.

Primary geometry gates, features, reference choice and model weights were fixed before examining unknown predictions. Alternatives are reported as sensitivity checks, never selected to improve label agreement. Original 0.8333 single-seed result is superseded in robustness interpretation by the existing five-seed study (control0.733,catalytic0.667 means); no GNN retraining occurs here.

Read FINDINGS.md/docx for evidence and limitations. Model scores are not calibrated experimental probabilities. Directory names ending _full do not prove full-length protein models. Scientific catalytic validation remains incomplete. Historical molecule provenance remains separate from supplier-SJF8240 equivalence.
