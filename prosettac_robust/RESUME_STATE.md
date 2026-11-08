# Resume state — retained-structure integration review

Read AGENTS.md, USER_REQUIREMENTS.md and INTEGRATION_REVIEW_RESULTS.md first. Earlier checkpoints and all material what/how/why notes remain in SJF8240_REMEDIATION_LOG.md and Git history.

## User scope

Use the exact structures, parameters and completed PRosettaC/MD/GNN results from the original experiment. Harden scripts and review rationale. Confirmed panel: 17 historical GNN evaluation targets plus prediction-only SIK3 and RPS6KB2. The user wants GNN, not separate conventional predictive models. No supplier structure replacement, wholesale study rerun, Ernie production or replacement biological GNN training is authorized. Helpful bounded diagnostics remain allowed; none should be restarted automatically.

Original /home/ayushregmibagale/protac-md, /mnt/d/PROTAC degradation Prediction and prior top-level work are read-only. Derivative is /home/ayushregmibagale/10PROTAC-MD/prosettac_robust; stage is /tmp/10protac_finish/prosettac_robust.

## Requested remaining review checks have executed

- All 174 production TPRs represented in saved MD graph metadata have expected protein identities for their ten MD targets. Compiled settings: nominal 20 ns, 2 fs, 100 ps output; 113 cases at 1.0 nm cutoffs, 61 at 1.6 nm. AXL is mixed. A TPR request does not prove actual trajectory duration/completion.
- Original VHL-named metric/static selectors actually select Elongin B. True VHL is present. All 141 core existing trajectories have corrected contact audits reproducing stored GROMACS distances within 0.000022 A. Do not repeat this completed batch.
- E2 full-frame sweep: all 141 core trajectories, 29,142 frames, all metric/trajectory time axes aligned. Full-current-trajectory saved E2 reproduction succeeds for 140; AXL/pose011's saved row instead matches its explicit deduplicated 0–20 ns window, resolving the numerical discrepancy. No saved table changed.
- Twelve cases have CA direct/PBC image sensitivity. The earlier 15-frame sample missed them. Virtual Gaussian-point MIC is a stress test, not an adopted physical correction. Explicit-window summaries cover 141 cases with 201 unique frames each: 129 clear at this diagnostic and 12 audit-only, never automatic training data.
- Nineteen static representatives: 18 expected target matches; CDK17-labelled representative matches CDK14/O94921. All five retained AF/domain/clean-protein source PDBs in that folder also match O94921. Do not replace/rename the old structure, graph or experimental label. All 19 old static VHL selectors identify Elongin B. This is a representative audit, not every pose.
- Current graph artifacts are expanded: static PT 46, dynamic PT 10; static NPZ only two prediction targets. Do not assume interchangeable exports or frozen original training inputs. Old scalar runner filters artifacts but not unknown labels; hardened runner already uses explicit reviewed labels/cohort and input hashes.
- Seventeen evaluation rows belong to the E2-scalar GNN. Two extension prediction rows use the schema declared by the separate dual-encoder predictor without E2 scalars. Retained prediction log and current table disagree about MD availability and rounded scores; neither as-run lineage nor independent biological validation is established.
- Fractional graph-index/label coercion fixed in packing; packet target/count/source-drift guards added. Explicit E2 window/duplicate/image guard and missing-state handling added. No model architecture or old prediction changed.

## Current verification

141 analysis tests passed. Docking overlay: 137 passed, four PyG-only tests skipped. software_checkpoint_v10.json contains logs and all active script hashes. Installed checks run during guarded save; inspect newest installation receipt/log and actual Git history rather than assuming a partial save completed.

E2 geometry reconstruction preserves root images, uses periodic bonded spanning-tree reconstruction and checks all bonds with original fallback. Four real coordinate comparisons max 0.000016 A; six original complete cases /1,206 rows max numeric descriptor difference below 0.000027. Original v1 and ineffective guarded v2 audit attempts were stopped and retained; v3 completed. Do not restart them. Exact source snapshots and termination records are preserved.

## Main artifacts

- INTEGRATION_REVIEW_RESULTS.md: completed checks, actual findings and remaining scientific limits.
- HISTORICAL_METHODS_AND_RATIONALE.md: actual old algorithm/protocol/metric descriptions and script changes.
- existing_result_audit/e2_full_frames_v3: full per-frame audit, source snapshots, AXL window resolution.
- existing_result_audit/e2_explicit_window_all_core_v1: 141 explicit-window summaries, including 12 audit-only cases.
- existing_result_audit/integration_19_v2: target/config/artifact reconciliation and CDK17 prepared-source identity supplement. Canonical O94921 JSON and representative proof are retained in integration_19_v1.
- existing_result_audit/integration_review_v1: test/execution logs and review_completion.json.
- workflow/scripts/02_analysis/09_reanalyze_existing_contacts.py: retained core-layout contact analysis, no new production receipts or supplier replacement required.
- workflow/scripts/02_analysis/10_summarize_existing_e2.py: hash-bound E2 frame packet summaries; explicit time/duplicate/image policies.

## Scientific decisions/evidence still outstanding

Unconditional biological/predictive sign-off is withheld. Code cannot resolve the CDK17/CDK14 experimental-target correspondence, assay endpoint/label missingness and discordance, absent exact historical fold weights/input snapshots, or independent generalization/calibration. Preserve the original historical labels and findings; do not silently exclude cases to improve results or relabel simulations as supplier-equivalent. No further original-study computation is inferred from these limitations.

The supplier connectivity difference is retained provenance, not a prerequisite to continue the requested retained-structure review. MAPK14's scalar paired/unpaired estimator discrepancy and investigation patch remain documented; do not execute original scripts that upsert old tables.

## Saving / runtime

Generate the guarded plan with python3 /tmp/update_isolated_prosettac.py; wait for completion, then apply with approved escalation. Hash drift checks/backups must remain enabled. Keep edits under the isolated derivative. Git repo Ayushregmibagale/10PROTAC-MD master; current durable checkpoint is determined by git log/status. User requested frequent commits/pushes.

Compatible historical TPR reader: /home/ayushregmibagale/.conda/envs/gromacs/bin/gmx (2025.4). The 2024.3 PLUMED binary cannot read the historical 2025.4 TPR format. Private Python 3.10 docking overlay and licensed binaries remain excluded from Git. Analysis environment has native PyG; docking overlay does not.

No project simulation or analysis batch remains running after the completed review. Checkpoint installation/tests may run during saving. Do not imply work continues after an inactive session.

## 2026-11-08 completed supplement

**2026-11-08 update:** completed the full 141-pose coherent E2 imaging sensitivity and retained stage/checkpoint audits. All 174 production checkpoints reach 20 ns. All 22,753 true-heavy bridging frames agree across three imaging conventions; 16 differences are nonbridging. One EM force-tolerance exception is disclosed. Analysis tests: 139 passed; docking: 135 passed / four PyG-only skips. See [COHERENT_E2_REANALYSIS.md](COHERENT_E2_REANALYSIS.md), [METHODS_DRAFT_RECONCILIATION.md](METHODS_DRAFT_RECONCILIATION.md), [RETAINED_GNN_MODEL_PROVENANCE.md](RETAINED_GNN_MODEL_PROVENANCE.md) and [PROTEIN_PANEL_19.md](PROTEIN_PANEL_19.md). Original results unchanged; unconditional biological/predictive sign-off remains withheld. Earlier dated counts below are historical checkpoints.


All current numerical batches finished. New audit roots: ligand_connected_e2_cohort_v1, coherent_e2_publication_v1, retained_stage_evidence_v1, retained_checkpoint_times_v1, retained_gnn_weights_v1. Do not repeat these completed sweeps. Read remediation section 38. Main historical outputs remain unchanged; corrected imaging has no automatic training release. Check actual git/installation receipts for durable-save state.

## Latest reporting/assay supplement

**Latest 2026-11-08 reporting check:** 141 analysis tests pass; docking 137 passed / four PyG-only skips. Subgroup AUROC now evaluates independently and records undefined single-class results. Retained E2 scores yield overall AUROC 0.800, MD-only 0.600; no-MD AUROC is undefined. The current dual-GNN TSV does not reconstruct the older 0.850 baseline. Accession checks confirm the CDK17 abundance row cannot be assigned to the retained CDK14 structure. See [GNN_REPORTING_AND_ASSAY_RECONCILIATION.md](GNN_REPORTING_AND_ASSAY_RECONCILIATION.md). No original data, labels or predictions changed.

Completed gnn_reporting_v1 and assay_accession_v1; do not repeat. Read remediation section 39. No new biological model was fitted.

## User experimental evidence — latest

**Latest experimental evidence:** user-supplied RPE1 LF/PS data have been inspected (500 nM, 12 h; heavy=newly synthesized pool). LF `_new` provides matched results for all19 targets; five historical positive targets show reduced abundance in that saved version. Multiple result versions and weak/ambiguous evidence remain explicit; no label release or model refit. See [RPE1_EXPERIMENTAL_EVIDENCE_REVIEW.md](RPE1_EXPERIMENTAL_EVIDENCE_REVIEW.md).

Original sources /home/ayushregmibagale/SJF8240_LF_and_PS/E250916_LF and E260106_PS are read-only. Both audits completed; do not automatically repeat raw multi-GB reports or run original R scripts. Conditions are confirmed by user. Detailed provenance and remaining version/peptide/identity issues in remediation section40. GNN results unchanged.
