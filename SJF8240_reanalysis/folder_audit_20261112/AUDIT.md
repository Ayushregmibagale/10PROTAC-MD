# SJF8240 folder-wide audit and reuse hardening — 12 November 2026

## Outcome

The folder was not uniformly robust. Critical regressions were present in the newer geometric/E2 exploratory scripts. Twelve unsafe exploratory entry points are retired by explicit guards; their exact original source and all previous outputs are preserved. Completed historical studies remain frozen and reproducible, not automatically reusable for another compound. The new public entry point is `SJF8240_reanalysis/RUN_ANALYSIS.py`, backed by `folder_audit_20261112/robust_analysis_v1`.

The reusable contact extractor reproduced all141 existing SJF8240 MD cases, 28,341 frames and16 ELOB/VHL target graphs. Fresh CPU inference replayed all148 reviewed model weights across the historical dual encoder, corrected v1/v2, weighted/E2 and mechanistic experiments. Maximum absolute prediction error was1.1920928955078125e-7. No historical score, label, structure, trajectory or weight was changed. No scientific model selection, new docking or MD was performed.

## Scope and evidence

The initial recursive inventory contains125 project-owned Python sources (11,573 lines):72 operational/library files,13 tests,3 one-off recovery scripts and37 archived/snapshot files. All parsed successfully. Exact-duplicate source groups were recorded; shared libraries and model families were reviewed together. Another657 installed runtime dependency Python files were inventoried and parsed, but are not claimed to have received a semantic third-party software audit. The125 count includes ignored historical backup/log source files, unlike Git's smaller tracked-source list. Four further E2 scripts appeared during the audit and were separately reviewed and captured in late_source_inventory.json; total reviewed pre-existing/late project sources:129.

Review covered role selection, source/frame/pose membership, periodic contacts, residue correspondence, graph tensor semantics, label evidence, unknown targets, grouped splits, fold-local loss/scaling, metric definitions, checkpoint/restart behavior, input/source hashes, runtime boundaries, replay and reporting. Static flags are review prompts, not automatic proof of defects. Tests and source review do not constitute a mathematical proof of correctness.

Six active frozen source manifests, containing80 file bindings, matched their recorded hashes. The four exact MDEncoder variants and the temporal architecture were preserved. Historical EPHA2=0 runs remain historical; the current confirmed-user map has EPHA2=1. SIK3/RPS6KB2 remain experimentally unknown outside the explicitly archived v2 assumption experiment.

## Findings and disposition

| ID | Severity | Finding | Resolution |
|---|---|---|---|
| F01 | Critical | geometric_gnn/01_extract_ca_positions.py and e2_geometry/02_pose_ensemble_geometry.py identify chain A residues1–103 as VHL. In the independently audited retained systems, this selects ELOB. | Retire entry points; preserve original sources. New extractor verifies exact model sequence and independent canonical UniProt sequence identity before accepting VHL/ELOB/POI roles. |
| F02 | High | geometric_gnn/04_build_per_pose_and_bridge_graphs.py guesses protein blocks from GRO residue jumps and computes direct nonperiodic contact distances. | Retire builder; require explicit sequence-verified selections and periodic distances. Full141-case replay matched reviewed extraction. |
| F03 | High | Same builder multiplies marginal contact frequencies and fabricates an edge for a graph with no qualifying contacts. This does not measure simultaneous bridging. | Replace with per-frame joint observations. Empty graphs remain empty. Anticorrelated-contact and empty-edge regression tests pass. |
| F04 | High | Exploratory trainers use pos_weight3 with degraded=1; one comment incorrectly says this weights positives less. It actually increases positive weighting. | Retire legacy entry points. New runner explicitly selects unweighted or training-fold-balanced loss. Tests verify class totals and held-out exclusion. Historical fitted losses/results remain unchanged. |
| F05 | High | Exploratory loaders skip missing targets/poses, broadly suppress warnings, and omit fold weights or source-bound independent replay. Some comparison values are hard-coded or printed as an unbound E2-GNN~0.85 reference. | Retire as reviewed workflows. Preserve outputs as exploratory/unverified. New runner requires every declared input and saves fold states, source/data bindings and replay receipts. Printed reference numbers are not audited evidence. |
| F06 | High | SJF-specific labels, absolute paths, chain IDs, ligand counts, cohort sizes and frame grids appear across historical pipelines. Copying these across10 compounds can silently change scope. | Leave frozen studies intact. New manifests supply compound identity, ligand, protein roles, exact time grid and label evidence; batch preflight rejects missing compounds. No SJF cohort default for a new compound. |
| F07 | High | Some exploratory dual/combined models insert synthetic zero graphs for missing modalities and mix non-MD negative targets with the labelled MD cohort. Per-pose augmentation changes optimization weights and sample count, not the independent protein count. | New supported runners require declared complete role inputs and one target-level packet; no hidden modality imputation or pseudo-replication. Joint multimodal/per-pose exploratory models are not approved by this audit. |
| F08 | High | e2_geometry/02_pose_ensemble_geometry.py uses pose-level observations in a class Mann–Whitney test, ignoring repeated poses of the same protein. | Retire inference path; retain per-target descriptive summaries and explicit independent-protein counts. No replacement significance claim. |
| F09 | High | Exploratory E2 geometry lacks a source-bound demonstration of calibrated E2/reference-frame alignment; N/C lobes are guessed by a40/60 sequence split. | Retire entry points. Only previously audited coherent-image lysine descriptors remain supported for SJF8240. Other compounds require explicit calibration and annotation review, not copied E2 assumptions. |
| F10 | Medium | mechanistic_v3/scripts/test_contracts.py uses a Windows staging sibling name and fails from the canonical project location. Its original four-test pass was in the staging workspace. | Preserve original frozen source and the observed failure log. Provide a hash-bound portable test amendment with the correct canonical fixture root; all4 portable tests pass. No model/data change. |
| F11 | Medium | mechanistic_v3 pair trainers lack a run lock and atomic checkpoint replacement; some resumability checks are weaker than older trainers. Several old launchers write a fixed canonical status and use exclusive log filenames. | Frozen research executables remain restricted. New runner locks output, atomically saves state, archives incomplete weights, validates source/input/fold/runtime on resume and uses caller-owned output. |
| F12 | Medium | Some scripts use assert for production validation, and old shared scalar helpers have zero/missing-value fallbacks. They are unsuitable as a general unknown-feature policy. | New reusable contracts raise explicit exceptions and reject nonfinite/missing graph features. Unsupported E2 scalar generation is not silently substituted. |
| F14 | High | Four late E2 scripts add unbound site lists/residue-number matching; CRL2 alignment swaps ELOB/ELOC/VHL and truncates arrays without residue correspondence; kinetics uses legacy ELOB distances as VHL and uncensored l×dt durations; activation-loop RMSD includes rigid motion and does not implement its advertised DFG classification. They also retain missing-input skips/pseudoreplication. | Retire all4 CLI paths with original snapshots. Preserve all existing TSVs as unverified exploratory outputs; do not present site annotations, clash, kinetics or loop summaries as validated mechanism. New comparable methods require verified annotation provenance/residue mapping, sequence-matched alignment and time/censoring contracts. |
| F13 | Scientific limit | Pooled LOPO ranking and matched-pair ranking are different estimators; pairs reuse the same2 negative targets. The0.750 VHL result has0% specificity. | Metrics explicitly name the estimator; report balanced accuracy/sensitivity/specificity and independent target counts. No score inversion, threshold search or biological certification. |

## Runtime validation

- Original parent suite:15 passed.
- Corrected MD-only v1:4 passed.
- Corrected hypothetical v2:6 passed.
- Balanced-loss version:8 passed.
- Balanced-loss + E2 version:11 passed.
- Original mechanistic canonical suite:3 passed,1 failed due to the documented path bug. This failure is retained, not suppressed.
- Portable mechanistic replacement:4 passed.
- New reusable workflow:28 passed, including short save/resume/CPU plumbing tests for both architectures, tamper rejection, unknown exclusion, sequence-role swap rejection, joint contacts, periodic translation and graph validation.

The current documented verification route therefore has76 passing tests; this does not mean the unamended historical canonical test passes. The one-epoch synthetic GNN tests are software checks, not new degradation results. Historical study seeds/epochs/weights were never tuned or replaced.

Fresh read-only CPU replay:32 historical dual +16 corrected v1 +20 hypothetical v2 +16 weighted +16 weighted/E2 +48 mechanistic =148 saved models. Replay JSONs retain every fold, checkpoint hash and error. The new extractor replay used two target-disjoint shards covering all141 cases, with checks against per-pose arrays and final graph tensors. All16 topologies match exactly; numeric tolerances and actual errors are recorded in real_extraction_verification.json.

## What to use for the10 PROTACs

Use the public launcher and per-compound manifests described in robust_analysis_v1/README.md. Validate every expected compound before fitting. The contact extractor and GNN runner are configurable and have been exercised end to end with SJF8240. The original and enriched GNN architectures are supported as separately declared candidates. Default fitting is never selected by whichever previously scored best.

For each compound, provide its own source/chemistry provenance, MD topology and trajectory hashes, ligand selection/count, independently identified proteins, expected frame grid, binding/degradation label evidence, unique protein groups and graph manifests. Known binders with unknown degradation status remain prediction-only. The two SJF examples deliberately validate just SJF8240; they do not invent labels or runs for the other9 compounds.

The generic raw-trajectory extractor currently produces original-MD joint-contact graphs for VHL/ELOB systems. Existing SJF8240 lysine-aware graphs can be used through the guarded runner. Generating equivalent enriched graphs for another compound still requires a reviewed E2 calibration, imaging/lysine correspondence, feature scales and four-window configuration. The frozen SJF-specific enriched extractor is not certified as a universal generator. CRBN/other ligases and joint cross-compound training are outside the implemented configuration; they must not be relabelled VHL. If repeated proteins are pooled across compounds later, hold out their biological identity groups across compounds.

The audit makes the reuse path stricter and reproducible. It does not promise similar results for all compounds, confirm experimental degradation, or establish productive E2–ubiquitin catalysis.

## Files

script_inventory.json and duplicate_groups.json contain the initial125-source project inventory plus installed runtime inventory. review_disposition.tsv provides per-file disposition. frozen_source_verification.json records immutable source checks. tests_*.txt retains every initial test run including the failure; portable_tests.txt/new_tests.txt record the passing replacements. replay_*.json contains148 fresh inference receipts. real_extraction_verification.json records the full trajectory/graph replay. Original retired sources are under archived_exploratory_sources/. quarantine_receipt.json binds all12 originals and guarded derivatives; late_source_inventory.json records the4 sources added while auditing. reviewed_entrypoints.json and audit_folder.py detect new or changed project code. AUDIT.docx is the readable Word copy.
