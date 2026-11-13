# SJF8240_renanalysis2 — ordered script guide

This is the clean, numbered entry point for the reviewed SJF8240 analyses. The original SJF8240_reanalysis folder and raw MD data remain in place. Supporting implementations are under _workflow/; normal use starts with the numbered scripts below. All new generated data go under outputs/. No retired exploratory script is included.

## Script order and purpose

| Script | What it does | Inputs | Outputs / when to use |
|---|---|---|---|
| 01_prepare_inputs.py | Verifies the original graph/label bindings and prepares the retained-trajectory configuration. | config/source_project.json; original reviewed manifests, labels and graphs. | config/trajectories.json and outputs/preparation.json. Start here. |
| 02_build_md_contact_graphs.py | Rebuilds original MD joint-contact graphs with sequence-verified VHL/ELOB roles and periodic distances. | config/trajectories.json; original TPR/XTC files. | outputs/md_contact_graphs/. Optional: skip when reusing existing graphs. |
| 03_extract_lysine_features.py | Measures per-lysine NZ solvent exposure/burial and E2 geometry at41 fixed frames per pose. | Existing SJF8240 trajectories and audited reference/contact evidence. | outputs/lysine_features/analysis/. Optional: use when rebuilding enriched GNN inputs. |
| 04_build_temporal_graphs.py | Combines the extracted lysine features with residue identity and contacts into29-feature graphs for4 chronological windows. | Completed03 feature packets. | outputs/lysine_features/inputs/temporal/ and a graph manifest. Run after03. |
| 05_prepare_gnn_study.py | Creates a study for original_md or lysine_temporal, explicitly selecting existing or rebuilt graphs. | Reviewed existing graphs, or completed02 /04 outputs. | config/studies/<model>.json with graph hashes and labels. Run once per desired model. |
| 06_validate_gnn_studies.py | Checks every prepared study: graph tensors/hashes, protein groups, labels, unknown targets and fold definitions. | config/studies/*.json. | config/batch.json and outputs/preflight.json. Complete before training. |
| 07_train_gnn.py | Fits leave-pair-out GNNs using the frozen study settings, saving every fold and its membership/scores. | Validated study; selected model, role and device. | outputs/gnn/<model>/<role>/. This starts real training. |
| 08_verify_gnn_cpu.py | Reloads saved weights in a CPU process and reproduces predictions. | Completed07 output. | cpu_verification.json beside each run. Use a CPU-compatible Python environment. |
| 09_compute_full_frame_e2.py | Recomputes per-lysine reach50 on all201 frames per retained pose, using4 legacy/corrected geometry/contact combinations. | Original MD trajectories, reference structures and audited contact evidence. | outputs/e2_full_frame/analysis/. Separate feature-validation branch; does not train a GNN. |
| 10_compare_e2_results.py | Aggregates the full-frame deltas by lysine/protein and computes exact protein-level statistics plus matched-pose controls. | Completed09 packets. | outputs/e2_full_frame/results/; the full report is collected by12. |
| 11_verify_e2_results.py | Freshly replays8 selected full trajectories with standard MDAnalysis unwrap, then rehashes all packets and source inputs. | Completed09 data and original source files. | outputs/e2_full_frame/validation/. Complete before exporting the E2 result. |
| 12_write_summary.py | Collects only analyses actually run in this new folder into a single readable summary. | New GNN completions/CPU receipts and/or new verified E2 results. | outputs/SUMMARY.md. If none ran, it explicitly says so. |
| 13_check_workflow.py | Checks packaged source hashes, copied-backend provenance, Python syntax and frozen E2 source manifests. | This folder's scripts/configuration. | outputs/workflow_check.json. Run at the end or whenever checking the installation. |

Every launcher accepts --help and --dry-run. A dry run displays its selected action without extracting trajectories or fitting a model. Do not use Python -O because the preserved temporal-feature builder contains assertion-based validation.

## Most common route: reuse existing GNN graphs

Run from this folder, using an environment containing the required scientific packages:

```
python 01_prepare_inputs.py
python 05_prepare_gnn_study.py --model original_md --graph-source existing
python 05_prepare_gnn_study.py --model lysine_temporal --graph-source existing
python 06_validate_gnn_studies.py
```

To start a new fitting run, choose the model and role explicitly. --role both runs VHL and ELOB sequentially. The examples default to300 epochs and seed42; the runner saves fold weights and does not silently change the threshold.

```
python 07_train_gnn.py --model original_md --role both --device cpu
python 08_verify_gnn_cpu.py --model original_md --role both
python 07_train_gnn.py --model lysine_temporal --role both --device cpu
python 08_verify_gnn_cpu.py --model lysine_temporal --role both
python 12_write_summary.py
python 13_check_workflow.py
```

For GPU fitting, use the existing compatible CUDA/PyG Python environment and --device cuda. Use the base CPU Python for08; its child process removes the GPU-only PYTHONPATH overlay. There is no automatic CUDA-to-CPU fallback. Training may take substantial time; preparing this folder has not started a fitting run.

## Rebuild inputs instead

For original contact graphs:01 →02 →05 --model original_md --graph-source rebuilt →06 →07 →08.

For enriched temporal graphs:01 →03 →04 →05 --model lysine_temporal --graph-source rebuilt →06 →07 →08.

Rebuilt studies use the newly generated graphs, not the old graphs. The retained141-case rebuild contains the8 labelled targets; unknown targets without rebuilt inputs are explicitly excluded, not assigned zero features or negative labels. Existing original-MD studies can retain SIK3/RPS6KB2 as prediction-only inputs.

The contact builder requires a fresh output directory and preserves failed partial extraction. The GNN runner can resume only when its study, input hashes, source and training runtime match. Reusing cached lysine/E2 packets first checks their bound external inputs. Final GNN completion should include CPU verification, not merely saved weights.

## Full-frame E2 check

This branch can run independently of fitting:

```
python 01_prepare_inputs.py
python 09_compute_full_frame_e2.py --workers 2
python 10_compare_e2_results.py
python 11_verify_e2_results.py
python 12_write_summary.py
python 13_check_workflow.py
```

03/04 reproduce the41-frame feature representation used by the enriched GNN.09/10 investigate the specific50-A reachability-change hypothesis on201 frames. These are different analyses: the0.750 historical GNN did not directly use per-lysine delta_reach_50A. The full-frame comparison is a feature check, not another GNN model.

## Labels, provenance and interpretation

Current user-labelled degraders: DDR2, EPHA2, EPHB2, MAPK14_full, MET, RIPK2_full. Non-degraders: ABL1 and AXL. SIK3_full/RPS6KB2_full have unknown experimental degradation labels. They are never made confirmed negatives by this workflow.

VHL is chain C and ELOB chain A in the retained audited SJF8240 compiled systems, but the contact builder verifies protein sequences instead of trusting those letters alone. Reference/imaging and ligand assumptions are SJF8240-specific. Another PROTAC requires its own configuration, labels, ligand and E2 calibration; this clean folder is not permission to reuse SJF labels or priors blindly.

The reviewed configurable GNN runner uses JSON-encoded pair-bound seeds, which differ from the historical mechanistic_v3 name-delimited seeds. A new fit is therefore a separately versioned run; it is not guaranteed to reproduce the old0.750 score. The copied previous reports under reference_results/ are explicitly historical, not outputs of new training. Matched-pair AUROC and pooled LOPO AUROC are different estimators. Frames/poses/lysines do not increase the independent labelled-protein count.

## Folder map

- Numbered scripts: the user-facing execution order.
- SCRIPT_GUIDE.docx / SCRIPT_GUIDE.md: this explanation.
- config/: study settings, source/data bindings and generated manifests.
- _workflow/: reviewed supporting implementations; normal users do not need to run these directly.
- outputs/: newly generated graphs, feature packets, weights, validation receipts and summary.
- reference_results/: snapshots/links for earlier findings; no bulk model or trajectory copies.

Keep the original source project and trajectories available: the new folder references those read-only assets. It is an organized script package, not a second copy of the MD dataset.
