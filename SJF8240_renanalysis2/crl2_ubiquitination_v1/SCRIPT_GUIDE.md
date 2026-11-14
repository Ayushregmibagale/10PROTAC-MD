# Ordered script guide

Use `/home/ayushregmibagale/miniconda3/bin/python3`. Dependencies: NumPy,pandas,SciPy,MDAnalysis,Biopython,scikit-learn,matplotlib,python-docx. Original paths are specific to this workstation. Use the numbered entry points and do not invoke Python with optimization flags. This is not a generic new-target production pipeline.

| Order | Entry point | Purpose |
|---|---|---|
| 00 | `00_fetch_references.py` | Verify downloaded immutable references; `--download` fetches missing files only and refuses changed remote content. |
| 01 | `01_prepare_scaffolds.py` | Construct observed closed and modeled open scaffold packets; explicitly model and check missing donor terminal coordinates. |
| 02 | `02_extract_md_compatibility.py --workers 4` | Process all141 existing cases and201 frames each. Uses source hashes,role audits,whole-component imaging and independent contact replay. `--limit 1` is a bounded diagnostic,not a completed cohort. |
| 03 | `03_validate_controls.py` | Analytic geometry/exposure controls and native MZ1 projection; scientific-control failure remains explicit. |
| 04 | `04_compare_targets.py` | Produce target rankings and operational metrics with equal-pose aggregation. Refuses incomplete141-case cohorts. All120 prespecified sensitivity settings are saved. |
| 05 | `05_independent_verification.py` | Replay six raw frames through independent reconstruction/alignment/distance implementations,check exposure and overlap,and verify AUROCs by pair counting. |
| 06 | `06_write_report.py` | Write findings Markdown/Word and target comparison PNG/PDF; requires independent verification. |
| 07 | `07_check_package.py` | Check immutable code/configuration,all141 packets,cohort counts and numerical result checksums. Does not grant scientific validation. |

All steps are already completed in this version. Regeneration of stage01 is unnecessary. Stage02 cached results require the same frozen source manifest; a method/reference/config change requires a fresh version and new results,not manually refreshing the old manifest. New target data require explicit cohort/protocol review. Never overwrite or relabel historical inputs. The eight labels are copied from the explicit user-confirmed mapping; SIK3/RPS6KB2 remain unknown and are excluded.

The extraction manifest binds the exact extraction code,protocol,copied helpers and scaffold packets. The software manifest separately binds later control/reporting/verification code. Raw TPR,XTC and reference PDB files are hashed against the prior reviewed receipts for all141 cases. The read-only XTC reader avoids offset sidecars in source directories. Cached extraction verifies source/result hashes; original input rehashes were completed for this run and must be repeated for a fresh input audit.

Primary score settings are frozen in config/protocol.json. Target-induced scaffold clashes count against the denominator. Missing atoms,modeled Gly76,mutant E2,missing sidechains,open-scaffold module transfer,static native-control rotamers and chemistry mismatch are documented in FINDINGS.md. Numeric PASS does not mean binder/degrader separation is validated. No automatic GNN-training or MD-production stage exists here.

Principal outputs: results/primary_target_scores.tsv,results/metrics.tsv,results/target_scores_all_variants.tsv,results/sensitivity_metrics.tsv,results/primary_frame_summary.tsv and results/poses/<target>/<pose>/packet.npz. Evaluation scripts never treat frames orlysines as independent labelled proteins.
