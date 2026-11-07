# Retained-structure integration review: completed checks and decisions

The requested code, E2 frame and 19-target methods review has been executed. **This is not unconditional scientific sign-off.** The review found real data/provenance limitations; changing code cannot establish missing experimental truth or turn a different protein into CDK17. All original structures, simulations, labels and predictions remain unchanged. No MD or biological GNN training was repeated.

## Checks completed

| Check | Result | Evidence |
|---|---|---|
| MD-to-analysis source roles/protocols | All 174 retained production systems represented in MD graph metadata have identified expected proteins; compiled settings are recorded per case | `compiled_md_protocol_v1`, `protein_roles_core_complete`, `prediction_compiled_roles_v3` |
| Full E2 metric/trajectory correspondence | All 141 core cases aligned over all 29,142 available frames; no positional truncation accepted as alignment | `e2_full_frames_v3` |
| Retained E2 numerical reproduction | 140 cases reproduce their saved rows over the entire retained trajectory; the remaining AXL case reproduces under its original 0–20 ns window | `e2_full_frames_v3/AXL_pose011_window_resolution.json` |
| Explicit common-window E2 summaries | 141 cases, 201 unique frames each over 0–20 ns; exact duplicate handling recorded | `e2_explicit_window_all_core_v1` |
| Dynamic image sensitivity | 12 cases flagged; up to 104.88 A difference between direct and periodic protein-CA minimum distances | Same full-frame audit; flags retained rather than converted to negatives or silently dropped |
| Static target identity | 19 representatives examined: 18 match expected POI; the CDK17-labelled representative matches CDK14 instead | `integration_19_v2`, CDK17 source-identity supplement |
| Legacy static VHL selector | Elongin B in all 19 sampled representatives | `integration_19_v2` |
| GNN packet handoff | Fractional index/label coercion fixed; target/count checks and packet source-drift hashes added; valid packet path tested | `04_pack_static_graphs.py`, regression tests |
| GNN missingness/cohort boundaries | Historical/current artifact differences documented; future runner already excludes unknown labels and requires reviewed cohort and input release | `integration_19_v2`, `gnn_contract.py` |
| Software regression | 136 analysis tests pass; docking has 132 passed and four GNN/PyG-only skips | `software_checkpoint_v8.json` and logs |

These are scoped checks, not an exhaustive review of every historical docking pose or experimental assay. Successful source reproduction does not validate a feature's biological interpretation.

## E2 resolution and new controls

The complete sweep includes the extended AXL/pose011 trajectory, rather than assuming all postprocessed files end at 20 ns. Its current full series has 1,002 rows through 100 ns. Removing one identical duplicate and selecting the explicit 0–20 ns window gives 201 unique frames and reproduces the saved E2 row within retained precision. This explains the numerical mismatch without changing any saved value; it does not prove when the old code ran.

The earlier three-case/15-frame probe missed the 12 image-sensitive cases. The full-frame audit now supersedes that limited coverage. Flag counts are ABL1 4, AXL 1, DDR2 2, EPHA2 2, EPHB2 1, MAPK14_full 2, MET 0 and RIPK2_full 0. Protein-CA checks identify image sensitivity, not every possible heavy-atom interaction. Virtual Gaussian-point minimum-image distances change substantially too; they are a stress test, not a justified replacement for the physical E2 model.

`workflow/scripts/02_analysis/10_summarize_existing_e2.py` now requires explicit time bounds, duplicate policy and image policy. It validates matching metric/trajectory frame times and hashed frame/metric receipts. `--image-policy reject` refuses an image-sensitive summary. `audit-only` preserves it with `IMAGE_REVIEW_REQUIRED_AUDIT_ONLY`; it never issues a training release. Missing bridge or nonbridge states remain missing, not zero. Conditional differences are explicitly paired within a pose. The old state mask is named ligand/Elongin-B/POI proximity, not relabelled as true-VHL pocket engagement.

All 141 common-window summaries were generated. Of these, 129 have no CA image flag and 12 are retained as audit-only. A clear flag is not proof that the Gaussian prior or catalytic interpretation is valid; absence of a flag is not complete coherent-assembly certification.

The slow original reconstruction and guarded attempt were stopped only for this audit, with partial results/source snapshots retained. The final bonded spanning-tree implementation preserves fragment roots, uses periodic covalent-bond displacements and verifies every reconstructed bond; inconsistent cycles fall back to the original operation. Four real reference frames differ from original coordinates by at most 0.000016 A, and six completed original cases (1,206 frame rows) have maximum numeric descriptor difference below 0.000027. Original and branched/multiple-fragment synthetic checks pass. This is a numerical implementation comparison, not a new physical model or independent biological validation.

## CDK17 structure identity

The reviewed source folder documents AlphaFold accession O94921. That accession is [CDK14](https://rest.uniprot.org/uniprotkb/O94921.json); human [CDK17 is Q00537](https://www.uniprot.org/uniprotkb/Q00537/entry). The retained representative matches O94921 at 100% sequence identity and observed-sequence coverage. The five retained AF/domain/clean-protein preparation files were also checked and match O94921. The historical graph's 350 POI nodes are consistent with that construct length, but this does not reconstruct every original graph-building execution.

Do not rename CDK17's experimental label to CDK14's, swap the structure, silently drop the target to improve accuracy or claim that this is a validated CDK17 structural example. Keep it explicitly marked as a historical CDK17-labelled/CDK14-structure case. Any biological evaluation change requires an independently chosen cohort/protocol; none was run here.

## Nineteen targets are not nineteen equivalent validation results

The 17 historical rows belong to the declared E2-scalar GNN runner. The two prediction-only rows use the output schema declared by the separate dual-encoder GNN predictor, without its E2 scalar branch. They are not two extra held-out experimental validation results. Current static PT contains 46 graphs and dynamic PT contains 10; these are expanded retained artifacts, not frozen historical 17-fold inputs. The static NPZ sidecar contains only the two extension targets, so it is not an interchangeable export of all 46 graphs.

The retained `predict_kinome.log` reports no MD for the two extension predictions and scores different from the current table, whose MD flags are true. The log cannot establish the exact execution provenance of the current table. No replacement predictions were fitted. The old scalar runner's cohort comprehension filters artifacts but not unknown labels; rerunning that old script on today's expanded packet is unsafe. The hardened runner's explicit reviewed cohort and unknown-label guards remain in place.

Legacy directory categories are not experimental label evidence. In particular, CDK17 remains in a directory called `low_affinity_degraders` although the saved evaluation label is zero. No label was inferred from that directory name or from docking/MD failure.

## Scientific questions not solved by software

Unconditional predictive/biological sign-off remains withheld because:

- The CDK17 structure identity does not match its intended target.
- Twelve E2 cases require a scientifically justified coherent-assembly/image policy before claiming physical reachability; this review does not silently adopt virtual-point PBC as a cure.
- Assay-specific missing/conflicting labels, including six absent rows in the reviewed table, remain unresolved; prediction-only targets remain unlabelled.
- Historical fold weights and exact frozen inputs are unavailable in the inspected artifacts; current prediction-table/log and one scalar-estimator discrepancies remain provenance limits.
- Selected targets/poses and correlated frames do not establish independent generalization, calibration or production convergence. Neither more passing tests nor retaining the old molecule can supply that evidence.

The scripts and rationale can now be reviewed on the retained-structure basis with these findings visible. This report is not an instruction to rerun the experiment. Stronger claims need an explicit scientific decision or additional evidence, not a fabricated PASS.
