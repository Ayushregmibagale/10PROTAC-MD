# Scientific review decision — confirmed 19-protein SJF8240 project

**2026-11-08 update:** completed the full 141-pose coherent E2 imaging sensitivity and retained stage/checkpoint audits. All 174 production checkpoints reach 20 ns. All 22,753 true-heavy bridging frames agree across three imaging conventions; 16 differences are nonbridging. One EM force-tolerance exception is disclosed. Analysis tests: 139 passed; docking: 135 passed / four PyG-only skips. See [COHERENT_E2_REANALYSIS.md](COHERENT_E2_REANALYSIS.md), [METHODS_DRAFT_RECONCILIATION.md](METHODS_DRAFT_RECONCILIATION.md), [RETAINED_GNN_MODEL_PROVENANCE.md](RETAINED_GNN_MODEL_PROVENANCE.md) and [PROTEIN_PANEL_19.md](PROTEIN_PANEL_19.md). Original results unchanged; unconditional biological/predictive sign-off remains withheld. Earlier dated counts below are historical checkpoints.


**Latest completed review:** [INTEGRATION_REVIEW_RESULTS.md](INTEGRATION_REVIEW_RESULTS.md) supersedes the earlier limited E2/static coverage below. Full core E2 timestamp/imaging audit and 19-target reconciliation have run. All 141 time axes align; AXL's saved row matches its old window; 12 image-sensitive cases are flagged; the retained CDK17-labelled source is CDK14. Software passes 136 analysis tests. Scientific release remains withheld; no original data or labels changed.

**Subsequent user scope clarification:** retain the exact old structures and focus on code/rationale. Supplier equivalence below limits chemical/biological claims; it does not block the retained-structure methods review or require replacement simulations. Current implementation/methods findings are in [HISTORICAL_METHODS_AND_RATIONALE.md](HISTORICAL_METHODS_AND_RATIONALE.md). This report's withheld decision concerns unconditional scientific claims, not authorization to keep the historical structures.

Decision: **WITHHOLD unconditional scientific sign-off.** This is an evidence review, not a certification of publication readiness or journal acceptance. Software checks support specific behavior; they cannot repair unverified chemical identity, experimental endpoints or missing historical execution provenance. The user confirmed the scope as the 17 historical evaluation targets plus SIK3 and RPS6KB2. No conventional predictive model was added.

## Scope and evidence

The historical GNN evaluation contains 17 target-level predictions, not 19 independently validated predictions. SIK3 and RPS6KB2 are prediction-only targets; their saved NON calls are not experimental negative labels. The historical cohort has five recorded positive and twelve recorded negative labels; assay review has not released those labels for a new scientific evaluation. The 174 MD graph pose records span ten targets. The completed corrected contact audit covers 141 poses from eight core targets. Targets lacking completed MD must be reported explicitly; they must not be represented as having undergone an equivalent MD protocol.

See `existing_result_audit/scientific_signoff_v1/coverage_19_targets.json` for the per-target roster, manifest stages, graph frame counts and assay status. Counts refer to retained artifacts, not independently reconstructed original-run receipts. Earlier D-drive documentation describes a smaller, older panel and contains historical preliminary interpretations; it is not evidence that all 19 targets completed the same protocol.

## Decision by claim

| Claim / component | Decision | Evidence and remaining requirement |
|---|---|---|
| Isolated software regression | Pass within tested environment | 125 analysis-environment tests, including GNN input/cohort guards and synthetic fold-weight replay. Installed log: `validation_controls/regression/installed_20261107T010353Z.txt`. Earlier docking log has 121 passed/3 GNN-only skips; do not relabel that as a current 125-test docking result. Tests do not cover every possible scientific input. |
| PRosettaC preparation/linker implementation | Bounded implementation checks pass; scientific inference restricted | Paired fixed-head controls and GSK215 reference control retained. The changed linker method did not generate the historical study. Finite search failure is not nondegradation; a favorable docking score is not thermodynamic stability or biological validation. Target/construct/pocket and pose-selection provenance need to match each reported result. |
| MD handoff implementation | Bounded checks pass | Separate minimization/equilibration/short diagnostic controls pass. They are not independent production replicates and do not validate historical convergence. The selected torsion-method control completed minimization only. |
| Actual simulated compound equals supplier SJF8240 | Not established; contradictory connectivity evidence | Historical Mol2 key MCJFROJTDORRBB-ZWUDBOSPSA-N versus manufacturer SCUZPTNDZHEHLB-ZWUDBOSPSA-N. Parameter and actual-production-TPR audits bind historical connectivity to retained runs. Same formula/atom count is insufficient. Preserve raw runs and describe their actual identity; never relabel a TPR by editing an input. A claim specifically about supplier SJF8240 needs compound equivalence evidence or an explicitly different study scope. |
| Old VHL/POI bridging descriptor measures true VHL | Fails for the 141 core cases audited | Old VHL index selects Elongin B. True VHL is present in chain C. Corrected sequence-identified contact measurements reuse existing trajectories. Do not discard simulations solely because an analysis selection was mislabelled. |
| Corrected historical contact calculation | Pass for stated arithmetic/coverage | All 141 core existing-trajectory GROMACS distance reproduction checks pass within 0.000022 A. Strict <5 A, PBC and heavy/all-atom variants are explicit. The chemical scope remains the historical simulated graph. |
| Contact occupancy proves stable binding/degradation | Unsupported | Contact is a geometric proxy. 4/5/6 A sensitivity and four equal 5 ns block summaries are saved. Poses/frames are correlated and selected, not independent biological replicates. No equilibrium, free energy, convergence or calibrated degradation claim follows. |
| E2 reference frame and catalytic interpretation | Partial evidence only | Static ELOC anchors match across 141 starts. Bounded dynamic probe: 15 frames/three cases show no direct/PBC protein-distance gap, proper rotations, max ELOC fit RMSD 2.1213 A. This is not an exhaustive coherent-imaging or timestamp audit. Gaussian points are an assumed positional prior; 13 identical arrays represent 51 unique positions, not 663 independent samples. Reachability is not catalytic probability. |
| Historical GNN inputs and scores reproducible exactly | Not established | Current files and declared builder lineage are audited, but exact frozen historical training inputs and per-fold trained weights were not located in the inspected `analysis/gnn` artifacts. Existing graph packets can change after evaluation. One MAPK14 scalar differs from the declared paired estimator; a retained patch supplies a plausible alternative but not an as-run receipt. Preserve the mismatch and predictions. |
| New GNN implementation can export/replay folds | Pass in bounded synthetic tests | Tensor-only fold weights, normalization reference rows, held-out scores, seeds and hashes are saved and tested. Synthetic replay does not reproduce or validate the historical biological study. Corrected graph definitions require a separate reviewed evaluation; old metrics cannot be carried forward. |
| Experimental labels valid for the intended endpoint | Withheld pending endpoint/source resolution | Primary workbook has rows for 11/17 historical targets. Six missing rows are not negatives. MAP4K5 differs for this endpoint; ABL1 has proteomics/western discordance. Compound identity, cell line, time, dose, endpoint and threshold must be fixed independently of predictions. Existing labels remain unchanged. |
| GNN generalizes / scores are calibrated degradation probabilities | Unsupported by current audit | 17 small, selected target-level held-out results are exploratory. Feature availability is associated with labels. Saved scores do not establish incremental graph value, calibration, robustness across seeds or external generalization. No conventional model is required to recognize these limitations. GNN-only branch/feature ablations, reviewed labels, fixed evaluation and independent evidence are needed for stronger claims. |
| All 19 proteins fully scientifically validated | Withheld | Two are prediction-only; nine historical targets lack an MD graph in the retained evaluation. Partial coverage must remain visible. The separate evidence for prediction-target compiled protein identities does not supply experimental labels or a corrected full trajectory/contact evaluation. |

## Descriptive contact sensitivity (existing simulations only)

The new audit uses true-VHL/POI simultaneous heavy-atom minimum distances, strict cutoff, the same 201 stored frames per pose over 0–20 ns, and equal weighting of poses. Four sensitivity blocks are (0,5], (5,10], (10,15], (15,20] ns with exactly 50 frames each; frame zero enters the full-window mean only. This is not an equilibration-discard choice. All choices were specified without model fitting or performance optimization.

| Target | 4 A occupancy % | 5 A occupancy % | 6 A occupancy % | Range of 5 A block means, percentage points |
|---|---:|---:|---:|---:|
| ABL1 | 73.49 | 79.81 | 86.13 | 4.00 |
| AXL | 81.67 | 86.15 | 89.25 | 2.11 |
| DDR2 | 71.01 | 74.07 | 78.52 | 7.78 |
| EPHA2 | 66.56 | 70.23 | 73.49 | 18.67 |
| EPHB2 | 75.57 | 77.90 | 81.41 | 8.74 |
| MAPK14_full | 89.97 | 94.50 | 96.49 | 4.67 |
| MET | 87.13 | 89.43 | 90.24 | 7.13 |
| RIPK2_full | 67.34 | 70.76 | 72.37 | 0.82 |

These are descriptive values, not confidence intervals. Leave-one-pose-out ranges and per-pose blocks are also retained. A small block range does not prove convergence; a large one identifies sensitivity requiring careful interpretation. The original MD GNN metadata reports 21 sampled frames for every one of its 174 poses; that sparse feature sampling is distinct from the 201-frame corrected contact audit.

## What can be reported now

A provenance-preserving retrospective methods audit: the historical simulated connectivity, protein identities and mislabelled analysis roles; numerical reproduction of old distance outputs; explicitly versioned corrected contact descriptors; descriptive cutoff/window sensitivity; and bounded software/control results. Historical GNN predictions may be shown as historical exploratory outputs with their original labels and limitations, not as newly validated supplier-SJF8240 performance.

## Requirements for unconditional study sign-off

1. Resolve compound identity against the intended experimental molecule. Do not silently substitute the supplier structure into the description of completed runs. Existing runs remain retained; any decision about additional computation depends on the chosen claim and is not authorized by this report.
2. Freeze assay-specific ground truth, compound correspondence and endpoint. Resolve missing/conflicting evidence independently of the GNN. Prediction-only targets remain unlabelled until experimental evidence exists.
3. Complete target/structure/construct/pose lineage and frame correspondence for the artifacts actually used in each reported feature. Keep corrected role/feature definitions separate from old evaluated inputs.
4. Establish the sampling/replicate evidence required by the actual claims. Descriptive occupancies require honest conditional scope; equilibrium/free-energy/convergence claims need additional evidence. Do not create independent samples by counting frames or selected poses.
5. For validated GNN performance, use an immutable assay-reviewed cohort/input release, frozen protocol and retained fold weights. Report GNN-specific ablations, uncertainty and independent validation appropriate to the claim. No conventional ML model or full MD rerun is automatically required or authorized.
6. Obtain a scientific author review of the complete methods, claims and limitations. This automated review does not substitute for author responsibility or external peer review.

The sign-off request does not override the instruction to preserve all originals or authorize a wholesale rerun. No Ernie campaign or new biological-model training was launched for this review.

## Primary sources

Manufacturer structure reference: [Tocris/Bio-Techne SJF8240 #7266](https://www.bio-techne.com/p/small-molecules-peptides/sjf-8240_7266), including the listed InChIKey and SMILES, checked during this review. This is a product reference, not proof of the identity in a historical simulation.

Assay source under review: [Bondeson et al., Lessons in PROTAC Design](https://pubmed.ncbi.nlm.nih.gov/29129718/) and its [supplementary workbook](https://ars.els-cdn.com/content/image/1-s2.0-S2451945617303549-mmc4.xlsx). Product literature also cites a separate [Burslem et al. RTK study](https://pubmed.ncbi.nlm.nih.gov/29129716/). Compound numbering across publications must be matched explicitly; related foretinib PROTACs must not be assumed identical from a shared warhead alone.
