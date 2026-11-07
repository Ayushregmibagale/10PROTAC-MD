# SJF8240 publication work: current status

Current priority: preserve/reuse the existing SJF8240 structures and completed PRosettaC/MD/GNN results; harden scripts and audit scientific rationales. Helpful bounded diagnostics are authorized, but a wholesale experiment rerun is not. Read `USER_REQUIREMENTS.md` and `AGENTS.md`. Scientific claims remain under review.
The original D: project, original protac-md, original MD/GNN results, and the
previous top-level 10PROTAC-MD workflow are unchanged.

**Validation:** 124 analysis-environment tests pass. Docking-environment checks pass 121 with three GNN-only skips (PyG absent). Tests verify software behavior, not publication claims.

## Findings established from retained artifacts

- All 255 inspected SJF8240 ITP files have the historical named heavy-atom connectivity. A separate element-labeled graph isomorphism check matches that parameter graph to the historical Mol2 and rejects the Tocris graph.
- All 174 actual production TPR files linked to the saved MD-graph pose metadata contain that historical connectivity. The 200-row random-screen manifest records 174 ANALYSIS_DONE, 19 FAILED and seven EM_DONE; all 174 completed rows are represented in graph-pose metadata. The trace does not infer completion for the remaining 26.
- The historical Mol2 key is MCJFROJTDORRBB-ZWUDBOSPSA-N; the [manufacturer's SJF8240 #7266](https://www.bio-techne.com/p/small-molecules-peptides/sjf-8240_7266) key is SCUZPTNDZHEHLB-ZWUDBOSPSA-N. Formula and heavy-atom counts match, but connectivity differs.
- Saved GNN evaluation contains 17 predictions. The current static packet contains 46 graphs, including prediction targets with label -1; it is not a frozen copy of the evaluated training cohort.
- The [original-study supplementary workbook](https://ars.els-cdn.com/content/image/1-s2.0-S2451945617303549-mmc4.xlsx) supplies rows for 11/17 historical targets. MAP4K5 differs from the historical binary label for this assay. ABL1 has proteomics/western discordance explicitly discussed in the [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC5777153/). Six targets lack rows in that table; absence is not a negative label. Assay choice and compound/endpoint correspondence must be frozen before training.

## Implemented repairs

Exact stereo anchor mapping, explicit RDKit thread allocation and histogram
provenance extend the existing repaired pRosettaC derivative. The actual
200-conformer SJF8240 before/after control produced identical SDF bytes,
histogram and range. The new 5×5×2 SJF8240/MET pilot completed with no valid pose;
its receipt and every rejected search are retained. It is not a biological label.

New graph tools use executed atom mappings, completed-production receipts,
compiled-graph chemistry checks, periodic minimum-image distances, stable node
correspondence, and measured same-frame co-contact. Graph dimensions are retained
for the architecture, but the corrected feature definition is versioned and
requires fresh evaluation. Missing graphs remain explicit missing data.

The isolated GNN runner requires explicit inputs/output, reviewed experimental
labels and protein identity groups, verified chemistry/schema and input hashes,
deterministic CPU fold seeds, training-only scalar normalization, immutable output
directories and source/input drift checks. Unknown labels do not enter training.
A native-PyTorch runtime/gradient check passes without loading the incompatible
optional binary extensions found in the installed environment. No global package
or original model was changed. The supported minimal dependency approach is
explained in the [PyG installation documentation](https://pytorch-geometric.readthedocs.io/en/stable/notes/installation.html).

The existing generic CE and separate pocket-comparison banks cover only one and four of the historical 17 targets, respectively; the known-target control bank covers MET. These banks are not a completed or adopted 17-target study preparation. Missing candidates are not biological negatives. The target-specific preparation review must be completed after the assay/cohort is fixed.

## What is still required

The historical trajectories cannot be converted into corrected chemistry by
editing their ITP files or resuming old checkpoints. Corrected systems and
chemistry-dependent features must be recomputed. Protein structures, software
architecture, genuine experimental records and historical raw results can be
retained with their provenance. Sampling and convergence still need evidence.

`sjf8240_review/sjf8240_rerun_manifest.tsv` lists every historical target and its
specific assay-review/MD status. `publication_readiness.json` records the remaining
release requirements. No training release or Ernie campaign has been issued.
See `SJF8240_REMEDIATION_LOG.md` for what changed, how it was checked and why.

The second, adaptive SJF8240-only diagnostic completed: 1,680 global trials,
50 accepted arrangements, 250 local poses and 250 explicit rejected linker
searches; no QC-valid full pose was found. Its receipt was independently checked.
The smaller no-pose run remains retained. This has not validated the corrected
SJF8240-to-MD handoff, and no biology label is inferred.

Focused controls show that each experimental head can be embedded independently;
the intended mapped global anchors and sampled distance interval are correct.
A guaranteed-feasible free full-molecule head geometry is reproduced on the first
attempt. That self-control does not use experimental head shapes or a protein
complex. The recorded global/local frames did not provide successful two-bound-
head linker geometry within this budget; further SJF8240 pose/sampling review is
still required before a scientific campaign.


### 2026-11-06 experimental linker/MD checkpoint

- The active isolated workflow now uses the centered, torsion-preserving fixed-head correction; the earlier default and both comparison variants remain in Git history and frozen controls.
- Paired 50-frame diagnostic: baseline 0 valid, centered-DG candidate 15 full-QC valid, identical local protein coordinates and 0 Å head-core displacement.
- GSK215 candidate integration check: 2.209 Å VHL RMSD, 0.75 native-contact fraction.
- Candidate-specific geometry tests: 3 passed in both runtimes. Active workflow now has 119 tests, including parameter-bundle completeness and modeled-loop interchain overlap rejection; the 105-test software checkpoint is historical.
- The broader-DG MD control passed repair, topology, EM (Fmax 172.45 <500), 200 ps NVT and 500 ps restrained NPT. The subsequent 500 ps stage with ligand restraints released is running. Exact selected-method handoff is being tested separately. No production convergence or biological accuracy claim is made.
- All immutable parameter-generation artifacts restored; validated docking overlay preserved privately inside the derivative folder. Base conda dependencies remain external.


### Existing-results audit after user scope clarification

The read-only cutoff/selection audit reconciles all 141 core completed MD cases to their original saved CSVs/headers. A separate derivative sensitivity audit uses common 0–20 ns windows and explicitly removes one identical duplicate in the extended AXL/pose011 record. Actual historical MD groups contain hydrogens. The sequence audit confirms that the recorded group labelled VHL is Elongin B in all 141 core systems; true VHL is chain C. The old table reproduces legacy-labelled group proximity, not VHL retention. Corrected contact re-analysis uses the same trajectories. Five angstroms is an operational cutoff, not a universal binding/degradation criterion. See `existing_result_audit/BRIDGING_RATIONALE_AUDIT.md` for exact definitions, cutoff tables, references and limitations. All original results remain unchanged.

The user-authorized selected-method minimization retry passed in 5,468 steps (Fmax 455.493 <500 kJ/mol/nm). No following MD stages were started for that retry. The earlier broader-DG control completed its staged MD plus a 20 ps diagnostic segment; neither control replaces the historical study.


### Completed existing-trajectory correction checkpoint

All 141 core trajectories were re-analysed without new MD. Sequence-identified true VHL is Protein_chain_C; every old group labelled VHL is Elongin B. All original GROMACS distance checks reproduce within 0.000022 Å. Separate legacy/all-atom/true-VHL/heavy-only definitions, cutoff tables and review figures are in `existing_result_audit/contact_core_v1/`. Original metrics and GNN predictions are preserved. One retained static structure per core target (17/17 representatives) has the same legacy Elongin-B selector; this is not an exhaustive per-pose static audit. The corrected E2 script's Elongin-C alignment anchor is consistent with the sequence check; its conditional metric semantics still need precise descriptions.

The retained E2/GNN scalar and availability audits are recorded in `existing_result_audit/E2_GNN_SCALAR_LINEAGE.md` and `GNN_AVAILABILITY_AUDIT.md`. They preserve old predictions, flag one scalar estimator inconsistency, and show that incremental GNN value needs availability controls. No model was fitted for these audits.
