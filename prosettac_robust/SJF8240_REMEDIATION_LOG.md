# Isolated SJF8240 remediation record

## 1. Scope and preservation

**What:** Prioritized the existing SJF8240 pipeline before other compounds.
**How:** All new changes are confined to prosettac_robust. The D: code/doc snapshot
has 105 SHA-recorded files. Original simulation, graph, prediction and parameter
files were read without alteration. Existing baseline workflow is retained in the
parent repository; derivative changes and backups are separately recorded.
**Why:** Preserve completed work and make every correction reviewable.

## 2. Chemistry and actual simulation lineage

**What:** Traced the constitutional-isomer mismatch beyond an input filename.
**How:** Compared canonical graph keys, element-labeled graph isomorphism, 255
ITP connectivity signatures, 174 actual production TPR graphs, manifest states
and the 174 per-pose MD graph metadata rows. Kept exact SHA-256 fingerprints and
machine-readable inventories. The GAFF-typed Mol2's direct RDKit parsing failure
is recorded as a format limitation, not evidence of corruption.
**Why:** Same formula and 79 heavy atoms concealed changed bond connectivity.
Published SJF8240 chemistry claims require the correct molecular graph.

## 3. Initial distance sampler

**What:** Replaced MCS anchor matching with complete disjoint stereochemical head
matching and a unique full-molecule anchor pair. Removed hidden four-thread use.
**How:** Added distance_anchor_indices, updated both adapter copies and the runtime
manifest, and retained an additional patch relative to the previous checkpoint.
Tests reject ambiguous, out-of-range and inverted-stereo mappings. A 200-conformer
paired control in the actual Python 3.10 docking runtime produced identical
SJF8240 conformer bytes, histogram and range. A per-run sampling ledger records
seed, input hashes, atom mapping, histogram, counts, range and thread settings.
**Why:** Avoid silent atom substitution and unexpected resource use while retaining
the working algorithm. This heuristic range is not exhaustive linker reach.

## 4. Current pRosettaC integration

**What:** Ran a bounded corrected SJF8240/MET control from frozen sources.
**How:** Five global poses, five local poses each, two requested conformers per
local pose, two workers, seed392 and 20,000 maximum global trials. All 25 local
poses reached explicit linker accounting; no QC-valid full pose was found. The
terminal receipt was independently verified and raw files retained. The static
builder's real empty-ensemble path was exercised with label -1.
**Why:** A reproducible zero-pose outcome is valid execution accounting and missing
prediction coverage; it must never become an experimental nondegrader label.

## 5. Graph feature corrections

**What:** Corrected periodic contact calculation, joint contact edges and atom
mapping in a separately versioned feature schema.
**How:** Production graph extraction requires source-bound preparation and stage
receipts, a mapped solute topology, actual TPR connectivity matching verified
SJF8240, unique heavy groups, stable residue correspondence and increasing times.
Uses triclinic minimum-image ligand-heavy/Cα distances and observed simultaneous
contact, averaged equally across poses. The old product of marginal contact
frequencies was a proxy, not observed co-occurrence. Tests distinguish perfectly
anticorrelated contacts and a boundary-crossing periodic contact. Static graph
construction verifies pose receipts, checks node identities across poses and
uses conditional contact distance means; graph packing checks tensor shapes,
finite values, indices, labels and per-packet hashes.
**Why:** Correct the physical calculation and prevent silently mixing old/new
schemas, changed node correspondence or incomplete simulation artifacts.

## 6. GNN cohort and execution

**What:** Isolated the evaluated cohort from mutable prediction packets and
preserved the model architecture for fresh corrected-data evaluation.
**How:** Recorded the current 46-entry static and 10-entry dynamic packets and
saved 17-target evaluation. The derivative requires reviewed labels, unique
protein identity groups, source-provenance hashes, explicit seed and output,
training-only normalization, finite values, deterministic CPU execution and
before/after source/input checks. It excludes unknown -1 labels and detects
aliases of the same protein crossing the LOPO boundary. Output directories
cannot overwrite historical work or prior evaluations.
**Why:** A changed graph packet or alias is not a new independent training target.
Scores are exploratory outputs, not calibrated degradation probabilities.

## 7. GNN runtime compatibility

**What:** Found an actual optional-extension import crash in the installed GNN
runtime; preserved the original installation.
**How:** Early direct imports warned about pyg-lib and aborted in compiled code.
The derivative configures the process to use native PyTorch operators before
importing PyG, excluding optional incompatible binary extensions. Verified both
command-line entry points and a synthetic full dual-encoder forward/backward
optimization step on CPU with a deterministic seed. Saved versions and control
receipts. No package installation or global environment mutation occurred.
**Why:** A reproducible Python dependency list is insufficient when optional
binary libraries crash at import. This control demonstrates execution only.

## 8. Experimental label audit

**What:** Compared historical labels with an explicit published assay and retained
uncertainty rather than deriving labels from predictions.
**How:** Downloaded the original publisher workbook, checked that it is a real
XLSX, recorded its URL/hash, and extracted the 17-target comparison. PMC first
returned an HTML download challenge; Europe PMC returned a non-open-access API
error. Those responses are retained as failed retrieval evidence; neither was
parsed as assay data. The publisher download succeeded. MAP4K5 is assay-discordant,
ABL1 has proteomics/western discordance and six targets are absent from this table.
Contradictory CDK17 statements in old notes do not substitute for assay evidence;
its measured abundance in this table supports a negative for this specific
compound-1 endpoint. Historical labels/results remain unchanged. Assay selection
has been requested from the user while independent software work continued.
**Why:** Cell line, concentration, duration, controls and endpoint define the label.
Changing a label because of geometry or model performance invalidates evaluation.

## 9. Validation, remaining work and publication limits

Regression logs, runtime versions, source fingerprints, paired sampler control,
compiled-chemistry rejection control and GNN gradient checks are saved under
validation_controls. Both GSK215 recovery seeds are retained: v26 missed the joint
criterion; adaptive v27 recovered one QC-valid pose at 1.889 Å and 9/16 contacts.
This single-reference result is not SJF8240 recovery or convergence evidence.

The corrected SJF8240 scientific campaign is not released. Assay/label resolution,
corrected pose coverage, rebuilt production simulations with independent replicas,
convergence assessment, new graph/E2 features and frozen evaluation remain required.
The original data is retained as historical data; no old trajectory is rewritten
or reclassified as corrected-chemistry MD. A small retrospective cohort does not
establish prospective cross-compound accuracy or guarantee journal acceptance.


## 10. Portable sequence data and exact MD frame binding

**What:** Closed a dependency gap and a pose/protein handoff gap before new MD.
**How:** Selected exact cached official UniProt entries for the 17 historical target
proteins plus VHL/Elongin B/Elongin C into a hash-bound 20-entry public cache,
retaining source batch URLs and hashes. Added canonical protein groups to the
reference cohort. Verified current SJF8240/MET chains A/B/C/D map with 100%
identity (140/101/85/278 observed Cα atoms). Protein completion now refuses to
overwrite an existing model/repair receipt and rejects ligand overlap below 1.5 Å.
MD preparation requires the repair receipt's source hash to match the exact
selected pose, in addition to validating the repaired protein's output hash.
A regression test rejects a repaired protein from a different arrangement.
**Why:** An untracked local sequence cache prevents independent reproduction;
a protein repaired for one pose cannot silently accompany a ligand from another.

The final source checkpoint passes 105 tests in each validation environment.
Graph and scalar packets must explicitly declare the corrected chemistry/schema;
matching tensor dimensions cannot silently authorize reuse of historical packets.
All tests and the synthetic GNN control are separate from scientific performance
validation. An additional SJF8240 50×5×2 bounded diagnostic is active from frozen
sources; no outcome or end-to-end MD success is claimed before it completes.


## 11. Closed the larger SJF8240 diagnostic without relabeling its failures

**What:** The 50×5×2 adaptive diagnostic completed with zero valid poses.
**How:** Independently verified its terminal receipt and immutable sources/inputs.
All 250 local poses reconcile to linker records. In a focused first-pose test,
three seeds succeeded for free embedding, three for VHL-only and three for
POI-only constraints, while all three two-head cases failed initial minimization.
A separate known-feasible full-molecule free-conformer geometry succeeded in one
attempt. Checked all 50 initial global poses against named anchors from parameter
provenance: their distances stay within the sampled 7.911–11.911 Å interval,
and the legacy serial-selected anchors agree with those names.
**Why:** These checks distinguish a finite failed protein-arrangement search from
an intrinsically unusable molecule, misselected anchors or a universal embedding
failure. They do not demonstrate a feasible experimental-bound ternary complex.
The SJF8240-to-MD handoff remains unvalidated by these runs; no artificial protein
complex or successful biological prediction is substituted. Keep raw failures,
focused controls and the sampling limitation visible before expensive production.


## 12. Paired bound-head linker-method diagnostic and first corrected-chemistry MD handoff

**What:** Isolated a constraint-embedding failure and tested a candidate method without changing the default production workflow.
**How:** Centered both bound heads before embedding; tested torsion-preference/basic-knowledge switches independently. A centered variant with both preferences disabled generated candidates while retaining exact stereochemistry, fixed heads, UFF refinement, original coordinate restoration, raw scoring and every existing structural cutoff. Across 50 local protein arrangements whose atom coordinates exactly match the baseline run, baseline ETKDG yielded zero valid poses; the candidate yielded 29 linker-feasible/scored poses and 15 passing full QC. Maximum head-core displacement was 0 Å. The candidate native GSK215 control retained 2.209 Å VHL RMSD and 0.75 native-contact fraction. These are integration controls, not independent prospective docking validation.

Three candidate-specific tests cover bound-head preservation, large coordinate translation and bounded rejection; each passes in both analysis and approved docking runtimes. A separate flag-isolation control suggests retaining experimental torsion preferences may also work when basic-knowledge terms are disabled; that alternative is not yet qualified or adopted. Candidate code, tests and compact evidence are under `linker_method_review/`; the default remains unchanged pending method review.

Selected the first lexicographically sorted QC-valid SJF8240 pose for an unlabelled MD integration control (not energy ranking). Sequence-informed repair matched all four chains at identity 1.0, modeled internal gaps, kept observed atoms unchanged and retained a 1.634 Å minimum protein–ligand heavy-atom separation. Corrected-chemistry topology preparation passed for 10,278 solute atoms. Local EM is underway, followed by the full staged handoff and a planned 20 ps diagnostic production segment; this duration cannot establish study convergence.

**Why:** Paired identical protein frames isolate the linker-method effect. Success must retain chemical identity and head geometry rather than use a free-conformer fallback. A validated software handoff is required before an expensive study; these controls cannot establish biological degradation or publication readiness.

## 13. Restore complete immutable ligand-parameter bundles

**What:** Fixed a clean-copy packaging omission affecting the MD handoff.
**How:** Restored 18 exact log/serialized files across six parameter-state bundles from the retained generation outputs, checked every file against its original validation hash, and added a distribution-level regression check. Explicit Git exceptions retain these bound artifacts; serialized files are hashed as bytes, never loaded. No parameter charges, topology, validation hashes or chemistry were changed.
**Why:** Existing receipts declare all generation artifacts. Omitting logs/pickles made reuse fail despite valid topologies. Keeping their exact bytes preserves the original evidence and lets a clean checkout satisfy strict verification. The expanded workflow suite passes 106 tests in the analysis runtime; the former 105-test checkpoint remains historical.


## 14. Check modeled protein gaps against every other chain

**What:** Added explicit post-repair cross-chain heavy-atom overlap rejection, in addition to peptide closure and ligand checks.
**How:** Measure each chain pair with a KD tree, reject any distance below the same 1.5 Å gross-overlap threshold, and include distances in future repair receipts. A regression counterexample puts a newly modeled atom inside another chain and must fail. Independently checked the selected repaired SJF8240 model: the closest chain-pair contact is B/C at 2.296 Å, so this additional guard does not change its outcome. Its original receipt remains unmodified; the retrospective check is separate evidence.
**Why:** A loop can avoid the ligand yet overlap a different protein chain. Sequence identity and closed peptide bonds alone do not rule out this modeling failure. This is a gross geometry guard, not an assessment of loop conformational accuracy.


## 15. Preserve the approved local docking runtime

**What:** Saved the validated 3.6 GiB dependency overlay privately inside the derivative project and added an explicit launcher.
**How:** Guarded installation copied every selected dependency file by SHA-256 into `.runtime/python310-packages`. The public launcher reads private configuration, checks that Python, RDKit and the PyRosetta database exist, then prepends this overlay for the command and child processes. Verified RDKit 2025.09.5 supports coordinate maps, NumPy is 2.2.6, OpenMM is 8.4 and PyRosetta imports from the persistent private directory. Licensed binaries/database remain Git-ignored. Base Python 3.10.22 and fallback conda packages are still external dependencies.
**Why:** The bare base interpreter has an older RDKit lacking the required API, and `/tmp` is not durable storage. This preserves the validated local runtime without altering original environments. It is not an Ernie deployment or self-contained environment release.


## 16. Prespecified comparison retaining experimental torsion preferences

**What:** Started a second isolated candidate comparison, keeping experimental torsion preferences while disabling basic-knowledge terms and centering the fixed heads.
**How:** Froze separate source snapshots, copied the same validated inputs, retained global/local seeds, counts and conformer budgets, and registered SJF8240 global10×local5×conformer2 plus native-frame GSK215 controls before execution. The first candidate, baseline and their completed results remain unchanged. Default workflow is not switched. Full active run directories are excluded from installation until terminal receipts exist.
**Why:** A successful broad change does not show that every changed preference is necessary. This comparison isolates whether a smaller change retains compatibility with bound head geometry. Outcomes are pending; no increase in biological accuracy is inferred from greater pose yield.

Parameter bundles use Git attributes that preserve exact bytes and display serialized/log artifacts as binary; their validation hashes remain authoritative. Frozen candidate source retains inherited upstream formatting. This prevents line-ending conversion from invalidating chemistry provenance and avoids treating immutable generated logs as hand-edited source.


The manually copied torsion-candidate v1 source trees omitted `source_manifest.json`. Independent `verify_receipt` correctly rejected the completed GSK control. Both v1 attempts are excluded from method qualification and retained unchanged, including the bounded active SJF attempt. Created fresh v2 attempts through the canonical `snapshot_sources` function and bound each pre-execution manifest hash into the campaign intent before launching. No historical manifest was retroactively invented. This is a diagnostic setup error, not evidence of a failed molecule or biological negative.

The benchmark recorder now validates any supplied source snapshot and campaign manifest binding before starting a worker. Added a counterexample proving that a missing source manifest cannot start a process or create an attempt record. This turns the observed setup error into an inexpensive pre-execution failure for future runs. Existing recorded controls remain unchanged.

The corrected-chemistry SJF8240/MET MD integration system passed full EM: 5,273 steps, Fmax 172.44859 kJ/mol/nm (<500), with a verified completion receipt and compiled-input/output hashes. Started the full configured 200 ps restrained NVT stage on the local GPU. This is the first corrected-chemistry SJF8240 EM handoff in the isolated workflow; it does not validate the modeled loop conformations, candidate docking method or production convergence.


## 17. Select the minimal qualified linker correction for local software validation

**What:** Adopted the centered, torsion-preserving variant in the isolated active workflow; scientific production is still gated separately.
**How:** Retain ETKDGv3 experimental torsion preferences, disable only basic-knowledge preferences, center exact bound heads before embedding, refine with fixed-head UFF and restore original coordinates. Global/local protein sampling, stereochemical identity, raw scores, clash cutoffs and head-preservation checks are unchanged. The properly snapshotted v2 SJF control independently verifies: same 50 local protein frames, 29 linkers/scored, 19 full-QC valid. The GSK native-frame control retains 2.209 Å RMSD / 0.75 native contacts. Both variants succeed for all three previously accepted global GSK arrangements at all three seeds (18/18 checks total); that retrospective selection is not independent validation. Added real SJF bound-head fixtures and two default-workflow regressions for native and translated coordinates.
**Why:** Centering removes coordinate-origin dependence. Retaining experimental torsion preferences changes fewer modeling assumptions than disabling both preference sets. Successful typed/stereochemical/fixed-head controls, not degradation labels or model accuracy, justify the software correction. The broader candidate and all failed baseline runs remain available for comparison. Exact-input scientific validation, independent sampling and microstate decisions remain necessary before production.

The ongoing MD control starts from the first broader-DG candidate; no claim is made that it is byte-identical to a torsion-preserving output. Its chemical identity and topology handoff are valid integration checks, while method-specific validation is tracked separately.


## 18. Corrected-chemistry static graph and actual encoder integration

**What:** Built and packed a graph from all 15 QC-valid broader-DG SJF8240/MET poses and exercised the real static encoder.
**How:** Independent docking receipt/chemistry verification precedes graph creation. The resulting graph has 497 nodes and 2,052 directed edges with the versioned 9-node/3-edge-feature schema. Label stays -1 (unknown). Hash-bound packing and the untrained native-CPU encoder produce a finite 1×128 embedding at fixed seed 7919. No loss, training, degradation probability or biological prediction is produced. Compact receipts are versioned; raw graph packets remain local.
**Why:** Synthetic tensor tests cannot alone establish that real corrected-chemistry structures produce compatible model inputs. This checks the actual handoff without converting a diagnostic into an accuracy claim or importing historical labels.

## 19. Audit the chemical states actually used by the MD control

**What:** Recorded protein charges, all titratable residues, 19 histidine tautomer choices, disulfide topology and ion counts from the immutable prepared topologies.
**How:** Verify the preparation hashes, parse actual per-atom topology charges/bonds and sum the complete system. Protein charge is -5 e; neutral SJF8240 plus 291 Na+ and 286 Cl- yields net zero; water count is 99,872. State assignment is the GROMACS standard-state/geometry-based default. No pKa calculation or pH-dependent microstate validation is implied.
**Why:** A neutral vendor drawing does not establish the dominant solution microstate, and automatic histidine choices must be visible for scientific review. This audit records the exact control without changing a running system; protein/ligand microstates remain part of the scientific production review.

The broader-DG MD control has completed 200 ps NVT and 500 ps restrained NPT with proper receipts. Its next 500 ps stage releases ligand restraints while retaining protein restraints. The final planned 20 ps unrestrained production remains a diagnostic, not a converged study.


## 20. Remove stale graph assumptions from the active MD feature producer

**What:** Removed unused historical residue-jump heuristics, original-project default paths and hardcoded assay labels; corrected feature documentation.
**How:** Execution still uses the mapped periodic implementation. The module now has no implicit screen directory or biological labels. Documentation states the actual dynamic definition: 5 Å ligand-heavy-to-protein-Cα proximity; VHL–POI edges use observed same-frame joint proximity and mean distances over sampled frames. Node features are three-dimensional, and equal per-pose aggregation is explicit. No graph algorithm, threshold or model architecture changed. A regression check prevents implicit source-directory or label selection.
**Why:** Historical comments incorrectly described six node features, all-heavy-atom contacts, product-of-marginals edges and thermodynamic stability. Those statements do not describe the repaired implementation and could cause incorrect manuscript methods or reuse of original data.

The first lexicographic QC-valid pose from the selected torsion-preserving method (`combined_10_1_0001`) passed the current sequence/gap/interchain/ligand repair checks and corrected-chemistry topology preparation. Its ligand–protein minimum separation is 2.435 Å and observed protein atoms remained fixed. A separate local EM is underway; this is tracked separately from the broader-DG MD system.


## 21. Named-column thermodynamic diagnostics and a caught audit bug

**What:** Added descriptive last-half temperature/pressure/density/potential diagnostics for completed control stages; fixed an error in the first audit implementation.
**How:** The first audit assumed GROMACS emits energy columns in request order. It actually sorts by internal energy index, yielding visibly nonphysical reported temperatures. Retained that result as `md_thermodynamics_initial_INVALID_column_order.json` with an explicit INVALID status. A new named-column reader requires unique indexed XVG legends, finite increasing-time data and positive temperature/density. Regression tests demonstrate reordered terms and reject unlabelled/nonphysical series. Regenerated in a fresh `energy_audits_v2` folder; no raw EDR or MD system changed.
**Why:** Finite numbers and matching dimensions do not prove semantic column identity. Publication analyses must bind names to data rather than infer order from a command. Valid final-NPT last-half values are mean T 299.969 K (sample SD 0.491), density 1001.326 kg/m³ (SD 0.907), pressure -3.623 bar (SD 48.445). Samples are correlated; these descriptive statistics are not confidence intervals or evidence of equilibrium/convergence.

The broader-DG control completed the second 500 ps NPT stage and began the prespecified 20 ps unrestrained integration segment. It remains excluded from any scientific training/release cohort.


## 22. User scope clarification and persistent requirements

The user clarified that the existing SJF8240 structures and completed PRosettaC/MD/GNN results should be retained and reused, with primary attention on robust scripts and the correctness of rationales such as bridging <5 Å. A wholesale repeat of the expensive experiment is not wanted. They requested that all instructions be saved persistently. Created `USER_REQUIREMENTS.md` and agent-readable `AGENTS.md`, covering preservation, location, Tocris reference evidence, scientific interpretation, compute, Git, detailed notes and communications.

Stopped only the new torsion-control local GROMACS minimization by SIGTERM when this clarification arrived; no original job or result was touched. The user then allowed helpful work/diagnostics to continue, asked about the minimization's purpose and explicitly said to continue. A fresh bounded minimization retry will test the same diagnostic system while preserving the interrupted attempt. This is local software integration validation, not a rerun of the completed experiment or a study release.


## 23. Audit the completed experiment's bridging definition without rerunning it

**What:** Added a read-only producer for existing distance-table, index-group and cutoff-sensitivity audits. Reused the existing SJF8240 outputs rather than generating new study data.
**How:** The initial full-manifest audit reconciled 140 cases, rejected a duplicate-time AXL series, and explicitly recorded 33 prediction-target cases without the expected local XVG files. A separate, explicitly scoped core audit reconciles 141/141 original CSV/header results, then removes only the exact identical duplicate in AXL/pose011 and analyzes a common 0–20 ns window. That AXL record actually extends to 100 ns. All cases use 201 stored samples at 0.1 ns spacing in the common window. Counterexamples reject timestamp mismatch, truncation, wrong units, conflicting duplicates and confusion between joint occupancy and a product of marginals. Prespecified cutoffs 3.5–6 Å are evaluated without labels or performance-based selection; 5 Å remains the inherited primary rule.
**Why:** The scientific rationale depends on what is measured, not just the numeric threshold. Actual saved ligand groups include 65 hydrogens among 144 atoms. The pose-screen, MD and dynamic GNN contact definitions differ. The manuscript must distinguish these, avoid calling proximity a degradation probability, and report threshold/duration sensitivity. Full explanation and primary-source references are in `existing_result_audit/BRIDGING_RATIONALE_AUDIT.md`. Original coordinates, trajectories, metrics and predictions remain unchanged. Reprocessing stored XVGs is not an MD rerun.

## 24. Close the explicitly authorized minimization diagnostic

The fresh selected-method EM retry completed in 5,468 steps with Fmax 455.493 kJ/mol/nm (<500). The interrupted attempt is retained separately. This confirms a bounded local preparation/minimization handoff for that diagnostic pose. No subsequent MD stages were automatically started for this retry, and it does not replace or require repeating any completed SJF8240 trajectory.

The resulting active software checkpoint passes 117 tests in each of the analysis Python 3.13 and approved-overlay Python 3.10 environments. This is software verification, not a scientific-release or journal-acceptance claim. Existing results and historical source snapshots remain separate from future-method controls.


## 25. Independent protein-sequence check reveals a legacy role-label error

**What:** Compared the proteins actually present in the first existing DDR2 production TPR against the cached canonical sequences, rather than trusting chain/index names.
**How:** Protein_chain_A has 101 Cα residues and matches Elongin B/Q15370 with 100% identity and observed-sequence coverage. Protein_chain_B matches Elongin C/Q15369 (85 Cα), while Protein_chain_C matches VHL/P40337 (140 Cα), both at 100%. The existing group named VHL selects Protein_chain_A. The POI chain matches DDR2 with 98.5% identity/coverage. Started a read-only compiled-TPR role audit across all 141 core cases before generalizing this result.
**Why:** Reproducing a metric proves implementation consistency, not biological label correctness. The all-atom contact table must currently be described as legacy-named group proximity. True VHL is present in the simulation, so corrected analysis selections can reuse existing trajectories. No original index, topology, trajectory, metric or model result is overwritten; no MD rerun is inferred.


## 26. Complete the 141-system protein-role audit and begin contact re-analysis

The first sequence audit identified 106 systems and left 35 unresolved because the auditor stripped `_full` from target names whose canonical cohort keys retain it. Preserved that incomplete report and its exact script; fixed exact-key lookup and processed only those remaining 35 cases. Combined hash-bound reports identify all 141 systems: every legacy `VHL` group is Elongin B/Q15370; true VHL/P40337 is Protein_chain_C at 100% identity and observed sequence coverage. All 141 POI groups match their expected targets. No simulated protein is substituted or removed.

A direct re-analysis of the first DDR2 trajectory reproduces old periodic GROMACS distances to within 0.000014 Å. At 5 Å, the old Elongin-B-labelled metric is 1.99%, while true-VHL/POI occupancy is 100% (also 100% for heavy-only true-VHL/POI at this cutoff). This is one case, not a cohort or degradation conclusion.

The initial trajectory check rejected time alignment at an overly literal floating-point tolerance: four XVG timestamps differ from XTC times by 0.000001 ns plus binary rounding. Preserved the diagnostic failure. The qualified reader uses the precision of seven-significant-digit XVG times (5.1e-6 ns tolerance for this 20 ns window), far below the 0.1 ns sampling interval. It still requires equal frame counts, increasing aligned times and finite data; no truncation or interpolation is allowed.

MDAnalysis's default reader attempted an offset/lock sidecar in the original directory; the read-only sandbox prevented the write. Added a read-only XTC reader that bypasses sidecar creation, including refresh mode, and tests that input hashes/directory contents remain unchanged. Future re-analysis uses that reader explicitly.

Froze a 141-case, two-worker, common-20-ns re-analysis plan and all executable/helper/config/sequence source bytes in `existing_result_audit/contact_core_v1`. This recomputes distances and proxy definitions from existing trajectories; it executes no docking, parameterization, minimization or MD. Active batch files are excluded from installation until terminal completion. Original results and models remain untouched, and corrected analyses will be separate outputs.


## 27. Static graph role audit and bounded interpretation

Read the original static graph configuration through Python AST/literal data only, never importing its executable top-level code. One retained combined PDB per core target (17/17) identifies the old chain-A/residue-1–103 selection as Elongin B/Q15370. These representatives support the role-label finding; no claim is made to have reprocessed every static pose or rerun docking. Source PDB/script/cache hashes and selection rules are retained.

Checked the existing corrected E2 script separately: it explicitly aligns sequence-verified Elongin C, chain B/85 Cα. That anchor choice agrees with the sequence audit. Do not infer that every old E2 coordinate transformation is wrong merely from the mislabeled chain-A contact descriptor. The conditional bridging masks and biological meaning of reachability scores still require explicit definitions. Saved `existing_result_audit/METRIC_GLOSSARY.md` to distinguish legacy scaffold proximity, true-VHL contact occupancy, pocket engagement, E2 geometry and degradation labels/probabilities.


## 28. Complete contact re-analysis on all 141 existing trajectories

**What:** Recomputed the legacy-named, true-VHL all-atom and true-VHL heavy-only contact definitions on the completed original MD, preserving all old outputs.
**How:** Frozen source/plan, two analysis workers, common 0–20 ns windows, 201 samples/case, explicit sequence identities and strict read-only XTC input. All 141 cases complete; all 28,341 sampled frame records are retained in the stated windows. Periodic stored-group distances reproduce GROMACS to maxima 0.00002093 Å (legacy Elongin-B label) and 0.00002114 Å (POI). Independently checked every per-case success flag before aggregation. Produced per-pose/per-target cutoff tables plus vector PDF/SVG and PNG review figures with input/source hashes. Runtime about 499 s; no docking or molecular simulation was executed.
**Why:** The user wanted defensible rationales without repeating the computational experiment. True-VHL and heavy-only definitions change occupancy substantially, so the old descriptor must be named and interpreted as E3-scaffold proximity. High sampled proximity is not a degradation probability. Corrected analysis results are explicit derivatives, not silent replacements for old metrics or model inputs.

Updated the runbook/readiness metadata to prioritize reuse of existing results. Retained the old corrected-chemistry campaign requirements as a historical/conditional proposal, explicitly not authorization to execute a campaign. Existing chemical-identity and assay-label concerns remain scientific findings that code tests or contact re-analysis cannot erase.


## 29. Future GNN fold replay without repeating training

**What:** Added per-fold model-state checkpoints and full-precision score retention to the isolated GNN runner; no existing GNN was retrained or overwritten.
**How:** Save tensor-only model state, architecture, held-out/training identities, fold seed, exact training-only scalar reference rows, normalized held-out scalars, full score and source/contract hashes. Loadable with `torch.load(..., weights_only=True)`. Numeric fold filenames avoid target-name path issues; original-project destinations, overwrites and held-out/training overlap are rejected. A training-intent file records inputs/settings before computation, and only the terminal receipt certifies a completed evaluation. Final receipts hash checkpoint artifacts as well as prediction/summary tables. Invalid hidden/head dimensions and nonfinite learning rates fail before training. Console display remains concise, while tabular scores are no longer rounded before AUROC/metrics.
**Why:** Saved prediction tables alone cannot replay a model or recover discarded score precision. Checkpoints let future fold inference be reproduced without repeating training, while preserving the old outputs. Regression tests restore a saved model and reproduce its exact score; a four-case, one-epoch synthetic check verifies each held-out fold and score/checkpoint binding. These tests use synthetic tensors, not the SJF8240 dataset, and establish no biological performance.


## 30. Durable resume checkpoint after the session limit

The requested final save had not completed before the interruption: the last pushed commit remained 3dfbdae, and later work was partly staged. On 2026-11-07, checked actual files/Git status and created `RESUME_STATE.md` with scope, findings, exact artifact/runtime locations, remaining tasks and the incomplete-save history. The three new GNN tests require PyG, which is absent from the docking environment. They now explicitly skip there rather than failing or implying unsupported validation. Analysis-environment tests still require them to pass. No dependencies were installed into the original environments, no study data was retrained and no molecular simulation restarted.


## 31. Trace retained E2 inputs and expose an aggregation mismatch

**What:** Added a read-only scalar-lineage auditor for the current retained GNN feature table and source E2/static-lysine tables. No GNN retraining or molecular simulation was executed.
**How:** Reconstructed the declared builder's four scalars for the 17 reference targets and checked saved precision, missingness, denominators and selected static ranks. Preserved input/source hashes and the exact auditor source. One field, MAPK14_full/e2_delta_reach_80A, differs: saved 0.022575 versus paired within-pose reconstruction 0.02569375. The alternative difference of available pose means is 0.02256875, compatible within retained precision; its bridge and nonbridge sets have 18 and 16 poses. This is a compatible estimand, not proof of a manual edit or historical execution. Preserved the mismatch and original data. Produced an all-target paired/unpaired sensitivity table without labels or performance selection; missing-value matches are counted separately.
**Why:** The old feature descriptions incorrectly call an unconditional distance mean bridging-conditioned, and a rank-selected static distance an across-trajectory nearest distance. The builder's `.first()` has no mean fallback. Corrected descriptions in the isolated future runner while preserving its numeric architecture/algorithm and old results. Missing-value imputation is now described as a policy, not a guarantee against confounding. Full explanation is in `existing_result_audit/E2_GNN_SCALAR_LINEAGE.md`; exact historical training input remains unproven without its original execution receipt.


## 32. Audit incremental GNN evidence against availability flags

**What:** Compared the saved 17-target GNN predictions to fixed MD/E2 availability flags, without fitting any model or changing labels.
**How:** Verified unique binary-labelled records, finite saved scores and binary flags; hashed the source before/after. Retained GNN AUROC is 0.800 (15/17 correct); the MD flag gives 0.875 (14/17 correct), and E2 flag 0.7917 (12/17). All five recorded positives have MD; all nine without MD are recorded negatives. No thresholds or hyperparameters were tuned. This is an explicitly retrospective descriptive audit, not new independent validation.
**Why:** A complex model's mechanistic/incremental claim needs comparison with simple signals already encoded by workflow selection. Availability can be a legitimate cascaded-stage feature, so this finding is not automatically leakage or proof the GNN is useless. With the small selected cohort, role/selection issues and uncertain label provenance, report the historical result as exploratory and retain the baselines/limitations. No original GNN or molecular simulation was rerun.


## 33. Verify E2 static reference correspondence and legacy patch pathway

**What:** Checked the E2 producer's static coordinate-frame assumptions without recalculating MD or altering model inputs.
**How:** Compared ordered 85-Cα Elongin-C coordinates from all 141 retained starting structures and the producer-selected reference PDBs against the common reference. Random-pose raw RMSD is 0 Å. The JSON's 13 arrays are byte-identical, each with 51 positions; flattening creates 663 records but only 51 unique coordinates. Equal duplication does not change the source minima/fractions, but is not independent ensemble diversity. The source's 20 Å Gaussian perturbation model is recorded as an assumption, not equilibrium sampling.

Also inspected, without executing, `16b_mapk14_e2_fix.py`. Its body explicitly rewrites MAPK14's 80 Å scalar using unpaired available-means and six-decimal rounding, providing a source-level pathway consistent with the discovered estimator discrepancy. Execution timing/inputs remain unproven. Its header says retain LYS211 while the worker excludes it and upserts original tables if invoked; keep it as historical investigation code rather than a routine pipeline step. No residue selection or original table was changed.
**Why:** Distinguish actual numerical/frame issues from misleading comments or inflated sample-count interpretation. Common static coordinates and duplicate points are verified facts; dynamic coherent imaging and the prior's physical validity still need separate checks.

## 34. GNN-focused input and imaging checks

What/how/why: user clarified GNN priority; no conventional model was added. The bounded imaging probe reads three retained trajectories with the sidecar-free reader, extracts only the old rotation function via AST, and hashes inputs before/after. Initial attempt selected nonexistent MET/pose001 and failed before writing output; the retained completed-case plan identifies MET/pose002, used for the successful probe. Fifteen frames show zero direct-versus-periodic protein-distance difference, proper rotations, and maximum ELOC fit RMSD 2.1213 A. This small probe does not settle rare imaging events or validate a virtual E2 periodic-distance replacement. No original file, simulation or prediction was changed.

Normalization previously silently omitted absent training identifiers through an index filter. It now requires explicit rows for all training targets and unique declared scalar columns; explicit missing values retain the documented training-only imputation. Regression covers missing-row rejection versus permitted NaNs. This prevents an unnoticed change in the GNN fold reference population.

## 35. Scientific sign-off review for the confirmed 19-target panel

User explicitly requested scientific sign-off and confirmed 17 historical evaluation targets plus SIK3/RPS6KB2. Read-only coverage reconciliation finds 174 MD graph pose records, each with 21 sampled frames, across 10 targets; 141 corrected contact cases cover eight core targets. The two additional targets have prediction calls, not validated biological negatives. The endpoint question is pending.

New audit hashes source tables, rejects frame/time drift, and reports true-VHL/POI heavy-atom simultaneous-contact occupancy at strict 4/5/6 A, four equal 50-frame 5 ns blocks (excluding initial frame for blocks only), per-pose variation and leave-one-pose-out mean ranges. These are descriptive sensitivity checks without confidence intervals, independent-replica assumptions or fitted models. EPHA2 has the largest block-mean range, 18.67 percentage points. No equilibrium conclusion is inferred in either direction.

Downloaded reviewed canonical SIK3/Q9Y2K2 and RPS6KB2/Q9UBS0 JSON from UniProt, retained exact bytes and merged with cached core sequences. Initial extension role audits v1/v2 failed because expected analysis index files are absent; preserved both reports and exact source versions. The compiled-only v3 identifies all four expected proteins in all 33 remaining TPRs without claiming verification of nonexistent legacy indices. No trajectories were modified or simulated.

SCIENTIFIC_SIGNOFF.md records the actual decision: unconditional approval withheld. Compound connectivity, assay definition, historical GNN exact replay, corrected feature evaluation, sampling/uncertainty and external validation remain scientific gates. Bounded software successes are explicitly separated from those gates. Automated syntax/coverage/cutoff-monotonicity checks passed. All originals and historical predictions remain preserved.

## 36. Retained-structure scope and scientific-method code fixes

User clarified again to use the original experiment's structures and focus on scripts/rationale. Persisted in AGENTS/USER_REQUIREMENTS/RESUME_STATE. Supplier comparison remains documented provenance, not an instruction to replace historical chemistry. Added HISTORICAL_METHODS_AND_RATIONALE.md and retained-data workflow entrypoint, without weakening prospective supplier-specific chemistry guards.

What/how/why: shared executed heavy mapping now rejects fractional/boolean/empty IDs, duplicate/overlapping groups, solvent IDs beyond the mapped solute and nonfinite/hydrogen masses. Contact analysis rejects nonfinite time controls/timestamps/coordinates, malformed periodic boxes and input drift; mapped GNN loading and contacts use sidecar-free XTC reads. Thermodynamics now reuses the indexed-legend parser instead of maintaining a weaker duplicate parser. Regressions test corrupt mappings, NaN trajectory time before discard and permuted GROMACS energy columns. Initial energy test fixture had too few samples per complete block; that failure is retained, and the fixture window was corrected without weakening the production check.

The new 09_reanalyze_existing_contacts.py supports retained core-layout trajectories directly, checking both stored VHL/POI time axes before combining their arrays. This fixes a timestamp-check omission in the older review helper. It records sequence/configuration/helper hashes, duplicate/time/window policy, historical chemical scope and immutable output. On DDR2/pose001 its 201-row contact table exactly matches the already preserved corrected table. The two extension targets lack the legacy index/XVG layout; this entrypoint fails explicitly rather than inventing mappings or replacing data.

Read-only compiled-input inspection of all 174 production TPRs used GROMACS dump -om, two workers, original TPR hashes before/after, retained input-record exports/logs, no grompp/mdrun. The first probe with local GROMACS 2024.3 failed on the newer TPR format; retained log explains failure. The original script identifies compatible GROMACS 2025.4 under .conda/envs/gromacs. All 174 input records specify 2 fs, nominal 20 ns and 100 ps compressed output. Two actual cutoff protocols exist: 113 cases at 1.0 nm, 61 at 1.6 nm; AXL has eight/ten respectively. Current template is 1.6 nm and cannot document all old runs. Methods now disclose this heterogeneity; it does not automatically establish numerical bias or require a rerun. Requested duration remains separate from actual trajectory duration (including the retained extended AXL analysis).

Exact retained source snapshots are frozen with original paths/hashes and were not executed. The local adapter is random rigid-body global sampling plus PyRosetta, not execution of PatchDock despite output folder names. Old source comments equating geometric fractions to productive/degradation probabilities are corrected in the methods companion, without modifying or retrospectively relabelling original outputs. GNN-only scope remains; no conventional predictive model was added.

Final retained-scope verification: analysis environment 129 passed; docking environment 126 passed/3 GNN-only skips. Full logs and all active script hashes are recorded in software_checkpoint_v7.json. The revised mapped contact analyzer also ran successfully on the already completed separate 20 ps diagnostic control, marked DIAGNOSTIC_ONLY, with no new simulation. This validates the real mapping/receipt/trajectory handoff without claiming it represents the historical study. New retained-data CLI on the original DDR2 case reproduces the old corrected contact table exactly. No GNN was trained or original project changed.
