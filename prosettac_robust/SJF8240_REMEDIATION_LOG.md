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
