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
