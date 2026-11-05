# 10PROTAC-MD remediation record

This document records the work performed, its rationale, the evidence used, and
its limits. It is a scientific and implementation record, not a production release.

## Scope, authorization and current decision

The user requested a publication-quality audit before committing months of Ernie
computation, then authorized repairs. The user subsequently specified Tocris as
the source of purchased compounds and requested notes of everything done, how,
and why. No production simulation, docking campaign, model prediction, or remote
write has been performed during this work. Production remains **NO-GO**.

Repository inspected: `/home/ayushregmibagale/10PROTAC-MD`.
Initial committed baseline: `282f30cf3f393e1b4fe3feea082a0df6392b334a`.
The original working tree had untracked `prosettac_inputs/`. Those inputs were
preserved. No git reset, checkout, commit, push, or deletion of user data occurred.

## 1. Initial audit: what was inspected and why

Inspected all five original project scripts, README, all three methods documents,
structure manifests, reference PDB headers, warhead SDF identifiers, and all 2,860
staged directories. Also read the local external pRosettaC runner, adapter,
resume script, PRosettaC main/parameter reader and relevant PROTAC-MD geometry
and XGBoost training code. External code was inspected without changing it.

The original audit tool read PDB/SDF coordinates and measured minimum heavy-atom
protein/head distances using NumPy; it also checked configs, chains and hashes.
Distance thresholds (>8 A remote; <1 A severe overlap) were diagnostic screens,
not validated binding criteria. Initial artifacts were written under
`/tmp/10protac_preflight_20261104` and copied to the Windows-accessible folder
`/mnt/c/Users/Ayush Regmi Bagale/10PROTAC-MD-preflight-20261104`.
The complete report is `READINESS_REVIEW.md`; row-level evidence is in
`input_audit.tsv` and aggregate counts are in `summary.json`.

| Group | Staged pairs | Confirmed issues from the initial scan |
|---|---:|---|
| GSK215 | 650 | All configs incomplete, all POI chains wrong, all heads >8 A away |
| JNJ1013 | 644 | All configs incomplete and POI chains wrong; 522 remote; 29 severe overlaps |
| SJF1521 | 644 | All configs incomplete and POI chains wrong; 574 remote; 15 severe overlaps |
| SJF8240 | 461 | All configs incomplete; 417 severe overlaps; 1 remote |
| SJF_delta | 461 | No configs; inherited geometry duplicates SJF8240 |

All 2,399 existing configs pointed to SJF8240's full SMILES. The local runner
requires `Protac`, `Anchor atoms`, and `Full`; the configs instead supplied
`Smiles` and omitted the other required values. The default cluster was PBS,
inconsistent with the intended LOCAL adapter. The 1,938 newly written POI
proteins used D while their configs requested A. Missing configs were silently
excluded by the old launcher.

## 2. Chemical and structural source verification

Why: named files and parent inhibitor labels were not reliable molecular identity.
Purchased PROTACs differ in linker, modified head, stereochemistry and exit vector.
The Tocris catalog was therefore made the reference for full molecular graphs;
deposited crystal structures supply bound coordinates independently.

Supplier molecular structures were transcribed as isomeric SMILES into
`config/molecules.json`. RDKit parsed each graph, calculated molecular formula
and stereochemical InChIKey, compared both to the supplier values, and rejected
unspecified stereocenters. Canonical full SMILES were written separately under
`molecules/<compound>/full.smi`, never by copying SJF8240 to another chemotype.
This validates the catalog graph; it does not certify the user's particular vial,
batch hydration state, purity, or solution protonation state.

| Compound | Tocris catalog | CAS | Source |
|---|---|---|---|
| GSK215 | 7818 | 2743427-26-9 | https://www.rndsystems.com/products/gsk-215_7818 |
| SJF8240 | 7266 | 2230821-68-6 | https://www.rndsystems.com/products/sjf-8240_7266 |
| SJF1521 | 7261 | 2230821-40-4 | https://www.rndsystems.com/products/sjf-1521_7261 |
| JNJ1013 | 8029 | 2597343-08-1 | https://www.rndsystems.com/products/jnj-1013_8029 |
| SJF_delta | 7267 | 2254609-23-7 | https://www.rndsystems.com/products/sjfdelta_7267 |

SJF_delta was initially unresolved and the user was asked for its catalog number.
Further search found the actual Tocris 7267 entry, so the identity was resolved
without assuming it was a linker-length variant. Its VHL ligand/attachment differs
materially from SJF8240. The first installed version records it as unresolved;
the subsequent staged version verifies its formula and InChIKey successfully.
The installed `molecules/verification.json` is the authority for current installed
state; the later section below records the second installation when applied.

Confirmed reference errors:

* Local 5JC7 header and RCSB identify chicken MDA5 helicase, not FAK.
* Local 6EGA header and RCSB identify IRAK4 with a type II inhibitor, not IRAK1
  with JH-I-25. The correct parent-ligand IRAK1 source located was 6BFN/DL1.
* The inherited VS4718 SDF gives CID 25151352. PubChem identifies that as
  pexidartinib, not VS-4718/PND-1186 or defactinib.
* GSK215 has a native ternary structure, 7PI4 with CCD component 7QB.
* SJF1521's full conjugated head differs from parent lapatinib. A parent reference
  alone does not establish the coordinates of every atom in the actual head.

Downloaded RCSB mmCIF sources and CCD chemistry files using HTTPS curl with
`--fail --location --max-time 30`. The first sandboxed request failed DNS; the
permitted network escalation succeeded. Downloads are archived rather than
replacing the inherited references. File SHA256s and URLs are in
`reference_sources/sources.json` after installation.

| Reference | Role | Treatment |
|---|---|---|
| 7PI4 / 7QB | Native GSK215–FAK–VHL ternary | Bound ligand reconstructed, full identity checked |
| 1XKK / FMM | EGFR–lapatinib parent complex | Parent reference only, no claim of exact PROTAC head |
| 6BFN / DL1 | IRAK1–JH-I-25 parent complex | Correct parent reference replacing the incorrect assignment |
| 3LQ8 / 88Z | MET–foretinib parent complex | Fresh deposited source; does not trust similarly named old fragments |
| 4W9H / 3JF | VHL/VBC–VH032-related ligand | Separate parent reference for the unmethylated VHL head |

Extraction method: BioPython reads label asym IDs (7PI4 author chain IDs are
multi-character); ligand atom names are joined to CCD atoms; CCD bonds/charges
are used to build an RDKit graph; deposited coordinates are assigned by atom
name; stereochemistry is assigned from 3D; missing heavy atoms cause failure.
Protein-only output excludes water and unrelated heteroatoms. Label-chain mappings
are explicit, including VHL A, Elongin B B, Elongin C C, and POI D in new outputs.
Hash and atom-index provenance accompany every extracted reference.

The reconstructed native GSK215 InChIKey exactly matches Tocris 7818. Its minimum
heavy contacts are 2.503 A to VHL and 2.973 A to FAK, without a <1 A overlap.
This is a preparation control, not a redocking result or prediction.

## 3. Repaired GSK215 head and coordinate placement

Why: the old algorithm paired the first min(N,M) C-alpha atoms regardless of
residue correspondence and transformed an unfitted PubChem conformer as if it
were already in the crystal coordinate frame. The inverse Kabsch expression was
not itself the principal error. The computed RMSD was also divided by sqrt(3)
relative to conventional per-atom RMSD and never enforced.

New GSK215 heads are cropped directly from the deposited full bound molecule.
The intended cuts are CCD **N4-C20** and **C21-N5**, leaving the central CH2-C(=O)
linker outside both heads. The POI head contains 36 heavy atoms, the VHL head 31,
and the omitted linker 3. Each head is checked as a stereochemically compatible
substructure of the exact full PROTAC. The anchors are one-based indices 34 (POI,
CCD N4) and 5 (VHL, CCD N5), not an unverified generic `6 1`.
During development an initial bond-name draft would have cut a piperazine ring;
it was corrected before installation. Tests now require ring preservation, the
expected partition and unchanged native coordinates. No incorrect intermediate
head was installed in the project.

The replacement placement script requires an explicit table of source chain and
domain boundaries. It aligns sequences with BLOSUM62 correspondence, then fits
a structural core iteratively, records the exact mapped residue pairs, and checks
pocket coverage and pocket RMSD separately. Reference head coordinates transfer
with `(reference_head - t) @ R` under `target @ R.T + t = reference`.

Current provisional filters: 100-400 domain C-alpha atoms; >=100 corresponding
residues; identity >=25%; retained core >=60% of aligned residues; fitted core
RMSD <=2.5 A; pocket coverage >=80% with >=10 mapped pocket residues; pocket
RMSD <=2 A; minimum head/protein heavy contact between 1.5 and 4.5 A.
These are conservative candidate filters, not benchmark-certified accuracy claims.

Outputs go to `prepared_inputs/`; legacy inputs are not overwritten or considered
complete merely because a warhead SDF already exists. The native FAK self-placement
control achieved pocket RMSD approximately 1e-14 A. Its preflight has only
`BENCHMARK_RELEASE_PENDING`. It must remain blocked until independent benchmarking.

## 4. Execution, parsing and analysis repairs

| File | What changed | Why |
|---|---|---|
| `scripts/lib/protac_qc.py` | Strict parameter parsing, atom selection, score parsing, geometry and hash receipts | One shared contract prevents preparation/execution/analysis disagreements |
| `01b_preflight.py` | Includes missing configs; validates chemistry, chains, 3D heads, contacts and file hashes | Invalid jobs must not be submitted or silently disappear |
| `02_batch_prosettac.py` | Removed unsafe nohup submission; read-only plans; locked atomic manifest union; changed-input invalidation | Old code could lose filtered rows, duplicate work or misclassify partial results |
| `03_record_attempt.py` | Isolated UUID attempts, immutable inputs, timeout/process-group termination, terminal receipts | A zero return code or score file appearing is not proof of scientific success |
| `01_extract_bridging_features.py` | Uses explicit local attempt index and validated PDB geometry; preserves units/missingness | Rosetta score.sc has no bridging distances; guessed columns manufactured false negatives |

The new parser skips SEQUENCE, accepts repeated identical SCORE headers, and rejects
truncated rows, nonfinite scores, duplicate/unsafe pose IDs and missing structures.
Every scored combined pose must match a generated PDB. A possible adapter -1.0
placeholder and shifted-score artifacts are rejected for explicit reconciliation.
This is intentionally conservative; -1.0 is not a universally impossible physical
score, but is ambiguous given this adapter's known fallback.

`completion.json` requires exit code zero plus successfully validated outputs and
records score and per-pose hashes. Analysis checks it and the original input hashes
again; post-completion modifications cause failure. Missing data never becomes a
zero-hit biological negative. Explicit VHL/POI/ligand selections prevent residue
number overlap from generating artificial bridging. Distances are heavy-atom
minima, with strict <5 A and <8 A comparisons; fractions are 0-1 and bridge_pct
is 0-100. The extracted geometry table is not claimed to be the full input schema
for the existing MD/E2-dependent trained model.

The old score-file polling path was wrong for the inspected runner and insufficient
for completion even if corrected. The old extractor still expected a two-stage
manifest after launch had become single-stage, and its rsync filters/caching could
miss nested outputs. The repaired path removes those assumptions rather than
pretending transferred score.sc alone is a complete result.

The exact remote engine remains unverified. Local inspected external code replaces
PatchDock with random rigid-body sampling, can copy unrelaxed outputs with a fake
-1.0 score, shifts scores in resume, and does not consistently seed all RNGs.
It has not been silently patched or described as validated upstream pRosettaC.

## 5. KLIFS acquisition/selection repairs (staged after first installation)

The old download success counter could increment without finding the filename it
expected from OpenCADD. It checked a different resume filename from the files
actually on disk. The staged repair isolates each API request in a temporary
folder, requires exactly one protein PDB, copies it to a canonical KLIFS-ID name,
and writes a hash/ID/chain/DFG sidecar. Existing outputs without matching provenance
are not silently accepted. Existing manifests are not overwritten by a new query.

The old metadata join used gene+PDB only, ignoring chains/alternate structures.
The staged replacement joins by exact sidecar KLIFS ID when available, or accepts
only an unambiguous gene+PDB+chain match. Ambiguous records are explicitly marked
and given no inferred DFG state. File hashes are verified against provenance.
The UniProt-only deduplication that could discard tandem kinase domains is replaced
by explicit UniProt+domain identity, with duplicate errors instead of silent loss.
No inherited manifest has been regenerated or overwritten with unverified metadata.

## 6. Installation and validation evidence

All mutations were first developed under `/tmp/10protac_remediation`, because
`/home/ayushregmibagale/10PROTAC-MD` is outside the writable sandbox. Installation
required filesystem escalation; it was approved. This is a filesystem restriction,
not a request to approve the science or launch a campaign.

The installer compares staged source hashes and expected original file hashes,
backs up replaced files, copies only the concrete repair set, regenerates supplier
and reference artifacts in the actual project paths, runs preflight, and runs tests.
It performs no SSH, simulation, git commit or push.

First installation backup: /home/ayushregmibagale/10PROTAC-MD/audit_baselines/before_remediation_20261104T130122Z
Machine installation record: `audits/remediation_install.json`.
Regression output: `audits/regression_tests.txt`.
Observed local package versions: `config/environment_observed.json` (not an Ernie lockfile).

Validation: 32 regression cases passed in the installed project. Coverage includes
catalog formulas/stereochemistry, GSK215 native identity, correct head partition,
native coordinate preservation, contact geometry, overlapping/absent atom
selections, Rosetta headers, malformed/truncated/nonfinite/duplicate scores,
missing PDBs, placeholder scores, missing completion receipts, output/input tampering,
nonzero exit, timeout, zero exit without output, manifest filtered-union behavior,
Kabsch inverse transform, domain correspondence, missing configs and refusal to
substitute another chemotype. The synthetic worker tests run tiny Python commands,
not molecular simulations. Native self-fit is not an independent accuracy benchmark.

Installed legacy preflight: 2,860 inspected, 0 READY, 2,860 BLOCKED.
Installed new candidate preflight: 1 inspected, 0 READY, 1 BLOCKED only by benchmark
release. These counts do not mean the legacy inputs were repaired: they mean the
new gate correctly refuses them while replacement candidates are being prepared.

## 7. Connectivity and interrupted work

A read-only SSH scheduler inspection was attempted with BatchMode and a 10-second
connection timeout. The sandbox attempt failed DNS. Escalated SSH permission was
requested, but the turn was interrupted before the remote inspection ran to a
confirmed result. No remote scheduler or engine identity was established, and no
remote mutation or job launch is claimed. The user was informed that staged work
was not yet installed; the first installation was completed afterward.

## 8. Remaining work; criteria for a production release

1. Complete mapped conjugated-head preparation for remaining chemotypes, including
   distinct VHL exits, not just parent inhibitor transfer.
2. Resolve exact target-domain inventory, KLIFS provenance/alternate states, AF2
   confidence, mutations and missing binding-site loops. Preserve exclusions.
3. Verify the exact Ernie adapter and scheduler; remove silent scoring fallbacks,
   pin seeds/versions, validate local and remote recovery, and keep immutable attempts.
4. Run independent known-complex redocking and held-out cross-docking, plus controls
   and sampling convergence. Do not count native self-placement as predictive evidence.
5. Validate ligand charges/torsions/topologies and freeze the MD protocol, replicas,
   equilibration, checkpoint handling, convergence and uncertainty analysis.
6. Freeze the exact model and schema; validate protein-family/compound generalization,
   missing features and calibration. Bridging is a structural descriptor, not a
   calibrated probability of degradation or ubiquitination.
7. Measure pilot runtime/storage and release a hashed snapshot only after these gates.

The project is not declared publication-ready merely because software tests pass.
No numerical production acceptance criteria were certified by this remediation.

## 9. Follow-up work after the first installation

Resolved SJF_delta using the current Tocris product route `sjfdelta_7267`. Its
SMILES produces exactly `C62H63F2N7O12S` and stereochemical InChIKey
`LYUUIIZKBMBNQV-CHOFKVETSA-N`. This supersedes the earlier unresolved status;
the initial unresolved record is preserved in the first-installation backup.
The earlier catalog-number question no longer blocks identification of the
five staged chemotypes, although the full ten-compound scope is still unspecified.

Fresh 3LQ8/88Z and 4W9H/3JF sources were downloaded and reconstructed from their
CCD definitions. The foretinib shared core is an exact 37-heavy-atom substructure
in both the deposited parent and Tocris SJF8240/SJF_delta. The unmethylated VHL
head is cropped from 4W9H by removing its three-heavy-atom acetyl cap, retaining
the 30-heavy-atom amine head. This matches SJF8240 and SJF1521 stereochemically,
and correctly does not match SJF_delta's different VHL ligand. No new coordinates
are invented for these fragments. SJF8240 placement uses the MET/foretinib domain
reference and the separately bound VBC reference; their heads are not assumed to
share a global frame before protein docking.

A new SJF8240–MET self-placement candidate passed with pocket RMSD about 1.1e-14 A.
Both native controls have only BENCHMARK_RELEASE_PENDING in preflight. Neither is
an independent redocking or cross-docking result.

Ran the corrected chain-aware join on the existing local KLIFS files without
changing their manifest: 496 filename-parsed files, 242 UNAMBIGUOUS_CHAIN,
254 AMBIGUOUS_METADATA. The remaining one of 497 PDB files was not silently
assigned metadata because it does not follow the parsed name convention.
The ambiguity table is saved as `audits/klifs_metadata_audit.tsv` after installation.
The 254 ambiguous records need exact KLIFS IDs or new provenance-backed acquisition.

Additional tests exercise chain-specific metadata joins, ambiguous-state rejection,
exact-ID disambiguation, tandem-domain preservation, sidecar hash verification,
foretinib/VHL bound-coordinate preservation and the distinct SJF_delta attachment.
Their final installed results and the second backup are recorded below by the
installer, rather than assuming that staged tests imply installed success.

## 10. Confirmed Ernie inspection and additional output safeguards

The later approved read-only SSH inspection succeeded. The earlier interruption
is no longer the final connectivity state. Ernie reports 96 CPUs; no scheduler
command was found in the checked PATH/common locations. The exact remote runner,
adapter, main, resume and LOCAL backend were copied into `engine_audit/ernie/`.
Four snapshots match separately obtained remote SHA256 hashes byte-for-byte;
LOCAL's hash is recorded from the retrieved content. No external engine file was
edited, and these snapshots are not imported by the repaired project workflow.

The remote adapter does have the feared -1.0 placeholder fallback. The remote
LOCAL backend also swallows failed child return codes and exceptions. The resume
script can substitute score zero after an exception and shift raw scores. Two
embedding paths disable chirality enforcement. These are now remote-confirmed
findings, not inferences from a local installation. Full details, versions, paths,
limitations and follow-up requirements are in `docs/ERNIE_ENGINE_AUDIT.md`.

The stage-report contract was strengthened because a surviving subset of PDBs
could otherwise look complete after hidden worker failures. Completion now
requires `engine_status.json`: explicit engine hashes and integer seed; raw score
role; each stage's attempted/succeeded/rejected/failed counts; arithmetic balance;
zero failed workers; linker-success/scored-pose count agreement. The exact current
Ernie engine cannot pass until it produces such accounting. The stage report is
hashed into the completion receipt and rechecked during extraction.

Added output molecular identity validation: select the full ligand explicitly;
reconstruct its PDB heavy-atom connectivity; assign bond orders from the expected
full Tocris graph; derive stereochemistry from the actual output coordinates;
compare stereochemical InChIKeys. Template assignment does not replace coordinates
or certify a pose whose chiral centers have inverted. A mirrored-ligand regression
case confirms this check catches identity errors that atom counts cannot detect.
Potential ambiguous bond-order matches are not used to choose a new pose; the
final graph/3D identity must match the supplier molecule or validation fails.

The inspected Ernie environment has Python 3.10.20, NumPy 2.2.6, SciPy 1.15.2,
RDKit 2025.9.5, pandas 2.3.3 and PyRosetta 2026.3+releasequarterly.5e498f1409.
No PyRosetta initialization or molecular calculation was performed for the version
check; distribution metadata were read. `/home` had about 1.4 T free, 81% used.
At 404,923 bytes per native-example complex, 2,860 x 10,000 such PDBs would occupy
about 11.58 decimal TB uncompressed, excluding all intermediates and MD. This is
illustrative arithmetic, not a prediction of accepted output count or runtime;
it demonstrates why a measured pilot and retention plan are mandatory.

42 staged regression cases passed after these additional safeguards. Final
installed test evidence is recorded in `audits/regression_tests.txt` and the
installation history. New cases cover failed-child accounting, missing stage
reports and mirrored output stereochemistry. No scientific simulation was launched.

## 11. Reproduction commands and artifacts

Run from the actual project directory after installation:

```bash
python3 -m pytest -q tests
python3 scripts/01_setup/00c_verify_supplier_registry.py
python3 scripts/01_setup/00d_prepare_bound_references.py --sources reference_sources
python3 scripts/01_setup/00e_prepare_gsk215_heads.py
python3 scripts/01_setup/00f_prepare_foretinib_heads.py
python3 scripts/01_setup/02_batch_prosettac.py --plan
python3 scripts/01_setup/01b_preflight.py --output audits/legacy_preflight.json
python3 scripts/01_setup/01b_preflight.py --inputs prepared_inputs --output audits/candidates_preflight.json
```

The preflight commands intentionally return 1 while production is blocked.
Candidate preparation is write-once: the supplied native/foretinib control tables
must be run in a fresh output snapshot if their candidate directories already
exist. Do not delete an existing attempt to simulate a resume. No command above
launches a docking or MD simulation or a remote job.

Exact code changes and original/source digests are in each backup's
`install_plan.json`. The original committed implementation remains available in
git and the pre-install backup. Supplier molecular SMILES and every source link
are in `config/molecules.json`; downloaded source checksums are in
`reference_sources/sources.json`; bound atom maps are in each reference's
`provenance.json` and `heads.json`; selected protein residue mappings and file
checksums are in each candidate's `preparation.json`; current preflight decisions
are in `audits/`; the raw Ernie evidence is under `engine_audit/ernie/`.

## 12. Modified conjugated-head modeling (candidate preparation, not a validated result)

Reason for this additional step: the remaining actual Tocris head graphs differ
from the parent inhibitor or VHL ligand in their crystal references. Transferring
a complete unrelated parent ligand would be chemically wrong; calling newly
constructed coordinates experimental would also be wrong.

`00g_model_modified_heads.py` now extracts the actual supplier head fragments at
explicit chemical attachment patterns and retains their atom-index/cut-bond maps.
SJF1521 uses its modified 34-heavy-atom kinase head; JNJ1013 its 29-heavy-atom
kinase head; SJF_delta its different 39-heavy-atom phenyl-linked VHL head.
The other head in each bundle comes from an exact bound fragment already described:
SJF1521 VHL from 4W9H, JNJ1013 VHL from 7PI4, SJF_delta foretinib from 3LQ8.
Every bundled head is checked against that compound's own full supplier graph.

For each modified head, RDKit MCS matches exact element/bond-order scaffolds with
ring-to-ring and complete-ring restrictions. The bound common scaffolds contain
27 atoms (SJF1521), 22 (JNJ1013) and 28 (SJF_delta). Complete-ring restrictions
were applied to both atom and bond matching after a development attempt exposed
an incomplete aromatic scaffold; no such invalid draft was installed.

Constrained embedding is repeated with integer seeds 1-12. Common-scaffold atoms
are restored to their deposited coordinates and fixed during UFF minimization
(maximum 2,000 iterations). Nonconvergence, nonfinite energy, scaffold movement,
changed 3D stereochemistry and minimum protein heavy contacts below 1.5 A cause
rejection. Among acceptable conformers, the lowest internal UFF energy is selected,
with seed order as a deterministic tie break. This energy is used only within the
same head; it excludes protein interactions and is not a binding affinity.
All 12 seeds passed the provisional checks for each of the three heads in this
local environment. That is 36 small local conformer/minimization calculations,
not docking recovery, MD simulation, a convergence study or an Ernie job.

Common-scaffold symmetry may yield multiple exact matches. The selected first
RDKit graph match, match counts, complete mapping, every seed, failed seeds,
chosen seed, internal energy, clashes, source hashes and review requirements are
recorded. Symmetry and novel-atom binding modes remain explicit scientific review
items; a chemistry match does not establish an accurate analog pose.

Modeled references are separated under `modeled_reference_candidates/` rather
than silently placed among experimental references. Their provenance sets
`experimental_coordinates_for_unmatched_atoms: false` and `production_ready: false`.
The placement command requires `--allow-modeled-candidates` to consume them;
otherwise it refuses these chemotypes. Each resulting preparation.json retains
`contains_modeled_head: true` and the source bundle hash. Preflight adds both
MODELED_HEAD_REVIEW_PENDING and BENCHMARK_RELEASE_PENDING.

Prepared three additional self-placement controls: SJF1521 on its EGFR reference,
JNJ1013 on its IRAK1 reference, and SJF_delta on the foretinib/MET template.
The SJF_delta–MET control is deliberately only a head-transfer preparation control;
it is not a claim that MET is its intended biological target or that this models
its p38delta ternary complex. Transfer to actual screen domains remains separate.
All five compound preparation controls pass mechanical chemistry/geometry checks;
none is released for production.

Additional regression cases verify the fixed bound scaffold, the supplier-compatible
3D stereochemistry, correct common-atom counts and mandatory modeled-head review.
The latest installed test result and candidate decisions are appended by the
installer and recorded under audits. All earlier counts in this document are
historical checkpoints, not claims that the earlier test count remains current.

## 13. Preparation invocation audit trails

The placement script now saves every invocation's results under a unique UUID in
`audits/preparation_attempts/`, including the exact input table hash, whether
modeled references were explicitly enabled, compound/target identity, and success
or rejection reason. `prepared_inputs/preparation_summary.json` inventories all
existing prepared controls rather than overwriting the inventory with only the
last filtered invocation. Each inventory row hashes its preparation metadata.
This prevents the MET controls for two different compounds from being confused.

The additional modified-head calculations are preproduction local preparation,
not an authorization to launch the Ernie campaign. All five compound controls
remain blocked; the three modeled cases have an additional explicit review gate.

## Installed validation receipt — 20261104T134519Z

Backup: `/home/ayushregmibagale/10PROTAC-MD/audit_baselines/before_remediation_20261104T134519Z`. All planned repairs installed and the suite passed:

```text
................................................                         [100%]
48 passed in 1.09s
```
No SSH, docking campaign or MD simulation was executed by this installer. Small local RDKit head-conformer calculations are documented above. Production release remains blocked.

## 14. Continued preparation and reference-project comparison

The user explicitly asked to consult `protac-md`, then clarified that its working
workflow should be made publishable rather than replaced. That project was read
only. `PROTAC_MD_COMPARISON.md` explains adopted settings and discrepancies;
`PROTAC_MD_REFERENCE_SNAPSHOT.json` hashes the inspected source files. The actual
MD code uses AMBER99SB-ILDN/GAFF2/TIP3P at 300 K, unlike its older paper methods.
The staged equilibration and seeded representative-pose sampling rationale are
retained. Old documentation is not evidence that a different force field ran.

An exploratory constrained-chemistry docking engine was developed in temporary
storage before this clarification. It is **not the production replacement**.
Its source and benchmarks are archived separately under experiments. The ordinary
batch launcher cannot invoke it. Existing protac-md results were not regenerated,
replaced or invalidated merely because this exploratory method was investigated.

## 15. Exact structure provenance and defined domains

`00h_resolve_klifs_provenance.py` reproduces the old selection rule only to nominate
an ID, then requires an identical atom/residue/coordinate/element fingerprint
from the ID-specific KLIFS response. All 496 conventionally named legacy files
were verified, resolving the former 254 ambiguous joins. Network timeouts were
retried through a validated cache, including LRRK1 14879 and LRRK2 15826. No
original coordinate file was overwritten. The additional `RSK3/2777_in.pdb`
uses a fallback filename and is outside that verified set. A request for ID 2777
returned an empty response; it remains excluded rather than assigned by guess.

Cached UniProt annotations and sequence, plus all 555 human KLIFS kinase metadata
entries, drive `00i_build_domain_inventory.py`. Request plans and file hashes are
recorded in annotation_download_provenance.json. A spurious legacy identifier
`Entry` was a table-header artifact, not a protein. The annotation download also
includes VHL, ELOB and ELOC for sequence-based protein repair.

Final defined inventory: 770 domains from 435 UniProt proteins (342 KLIFS and
428 AF2). Explicit exclusions cover absent/non-kinase annotations, sequence
mismatch, insufficient coverage, oversized domains and low-confidence AF2 domains.
Structures with modified residues are excluded pending an explicit chemical-state
model instead of silently dropping phosphoresidues. Domain boundaries come from
UniProt annotations, preserving tandem domains. Confidence and DFG source strata
are retained; no missing DFG state is guessed. TTN residue numbers exceed the PDB
four-column limit: an explicit domain-relative offset preserves the UniProt
correspondence without producing malformed PDB records.

`01c_prepare_domain_screen.py` tested 3,850 compound/domain/conformer combinations.
Sequence correspondence alone yielded 35 candidates and rejected 3,815. That
attempt is preserved. A separate CE structural-correspondence attempt retained
112 candidates (GSK215 8, SJF8240 5, SJF1521 49, JNJ1013 45, SJF_delta 5), with
3,738 explicit rejections. Identical pocket RMSD, coverage and clash requirements
were retained; thresholds were not loosened to admit more kinases. This limited
coverage is a methodological limitation requiring resolution/justified scope,
not a full-kinome success claim. A rejected transfer is not a nondegrader label.
Unique run names and manifest keys preserve domain and source-conformation IDs.

## 16. Compound-specific ligand parameters

`04_parameterize_ligands.py` generated five separate GAFF2/AM1-BCC sets using the
exact supplier stereochemical molecular graphs, explicit H, neutral catalogue
reference microstates, ETKDG/UFF starting geometries, ACPYPE and AmberTools.
Each calculation recorded its input, command, logs, name map, connectivity and
artifact hashes. No alternative charge method was substituted after failure.
Output atom counts are GSK215 129, SJF8240 144, SJF1521 137, JNJ1013 118 and
SJF_delta 147. All topology partial-charge sums match the formal charge, masses
are finite/positive, topology bonds match the mapped source graph, frcmod files
contain no unresolved revision markers, and the quantum-optimized coordinates
retain the supplier stereochemistry.

The first delta conformer embedding failed; bounded explicitly seeded retries
resolved preparation without changing the chemistry. SJF1521 finished its charge
calculation but initially failed name validation because AmberTools changed Cl65
to CL65 in the GROMACS topology. An explicit bijective chlorine-name alias was
recorded, the original mapping preserved, and connectivity plus optimized 3D
stereochemistry rechecked. Its expensive quantum calculation was not repeated.
Neutral catalogue microstates are reference parameter sets, **not evidence that
all five compounds are neutral in the intended assay**. Protonation/tautomer
assessment and sensitivity remain a production design requirement.

The initial isolated runtime had an older RDKit lacking SetCoordMap. Its failed
engine attempt was retained. A temporary Python 3.10 runtime then used the same
RDKit/NumPy/SciPy/pandas and PyRosetta package versions observed on Ernie. The
user's existing PyRosetta wheel was extracted locally; its binaries/database are
not redistributed in this project. Additional packages were installed only into
the isolated temporary runtime. PDBFixer source commit:
5e658a4fb8d2b90d65ca8408c845cbf104ce7099. Official Rosetta ligand-generator source
commit: eddc15bbac60cf91302d3faf76895bd6d215cc50, with source hashes/license headers
preserved in the experimental snapshot. This is separate from the baseline
adapter's legacy parameter generator, which still requires a compatibility audit.

## 17. Protein completion, mapped topology and bounded MD diagnostics

`05_complete_protein.py` aligns each chain to cached UniProt sequence, fills
internal gaps and missing atoms with seeded PDBFixer, preserves observed heavy
atoms, and validates peptide bonds and ligand clashes. It does not invent
unobserved terminal tails or use pdb2gmx -missing. Native GSK215 required internal
ELOC and FAK loops; every added residue is recorded in protein_complete.pdb.repair.json.
Observed heavy atoms moved by 0.0 Angstrom; closest protein/ligand heavy distance
was 2.45688 Angstrom. These added loops are modeled hypotheses.

`06_prepare_md.py` transfers ligand coordinates by a stereochemical chemical map,
places/minimizes H with heavy atoms fixed, and checks topology name order. It uses
AMBER99SB-ILDN/GAFF2/TIP3P, 1.5 nm dodecahedral padding, neutralization plus 0.15 M
NaCl and 300 K, with explicit seeds. Solvent is selected by name. Whole-system and
chain-specific/ligand heavy-atom groups replace numeric solvent and ambiguous
chain/residue-cutoff selections. Force-field files are copied into new systems
for portable, hashed downstream inputs. No blanket -maxwarn is used.

The staged plan is EM (force tolerance 500), 200 ps restrained NVT, 500 ps fully
restrained NPT, 500 ps ligand-release NPT, then configurable 20 ns production.
NVT creates velocities once; later stages continue checkpoints. C-rescale is used
for equilibration, Parrinello-Rahman for production. Solute-only compressed XTC
reduces storage. `07_run_md_stage.py` records immutable stage inputs/completions,
checks artifact hashes, and resumes with checkpoint checksum/append semantics.
Ordinary EM is rejected if it does not reach its force tolerance. No step-limit
exit can be misreported as converged minimization.

Actual local diagnostics: strict topology preprocessing passed for a 229,735-atom
solvated GSK215 system with a 10,253-atom solute. A 1,000-step minimization and
20-step (0.04 ps) restrained NVT ran successfully. These are **mechanical smoke
checks**, not completed equilibration. The short EM retained Fmax ~6,171 kJ/mol/nm,
above the scientific 500 target. The first NVT preprocessing found a missing
System index group. That failure is preserved in system/; system_v2/ explicitly
adds the group, and its EM/NVT checks passed without warning overrides or LINCS
warnings. Production durations were not executed.

`02_analyze_md_contacts.py` uses mapped heavy-atom groups and triclinic
minimum-image distances, checks trajectory/preparation receipts, distinguishes
Angstrom from nm, and requires explicit discarded time and block duration.
The solute-only smoke trajectory was analyzed using a matching reduced TPR,
whose parent hash and conversion command are recorded. Smoke output is labeled
DIAGNOSTIC_ONLY. Block variation is not presented as independent-replica
uncertainty or proof of convergence.

## 18. Exploratory docking evidence and why it is quarantined

The exploratory engine's native GSK215 control produced one genuine scored pose,
with raw ref2015 +124.813187 (no forced negative shift). A same-seed repetition
had identical atomic records and score. Raw PDB bytes differ because Rosetta
writes the absolute output path into energy-table comments. A small perturbed
control produced eight chemically valid poses with POI C-alpha RMSDs after VHL
alignment of approximately 1.29–6.39 Angstrom. This is a local integration test,
not an independent global recovery benchmark.

An independent global diagnostic sampled 300 chemical proposals, admitted one
protein orientation, and rejected its full-linker construction: terminal outcome
COMPLETE_NO_VALID_POSES. That is an unsuccessful sampling budget, not a biological
negative. It is a reason not to deploy the experimental engine. Geometry yield
now retains all attempted protein poses in its denominator. Conditional contacts
among valid head-constrained ternaries are not degradation probabilities, and
old trained ML weights are not automatically applicable to a changed ensemble.

## 19. Minimal repairs of the existing baseline adapter

`baseline_engine/` contains an isolated, reviewable patchset against the archived
Ernie adapter/LOCAL/resume sources and a separately identified local clustering
source. Original and patched hashes and a unified diff are retained. It keeps the
existing adapter's docking/conformer algorithms and ref2015 side-chain repacking.
Changes propagate failed child processes/repacking/scoring, remove fake -1/0
score fallbacks, seed Python/NumPy/RDKit/PyRosetta, enable chirality enforcement,
remove count-only resume acceptance, and allow finite positive raw scores during
clustering instead of shifting scientific score.sc. Worker default is bounded.
Seeds derive from the explicit base seed and task path relative to the declared
attempt root. The patchset has not overwritten the original projects or Ernie.

The baseline still needs its complete runtime dependency snapshot, an end-to-end
smoke/pilot using its exact ligand/virtual-atom conventions, checks of both bound
head geometries in final outputs, and a parameter-typing audit. Missing evidence
is not converted into a production-release certificate. The batch launcher
continues to block production and explicitly excludes the experimental engine.

## 20. Baseline repacking component test and ligand registration

After installation of the first publication-repair package, the baseline
adapter's actual repacking function was tested locally on the native GSK215
complex. It returned a finite unshifted score (+127.829 at recorded precision),
retained all protein heavy-atom identities, and retained the full supplier
molecular graph and 3D stereochemistry. No replacement docking algorithm was
involved. This test used the explicitly recorded official PTC params from the
native diagnostic, so it does not validate the older parameter generator.

Inspection also found that `_init_pyrosetta` returned immediately after its first
initialization, even when a later stage requested a different ligand parameter
file. With ignore_unrecognized_res true, that can allow a later residue type to
be omitted. The isolated baseline patch now uses strict unrecognized-residue
handling and explicitly registers each relaxation stage's nonstandard residue
set in that pose before loading its PDB. This follows the installed PyRosetta
API's generate_nonstandard_residue_set documentation/example. A second real test
initialized globally with PT0 first, then repacked a PTC-containing native complex
in the same process. It passed graph/stereo and protein-inventory checks. Inputs,
outputs, logs, genuine score and component validation receipts are retained under
benchmarks/baseline_repack_gsk215_v1 and v2. The complete baseline dependency and
virtual-atom pipeline test remains a separate requirement.

Derived domain/preparation paths were made relative to project_root for transfer.
Original pre-portability metadata is archived in audits/pre_portable_metadata.zip
with before/after hashes. Historical execution receipts retain their actual
absolute execution paths rather than being rewritten to pretend a different
location was used. Installation receipts and checksummed backups record all
changes to 10PROTAC-MD; neither protac-md nor Ernie was modified.

## 21. Faithful baseline dependency integration and real stage diagnostic

Read-only SSH retrieved the existing 344-KiB PRosettaC source tree and the three
adapter helpers (sdf_to_params.py, run_docking.py and constraint_generation.py).
The original tar archives, hashes and retrieval scope are retained in engine_audit/ernie.
Ernie's OpenBabel library reports 3.1.0. A pinned OpenBabel wheel 3.1.1.22 was added
to the temporary local validation runtime; it also reports library 3.1.0. Exact
local package/interpreter versions are in baseline_engine/environment.json.
No remote source, environment or simulation state was changed.

An isolated runtime is assembled under baseline_engine/runtime from the archived
sources plus explicit patches. Its file hashes are in runtime_manifest.json.
Each rebuild archives the previous runtime in audits/baseline_runtime_history.
Paths are supplied explicitly instead of falling through to the user's original
PRosettaC installation. Existing global/local docking and free-conformer sampling
remain the baseline algorithm. A separate diagnostic driver supplies the native
GSK215 global frame, one local docking result and two full-molecule conformers.
These deliberate budget reductions are labeled diagnostic and are not production
settings or an independent global docking benchmark.

The parameter audit confirmed that the old helper stripped explicit H, used its
own atom types/Gasteiger charges, and wrote only the first conformer despite a
multi-conformer call. The adapter now has a compatible official Rosetta parameter
helper, preserving explicit H, the supplied coordinate frames, an all-atom name
map and the supplied conformer library. The official tool and hashes/license
remain included. This repairs ligand parameterization inside the existing
workflow; it does not introduce a different protein docking search. Parameter
and chirality repairs can change accepted poses/scores, so old trained models
cannot be assumed calibrated to the repaired outputs.

Actual native-frame integration attempts, all preserved:

* v1: stopped at local docking because the original caller prepended its working
  directory to already-absolute returned parameter paths. The adapter return-path
  contract was corrected while retaining absolute internal tool paths.
* v2: reached constraint generation, where the old head reader counted explicit H
  as heavy atoms and failed. The subsequent cleanup tried to remove a virtual-atom
  file that had never been created. Head extraction now uses verified parameter
  name/index provenance rather than PDB serial order or H-inclusive counts.
* v3: reached full-PROTAC parameterization and repacking, exposing a ligand-rotamer
  atom-name mismatch. Rosetta's exact four-column PDB atom-name fields are now
  preserved, including spaces, in every conformer; stripped display names alone
  were insufficient for the rotamer-library reader.
* v4: generated/scored a real full PROTAC but the original parent returned zero
  while its clustering helper was missing. The original unchecked os.system call
  hid that failure. This attempt is explicitly marked FAILED_OUTPUT_STAGE.
* v5: completed initial head preparation/relaxation, anchor sampling, the specified
  native global frame, original local docking, free-conformer generation, official
  full-ligand parameters, repacking/scoring and clustering. Checked subprocess
  calls and the explicit original-pipeline clustering path replace the failing
  stage interface. The raw positive score remains unshifted. The result summary
  correctly says finite raw scores rather than energies below zero.

The v5 pose retains the correct molecular graph and 3D stereochemistry but has
protein-ligand minimum distances about 1.224 A to VHL and 0.493 A to the POI.
Thus a naive <5-A contact rule calls it bridging even though it has a severe
steric overlap. baseline_engine/audit_outputs.py records this as rejected in a
separate QC report, preserves raw scores/coordinates and emits no biological
label. The complete *mechanical* diagnostic therefore passed its stage interfaces
but produced zero QC-valid poses in that tiny budget. It is not a publication
release or evidence that the original large-ensemble sampling never finds a
valid pose. Appropriate quality filtering and a representative pilot remain
necessary before spending months on the campaign.

Concurrent score writes now use a file lock and retain one header and all raw
rows, at six-decimal precision. Added regression cases test concurrent writes and
name-mapped head extraction with reordered PDB records and explicit H. The latter
prevents a superficially successful format conversion from corrupting head
coordinates. All patch sources, runtime-build sources, failures and hashes are
retained for review and reproduction.

## 22. Binding-mode clarification requested by the user

GSK215 uses a VS-4718/PND-1186-derived FAK-binding head in the ATP pocket; its
reported highly cooperative ternary-complex formation is distinct from binding
an allosteric pocket. The actual GSK215/FAK/VHL structure 7PI4 remains the bound
reference, rather than a generic ATP-site pose. Sources checked:
https://www.rcsb.org/structure/7PI4 and the original study
https://pubmed.ncbi.nlm.nih.gov/34416073/ (Discovery and Characterisation of Highly
Cooperative FAK-Degrading PROTACs). This clarification did not change the molecule
or substitute a guessed ligand pose.

## 23. Larger bounded baseline diagnostics, task seeds and output contracts

The native-frame diagnostic driver now accepts explicit bounded local-pose and
conformer budgets; production defaults are unchanged. v7 used three local outputs
and the existing baseline's full 100-conformer budget. All three outputs retained
supplier graph/stereochemistry. At the historical coarse 1.0-A overlap screen,
one passed and two were rejected. These numbers are diagnostic, not a recovery
rate or a biological performance claim.

Inspection of task seeding found that child processes shared a working directory.
Seeding from that directory alone gave different tasks the same initial seed.
LOCAL now assigns a deterministic task key from normalized script/arguments,
excluding absolute interpreter/install/run prefixes. The adapter derives stage
seeds from the explicit base seed, relative working directory and that task key.
Thus worker scheduling does not define the seed, while different local/constraint
tasks receive different streams. Each worker's full stdout/stderr and normalized
command, task key and exit code are retained in LOCAL_logs. Regression checks
verify same-task repeatability and changed streams for a changed task key.

v8 reran the three-local/100-conformer diagnostic with this final task-seed
handling. All original stages completed, and all three final graphs/stereochemistry
passed. Two passed the historical 1.0-A gross-overlap screen and one failed. The
scores remain genuine, positive and unshifted. There is no claim that a gross-clash
pass proves a correct binding pose or favorable energy.

The head-preparation path has always rejected contact distances below 1.5 A.
New selection contracts now explicitly carry that same minimum heavy-atom
contact cutoff, units and rationale, rather than relying on an implicit weaker
1.0-A output default. Historical executed selections/receipts remain unchanged
and retain their historical default. A separately named 1.5-A sensitivity report
for v8 rejects all three poses (the two historical coarse passes have minimum
contacts ~1.286 and ~1.478 A). This is fully recorded, not hidden by overwriting
the earlier report or redefining it as success. New prepared controls/candidates
have updated selection hashes and preparation-summary hashes; prior versions
are archived in audits/before_explicit_qc_selection.zip. There is no automatic
biological label when the valid subset is empty.

After the remote dependency snapshot became available, the current clustering
patch was rebuilt against the exact Ernie original, replacing the initial local
variant. The shared changes preserve the original clustering mathematics while
fixing absolute source/summary paths, keeping finite positive raw scores and
raising for inconsistent atom selections. A separate v6 clustering check on
identical archived input poses passed. Original local/remote hashes and all
intermediate runtime snapshots remain available; provenance does not conflate
those source versions.

The final prepared-screen preflight finds no mechanical input errors: its only
blocks are BENCHMARK_RELEASE_PENDING for all 112 candidates and the additional
MODELED_HEAD_REVIEW_PENDING for 99. These gates are deliberately still present.
Software checks passing cannot replace representative geometry-valid sampling,
physical equilibration, binding-mode assessment or independent experimental
validation. No such missing evidence is labeled complete.

## 24. Measured storage projection and lossless structural archiving

`02a_resource_plan.py` measures actual baseline output/intermediate bytes and
requires the pilot's recorded conformer budget to match the projected budget.
For 112 candidate pairs, 200 global by 50 local poses and 100 conformers per local
pose, v8's measured files project to about **7.24 TB raw**, or **1.40 TB if every
intermediate is individually gzip compressed**. Even gzip lacks a 30% reserve
against the earlier observed ~1.4 TB free on Ernie, before MD, retries and backups.
This is a small GSK215-based projection, not a disk reservation or a measured
runtime/memory guarantee. Protein sizes and other chemotypes can change it.

`03_archive_pose_ensemble.py` implements a lossless, portable alternative for
scored, local-docking and global-docking PDBs: chunked NumPy int32 coordinates at
exact PDB precision (0.001 A), compressed JSON metadata, shared atom-field
templates, sparse field changes and original source SHA256. Histidine hydrogen
names can change after packing, so sparse atom-field changes are preserved rather
than assuming a single identical topology. Negative zero, atom order, non-atom
records, all coordinate fields and final-newline state are retained. No pickle
payloads are used. It does **not delete source files**.

The initial prototype duplicated entire templates for small atom-field changes
and was larger than its seven input PDBs. That attempt and restore receipt remain
preserved. Sparse templates plus compressed manifests fixed the overhead. The
final v8 archive holds seven source PDBs: 6,231,904 original bytes became
1,256,982 archive bytes. All seven were restored and their complete SHA256 hashes
matched the originals. Regression tests cover changed non-coordinate fields,
negative-zero coordinates, missing final newline and corrupted-chunk rejection.

Retaining the measured lossless pose archives plus compressed scores, parameters,
commands, code/environment/seed provenance and logs projects to ~477 GB for the
112-pair example. This assumes a separately validated policy for discarding
reproducible temporary conformer/intermediate files after completed verification;
that cleanup policy was **not executed**, and raw diagnostics are all retained.
It excludes extra MD trajectories, retry budgets and archive-creation overlap.
The archive utility demonstrates a safe structural retention building block; it
does not itself certify that the full campaign fits or launch any jobs.

## 25. Open protonation-condition question

All five verified Amber/GAFF2 sets are explicitly labeled neutral supplier-reference
microstates. The user was asked which assay/physiological pH should govern the
remaining protonation assessment (pH 7.4 assumption, provided assay pH, or 6.5–8.0
sensitivity). An unanswered question is not treated as an instruction to silently
choose a different chemical state. Independent code/documentation/installation
work continued while that clarification was pending.

## User storage preference — full output retention

The user explicitly accepted the approximately 7.24-TB docking-storage estimate.
Full raw output retention is therefore the planning default; the lossless archive
utility remains an optional data-packaging tool, not a required reduction in
scientific outputs. The estimate covers the measured docking/intermediate model,
not additional MD, retries and backups. Acceptance of the storage budget does not
establish that the currently inspected Ernie filesystem has that free capacity;
mounted storage is being checked read-only before a campaign path is selected.
No files were deleted and no simulation campaign was launched.

## 26. Both-head linker correction within the original baseline

Native-frame tests of the inherited free-conformer rigid fit showed full ligand
heads drifting approximately 6–16 Å despite apparent bridging contacts. A rigid
fit cannot satisfy two bound heads if the free conformer's internal geometry is
wrong. The original global/local docking stages were preserved. The linker step
now uses stereochemical, disjoint head mappings, two-head ETKDG coordinate
constraints and fixed-head UFF refinement. Both full-head initial RMSD/maximum
errors and post-refinement coordinate shifts are checked. Impossible geometries
are explicit sampling rejections; no free-conformer fallback is allowed.

The direct native linker control generated two feasible conformers; actual
PyRosetta repacking produced a finite raw score 188.908562 with minimum heavy
contacts 2.4569/2.6345 Å and zero head-core displacement. The integrated v9 native
frame/local docking diagnostic accepted one pose, raw score 143.080191, from
three local poses with two conformers per feasible local geometry. Two local
geometries were infeasible. Old failures v1–v8 and their source snapshots remain
preserved. A native-frame control is not independent global redocking recovery.

## 27. All-compound chemical components and receipt-compatible producer

All five full ligand graphs passed official Rosetta parameter generation,
PyRosetta loading and stereochemical verification: GSK215 129 atoms, SJF8240 144,
SJF1521 137, JNJ1013 118, SJFδ 147, including explicit hydrogens. This validates
chemical components, not their binding modes. Atom naming retains the exact
four-column PDB fields, and each conformer is restored to the supplied frame.

baseline_engine/run.py wraps the original stages in the receipt contract. Raw
engine output remains in engine_raw. Canonical output contains only separately
QC-accepted poses with the same unshifted raw scores. Global sampling, local
docking, linker modeling, scoring and pose-QC counts are separately recorded.
QC checks supplier graph/stereo, gross overlap and docked head-core displacement;
these are necessary geometry checks, not proof of physical binding.
The v10 native-frame run completed with one accepted pose. The v11 one-trial
global run completed with no valid poses and zero downstream attempts, recorded
as COMPLETE_NO_VALID_POSES without inventing a biological negative.

## 28. Durable bounded controller, interruption and resume

The functional foreground controller now requires explicit Python/seed, enforces
pilot budgets, hashes source/input/protocol identities, uses fresh UUID attempt
directories, locks manifests and preserves attempt transitions. Completed runs
are skipped only after receipt revalidation for the same key; failed retries
retain earlier attempts. Global sampling records actual trials, including zero
successes. Full mode requires versioned evidence, exact source/input hashes and
a storage budget; no release file was fabricated. Timeout/interrupt handling kills
the whole process group and writes terminal receipts. SIGTERM/SIGHUP are handled
as interruptions. Regression coverage includes child termination and concurrent
score writes. No scheduler or SSH campaign is launched by this controller.

## 29. Explicit protonation assessment and charged GSK215 parameters

With no assay pH supplied, the stated provisional assumption is pH 7.4 with a
6.5–8.0 sensitivity assessment. Dimorphite-DL 2.1.0 retains mapped heavy atoms and
requires exact uncharged supplier-parent stereochemistry. At pH 7.4 it enumerated
two GSK215 states (neutral/+1) and a neutral state for each other compound. The
wider conservative interval enumerated eight GSK215 states and two SJF1521
states; it does not assign populations or establish that every state is populated.
Tautomer coverage and protein-environment pKa shifts are not claimed complete.
The full tool/data/source hashes are retained in the assessment.

The GSK215 +1 state was separately parameterized using GAFF2/AM1-BCC: 130 atoms,
net charge +1, correct optimized-output stereochemistry. The five neutral sets
remain independent reference states. Parameter validation now binds the exact
atom map and parameter-state InChIKey. The MD builder requires the selected state
explicitly and propagates these hashes; no state was silently made the primary
production choice. The charged system passed actual Gromacs topology generation.
A force-field include path issue found by this build was corrected using an exact
unique include match; the failed first build remains preserved.

## 30. MD execution details and limits

The MD stage runner now records auto/cpu/gpu device choice and optional GPU ID.
Existing CPU diagnostic receipts remain historical facts, not rewritten runs.
Nonempty preparation directories are rejected to preserve failures. Neutral
GSK215 EM (1,000 steps) and NVT (0.04 ps) passed mechanical checks, but EM force
remained above convergence and these are not production equilibration evidence.
No charged dynamics, full replicas, scientific trajectory or convergence claim
was generated. Periodic mapped-contact analysis was verified only on the tiny
neutral diagnostic trajectory.

## 31. SJFδ known-target coverage correction

The 112-candidate CE transfer screen excluded the established SJFδ MAPK13 target.
This is a structural-transfer exclusion, not biological inactivity. The original
Smith et al. study (doi:10.1038/s41467-018-08027-7; PMC6328587) describes foretinib
from EPHA2 5IA4 overlaid on p38δ 4EYJ, followed by OPLS3/Jaguar treatment.
5IA4 is not a p38δ complex and 4EYJ does not contain experimental foretinib.

00j_prepare_delta_p38_control.py reproduces an inspectable template-pair overlay
using official KLIFS 7189/2949 pocket residue correspondences. It fits 83/85
observed Cα positions, pocket RMSD 1.607829 Å and minimum head/protein heavy
contact 1.794427 Å, passing the existing 2.0/1.5 Å limits without relaxing them.
The exact supplier foretinib subgraph and separately modeled VHL head are
preserved. Source structures, CCD, pocket maps and original paper XML are hashed.
The new MAPK13 preparation is supplementary; the 112-pair screen and all 3,738
exclusions remain unchanged. Its evidence explicitly labels both modeling steps,
and does not claim reproduction of OPLS3/Jaguar torsional correction or an
experimental binding pose. Limited broad-kinome transfer coverage remains a
scientific design issue to resolve before general kinome-selectivity claims.

## 32. GitHub checkpoints and data scope

The user requested ongoing GitHub commits. Source, methods, compact verified
reference/control inputs, ligand parameters and selected small validation receipts
are prepared as reviewable commits. Full diagnostics, original archives, legacy
inputs and the structure inventory remain locally retained outside Git. Rosetta
helper hashes and an official-source fetcher are published instead of republishing
separately licensed helper source. PyRosetta binaries/databases are excluded.
Git exclusions do not delete scientific data. The user-approved full raw retention
preference remains in effect. Software checks pass 68 regression tests after the
MAPK13 control correction; installed validation is recorded separately by timestamp.

## 33. Explicit pocket-correspondence comparison

346 official KLIFS 85-position maps were retrieved with four bounded workers,
validated and cached with URL/time/SHA256 receipts. Experimental residue numbers
were mapped through the verified source/UniProt sequence correspondence into the
explicit domains. AF2 domains inherit positions only from consistent experimental
maps of the same UniProt domain; conflicting maps are rejected rather than voted
on. 556/770 domains have explicit evidence; 214 remain unavailable. The latter
are data-coverage limitations, not negative binding labels.

The independent comparison uses explicit pocket least-squares alignment plus an
independent domain-core check. It retains >=100 core residues, >=60% core,
<=2.5 Å core RMSD, >=80% head-neighborhood coverage, <=2 Å pocket RMSD and the
1.5–4.5 Å minimum head-contact interval. Both the complete mapped pocket and the
head neighborhood must pass the RMSD threshold; correspondence must be bijective.
The earlier 112 CE candidates and all their rejections remain unchanged.

The comparison prepares 161 candidates: GSK215 34, SJF8240 12, SJF1521 42,
JNJ1013 63 and SJFδ 10. It records 2,619 geometry rejections and 1,070 unavailable
compound/domain correspondences. All five known-target controls are represented,
including MAPK13/SJFδ. Only 41 candidates overlap the earlier screen, so the two
methods are not interchangeable and are not silently pooled into a larger study.
All 161 pass mechanical preflight; the remaining blocks are scientific release
and, for 115, modeled-head review. The alternate screen is comparison-only and
has not changed the production controller default or the user's 112-pair budget.

## 34. Global baseline pilot and accurate sampling counters

A bounded GSK215 global pilot (three requested global poses, three local poses
each, two conformers, at most 20,000 global trials, explicit seed and 600-second
limit) completed. The authoritative ledger records 89 random trials, three global
poses and nine successful local dockings. None of the nine local geometries
passed fixed-head linker generation, so no fabricated scores or poses were
created. Terminal status is COMPLETE_NO_VALID_POSES. This small sample cannot
establish inactivity, global recovery accuracy or production feasibility.

The older human-readable sampler line reported 90 trials because it counted an
extra loop entry that immediately broke. The ledger already counted the actual
89 samples. The runtime builder now makes the printed count use that same actual
counter; the historical pilot log and receipt are unchanged and preserved.

## 35. MD early-stop, TPR and preparation binding

A completion-marker-only check could accept a graceful Gromacs early stop. The
runner now verifies the final logged step against finite nsteps/init-step, and
requires a positive force-convergence statement plus finite Fmax below emtol for
scientific EM. Smoke runs remain explicitly diagnostic. TPR reuse now requires a
hash-bound generation receipt for the exact inputs/protocol; missing, modified or
unbound TPRs fail. Parent-stage logs, output hashes and smoke-EM dependencies are
validated before progression. Topology, local force-field files, index, initial
system and stage MDP are verified against the preparation receipt.

The topology builder's validation TPR is named preflight_em.tpr, leaving stage
TPRs to the stage runner and their own generation receipts. A diagnostic cloning
attempt found that the older neutral system used Gromacs' external force-field
lookup and had no local force-field directory. That attempt's 20-step output is
preserved as a historical mechanical test, not promoted as complete force-field
provenance. The current runner refuses such older systems. A fresh neutral system
was built using the current portable force-field preparation for the new checks.
An initial diagnostic command named a nonexistent pose; it failed before writing
system files and was corrected to the actual combined_0_0_0001.pdb path.

## 36. Requested Nature standard

The user explicitly requested Nature publication quality. NATURE_READINESS.md
links current Nature reporting/code/data guidance and the Communications Biology
MD reliability perspective, separates policy from project-specific scientific
choices, and maps actual evidence to outstanding requirements. Software tests and
native-frame integration are not represented as independent biological validation
or a guarantee of journal acceptance. Relevant experiments, sampling/convergence
and defensible claim scope remain substantive requirements.

The new portable v4 system passed actual 1,000-step EM and 20-step NVT mechanical
diagnostics with TPR generation, prepared force-field and actual final-step checks.
The diagnostic EM has not met emtol and is labeled SMOKE_COMPLETE only. A full
50,000-step-bounded EM with the actual Fmax<500 requirement is running separately;
a five-compound bounded global pilot also remains in progress. Active raw
benchmarks are deferred from this installation batch to avoid copying changing
files. Static smoke receipts and the comparison evidence are included now.

The MD runner also takes a nonblocking lock across the entire system directory,
preventing concurrent stages from overwriting shared outputs. Independent replicas
use separate systems. Regression coverage exercises lock contention. The new
software suite contains 76 tests. The already-running full EM diagnostic predates
this lock addition; no competing stage is started in its directory.

## 37. Quantitative native-interface recovery analysis

04_native_recovery.py aligns the POI Cα coordinates and measures mobile VHL Cα
RMSD and native protein residue-contact recovery at an explicit heavy-atom cutoff.
Compatible chain/residue/insertion/name coverage is required and duplicate Cα
identities are rejected. These are project metrics, not DockQ/CAPRI classes or
binding energies. A translated whole complex gives the same metrics; a VHL-only
30 Å decoy is rejected by the geometric comparison in regression coverage.

The accepted v10 native-frame/local-docking control has POI Cα RMSD 0.000508 Å,
VHL Cα RMSD 2.209 Å after POI fit, and recovers 12/16 native contact pairs (0.75)
at 5 Å. All reference Cα residues are mapped. This quantifies local native-frame
consistency, not independent global recovery. Raw v13 local protein geometries
are also measured separately; they are not full-ligand QC-accepted poses. The
software suite now contains 77 tests. No acceptance threshold or global accuracy
claim was invented after observing these values.

## 38. Explicit internal embedding budgets and performance repair

The five-compound v14 pilot completed for GSK215 and JNJ1013 without valid poses,
produced one QC-valid SJF1521 pose, then timed out on SJF8240. The timeout correctly
killed the process group and remained FAILED. RDKit 2025.9.5 uses a zero default
maxIterations to permit ten times the all-atom count of internal embedding retries
per call. Twenty outer seeded attempts multiplied this hidden cost. The installed
Python binding has no timeout property even though newer RDKit documentation does.

The repaired search explicitly uses 20 internal retries per outer attempt and
records failure counters. Exhaustion is not biological inactivity. Bounds,
stereochemistry, both-head RMSD and UFF convergence were not relaxed. Producer
settings pin the outer budget, resume worker count and parent task-key state so
inherited environment variables cannot silently change a cached protocol. The
SJF8240 retry completed within the small pilot limit. These budgets remain sampling
choices, not evidence of exhaustive conformational coverage.

## 39. Invariant ligand-rotamer orientation frame

Further multi-conformer controls found head displacement after Rosetta packing,
including approximately 0.636/1.605 Å on one SJF1521 pose and much larger shifts
on another. Input Cartesian conformers had preserved both heads. Inspection of
the actual PyRosetta residue type and pinned official ResidueType.cc showed that
orientation uses NBR_ATOM and two bonded heavy neighbors. The old center was a
linker attachment atom; a movable neighbor rotated the stored conformer.

The neighbor center is now a ring atom wholly inside the fixed VHL head with its
heavy neighbors fixed too. Parameter generation verifies frame invariance across
the complete conformer library, rather than disabling ligand conformer selection.
An actual retained failing SJF1521 ensemble was reparameterized and repacked.
The v19 control gives zero displacement for both head contact cores and minimum
contacts 2.435/2.818 Å without gross overlap. The first v18 harness omitted the
H-added head SDFs required by provenance matching; that harness failure and its
scored artifact are retained, and the missing provenance inputs were copied in
v19 without altering the original ensemble.

## 40. Frozen execution packages and input-drift checks

Staged source was edited while two diagnostic pilots were running. Those runs
are retained as diagnostics, not claimed as frozen-version validation. This
exposed a general reproducibility gap in workers reading the live working tree.
The controller now creates an atomic, verified campaign source/config package and
copies it into each attempt. Both recorder and producer execute from that package.
Later source edits cannot change workers or later pairs of the same campaign.
Python, JSON, XML, flags and parameter configuration hashes are included; no
PyRosetta binaries/database are copied into the package.

Input hashes are checked before and after copying each pair. Completed source
snapshots and their manifest are rechecked during receipt verification and skip/
resume decisions. Tests mutate the working tree between two simulated recorder
calls and confirm both use the original campaign code; separate tests catch
snapshot tampering and changed inputs. The real v21 controller completed a bounded
search and an identical invocation printed Already validated without rerunning.
Exact public source for the v20 five-compound pilot is archived in
validation_evidence/execution_sources/; separately licensed helper files are omitted
from Git and recoverable by their pinned official-source hashes.

## 41. Frozen five-compound end-to-end control

The v20 run uses explicit source snapshots, three global poses, three local poses
per global pose, two conformers, one worker, at most 20,000 global trials and a
600-second limit per compound. All five completed without execution errors.
GSK215, JNJ1013 and SJF8240 each produced no valid pose in that small search.
SJF1521 scored three feasible local geometries, retaining two after QC; SJFδ
scored three and retained all three. There were 45 local protein poses in total.
No-valid-pose results remain sampling outcomes, not inactive-target labels.
Independent receipt verification checked code/input/output hashes. This control
validates multi-chemotype execution, not general predictive selectivity.

## 42. Completed staged local MD validation

The current portable neutral GSK215 system converged to maximum force
472.168 kJ/mol/nm (threshold 500); Gromacs reports 4,222 EM steps. It completed
100,000 NVT steps (200 ps), 250,000 restrained NPT steps (500 ps), and 250,000
unrestrained NPT steps (500 ps), with actual final-step and TPR/input checks.
The initial GPU attempt failed because sandbox device access reported zero GPUs;
its log is retained. The authorized retry used the idle local RTX 5090 and the
existing CUDA 12.9/sm120-capable Gromacs build. No Ernie simulation was run.

The last 250 ps of unrestrained NPT has descriptive saved-sample means 299.989 K,
8.964 bar and 1004.720 kg/m^3 density. Gromacs accumulator statistics over the full
stage are retained separately and differ from sampled-window statistics. No
confidence interval or independent-replica count is fabricated. The 26 retained
solute-trajectory frames maintain contacts within 5 Å of both proteins; this is
short native-control evidence, not a binding free energy or convergence result.

05_md_thermodynamics.py binds energy extraction to stage receipts and records its
window/block definitions. 06_reduce_md_topology.py verifies the parent TPR,
preparation and exact ordered SOLUTE index, converts by group name and checks the
atom count. Equilibration contacts can now be analyzed explicitly as diagnostic
while production-only scientific summaries remain gated. Source/runtime metadata
and failed intermediate diagnostic setup attempts remain retained.

## 43. Real 100-conformer-budget controls and storage calibration

The v22 component tests regenerate 100 requested conformers for known feasible
SJF1521 and SJFδ geometries, then run real parameter generation and Rosetta packing.
They generated 74 and 97 valid conformers respectively; both preserve the two
heads exactly and pass gross-clash checks. Local measured component cost is about
44 seconds and about 1 GiB process peak RSS. These are local component measurements,
not aggregate parallel memory or promises about Ernie performance.

The frozen v23 native-frame GSK215 workflow requested 100 conformers per local
geometry, generated 100 in its feasible geometry, rejected two other geometries
and retained one scored/QC-valid pose. All original stages except independent
global search were exercised at that per-local conformer budget. The storage
planner was repaired to compare requested budgets using sampling ledgers, while
checking actual generated counts against parameter provenance; it no longer
incorrectly rejects a legitimate 74-of-100 outcome.

Historical upstream code removed empty conformer/reference files. The planner
handles such zero-output ledgers without inventing missing structures. The
current isolated runtime retains those files too, honoring full raw retention;
the v24 real control verifies all three zero-conformer markers and references.
This changes no molecular geometry or scoring. Inherited files and earlier runs
are never rewritten to pretend those empty files were present historically.

The new small native pilot projects ~3.37 TB raw for the same 112-pair,
200-global/50-local/100-requested-conformer example. The earlier ~7.24 TB was based
on the older free-conformer output pattern. Different feasible fractions and
protein/chemotype sizes change these estimates, so the accepted ~7.24-TB planning
allowance is not silently reduced. Neither estimate includes MD, retries, backups
or a guarantee of capacity; both exceed currently observed /home free space.

## 44. Full-complex and intraligand clash coverage

Clash QC now covers ligand versus all protein ATOM records and the full E3 complex
versus POI, including elongins outside the VHL binding chain. It also checks
intraligand heavy-atom distances for pairs separated by more than two bonds.
These are gross-overlap exclusions at the existing 1.5 Å preparation threshold,
not universal journal rules or affinity estimates. Tests introduce ligand/elongin
and E3/POI overlaps that a VHL-binding-chain-only check would miss; covalent and
angle-related ligand pairs are excluded from the self-overlap test. All five
previously accepted frozen-pilot poses pass the expanded checks. Their E3/POI
minimum heavy distances are approximately 2.81–3.04 Å.

## 45. Actual-runtime compatibility testing

The suite first passed in Python 3.13/Biopython 1.85. Running it in the actual
Python 3.10 docking environment exposed two CE-interface failures: Biopython 1.84
returns paired index lists; 1.85 returns named results with .path/.length. The
adapter now normalizes both forms without changing correspondence or thresholds.
A regression test verifies identical indices for both APIs. All 91 tests then
passed in each environment. Pytest 8.4.2 and its pure test dependencies were placed
under /tmp for the 3.10 run; no user environment was modified. The failing log,
passing summaries and exact observed runtime versions are retained.

## 46. Scope, Ernie readiness and independent recovery benchmark

Fresh read-only Ernie probes show 96 CPU affinity slots and about 1.50 TB free in
/home. No Gromacs executable was discovered in PATH or the few standard locations
checked; this is not an exhaustive installation search. No scheduler or quota tool
was found in PATH. Limited share-directory inspection did not identify an allocated
personal writable directory. No remote files or jobs were created. The user was
asked for the allocated output path and whether to retain the five verified
compounds or supply the remaining Tocris catalogs for the original ten-compound
aim. Elapsed time is not treated as an answer to either required study input.

The original ML/reference project was inspected read-only. Its small-panel LOPO
results and MGTbind structural labels do not establish cross-compound transfer
for this new generator. ML_TRANSFER_SCOPE.md records the required schema,
independent-unit, label/assay and holdout checks. No unknown labels were recoded,
models retrained, old results modified or new probabilities fabricated.

A separate frozen GSK215 global-recovery diagnostic (v26) was started locally:
50 global poses, five local poses each, two conformers per local pose, two workers,
20,000 maximum global trials and a one-hour timeout. Its necessary QC-valid-pose
check and exploratory native-interface criterion (VHL Cα RMSD<=5 Å after POI fit
and >=0.5 native contact recovery at 5 Å) were recorded before execution. These
are project diagnostics, not DockQ/CAPRI classes or proof of general accuracy.
Native-frame controls are not substituted for independent global recovery.
The active benchmark is excluded from installation until its files stop changing.

## 47. Final local validation and preservation details

Full-complex clash coverage and explicit intraligand nonbonded checks were
rechecked against the completed five-compound pilot; all five previously accepted
poses still pass. A final frozen native-frame control exercises the latest
producer/QC. The actual-runtime CE compatibility fix raises the suite to 91 tests,
passing under both Python 3.10/Biopython 1.84 and Python 3.13/Biopython 1.85.
Exact versions, the initial two-failure log and passing summaries are preserved.

The resource planner additionally handles historical zero-conformer SDFs removed
by upstream cleanup without pretending nonempty missing ensembles are valid.
Current runtime cleanup is disabled for those empty markers and reference files,
and a real all-rejected local control verified their retention. Existing
quarantined experiment copies are frozen rather than receiving later helper
updates through the installer. Completed raw diagnostics and source packages are
installed into the project with drift checks/backups; the active global-recovery
benchmark is deferred. No simulation was submitted to Ernie.

CONTROL_ROLES.json makes explicit that SJFδ/MET is a transfer-geometry control,
not an asserted degradation-positive target. The MAPK13 published-template model
is the separate known-target preparation. No completed input/receipt was rewritten
to change its historical control role. Remaining compound count, allocated storage,
independent recovery and scientific-validation decisions are not silently resolved.


## 48. Installed and Git-only checkpoint validation

**What:** Installed the completed repair set and stopped benchmark artifacts with the active global-recovery benchmark explicitly excluded. Published checkpoint `0eae588` to GitHub master.

**How:** The installer checked source and destination hashes, backed up overwritten files under `audit_baselines/before_publication_repairs_20261105T182837Z`, atomically installed 3,810 changed/new files, and verified installed hashes. The installed regression suite passed all 91 tests. A fresh Git-index export in `/tmp/10protac_git_checkpoint_v3` also passed all 91 tests, and all 150 packaged evidence checksums matched. No original protac-md files or Ernie simulation outputs were changed.

**Why:** Tests in a staging tree alone can conceal untracked dependencies or missing public assets. A clean Git-only export checks that the published checkpoint includes the files needed by the regression suite. Licensed external docking tools remain separately required and explicitly documented.

The final frozen native-frame diagnostic v25 completed with one chemically/QC-valid pose. Receptor-aligned VHL Cα RMSD was 2.209166 Å and 12/16 native residue contact pairs were recovered (0.75 at 5 Å), with full compatible Cα coverage. This initializes from the experimental protein arrangement and is not independent global recovery. The separate v26 global-search diagnostic is still active; it has not been assigned a successful or failed recovery outcome.


## 49. Completed global recovery diagnostic and retained the missed criterion

**What:** The v26 global diagnostic completed, but did not meet the joint recovery criterion recorded before execution. No thresholds, seeds, outputs or failed searches were retrospectively changed or removed.

**How:** Independently revalidated the completion receipt, exact inputs, frozen execution sources, scores, final pose hashes, chemistry and full-complex geometry. All 50 accepted global arrangements (1,071 trials), 250 local poses and 250 linker ledgers reconcile. One local geometry yielded two conformers, one scored pose and one QC-valid pose. Receptor-aligned VHL Cα RMSD is 6.259669 Å; 8/16 native contacts are recovered at 5 Å. The criterion required both RMSD ≤5 Å and fraction ≥0.50, so this is a miss. The real positive raw score 453.929638 is retained without shifting or sentinel conversion. Both final bound-head core RMSDs remain zero.

Analyzed every raw local protein pose against the same reference and selections: none of the 250 met both thresholds. The closest by VHL Cα RMSD is 4.168323 Å with contact recovery 0.375; it generated no linker conformers. Raw protein contact fractions may differ from final packed poses because side chains are repacked. Thus linker rejection alone cannot explain the failure, and protein geometry alone must not be counted as chemically validated recovery.

**Why:** Successful software execution and plausible chemical geometry do not establish recovery accuracy. Keep this missed benchmark visible and withhold campaign release. It is not evidence of biological nondegradation, nor proof that the sampler can never recover the native structure.

**Next bounded diagnostic:** Prespecified adaptive follow-up v27 uses seed 155921 and exactly the same frozen source, input, 50 global ×5 local ×2 conformer budget, 20,000-trial maximum, two workers and 3,600-second cap. It was initiated after reviewing v26, so it is exploratory rather than confirmatory. Both outcomes will be retained; no Ernie job is submitted.
