# Remediation status

**Identified local input/software faults have been repaired and validated. The
scientific campaign release remains pending. No Ernie simulation was launched.**
The original protac-md project was read-only. Repairs are in 10PROTAC-MD; its
working global/local docking sequence is retained. The exploratory replacement
engine is quarantined, and existing archive copies are not overwritten by updates.

| Area | Current evidence |
|---|---|
| Purchased chemistry | Five exact Tocris stereochemical graphs, compound-specific heads and anchors verified |
| Legacy inputs | All 2,860 preserved and excluded from campaign release |
| Source/domain identity | 496 exact KLIFS fingerprints; 770 kinase domains across 435 UniProt proteins; explicit exclusions |
| Preparation comparison | Original CE screen: 112 candidates. Separate explicit-pocket comparison: 161, including all five known-target controls; neither is a complete proteome screen |
| Baseline linker/scoring | Both heads constrained; explicit bounded embedding; invariant rotamer orientation frame prevents scoring-stage head displacement |
| Full-complex QC | Ligand identity/stereo, ligand internal nonbonded clashes, ligand versus all protein chains, E3-complex versus POI overlap, and fixed-head displacement checked |
| Frozen global pilot | All five compounds complete without execution failure: SJF1521 two and SJFδ three QC-valid poses; GSK215/JNJ1013/SJF8240 no valid pose in the small budget |
| Full conformer budget | Native-frame GSK215 control: 100 requested/generated conformers in its feasible local geometry and one QC-valid final pose; not independent global recovery |
| Component resource controls | 100 requested conformers: SJF1521 generated 74 and SJFδ 97; both pass actual repacking/head checks; local process peak roughly 1 GiB and component runtime roughly 44 seconds |
| Ligand MD parameters | All five neutral GAFF2/AM1-BCC sets and GSK215 +1 set verified, including QM-output stereo and atom maps |
| Chemical-state scope | Provisional pH 7.4 plus 6.5–8 sensitivity; no populations or finalized primary state claimed; protein states require assay-aware review |
| MD protocol validation | Neutral GSK215 EM converged (Fmax 472.168<500), followed by complete 200 ps NVT + 500 ps restrained NPT + 500 ps unrestrained NPT |
| MD diagnostics | Receipt-bound energy, ordered-solute topology and periodic contacts verified; 26 sampled final-window frames retain both protein contacts; no production/convergence certificate |
| Execution integrity | Frozen campaign/attempt code, atomic manifests, explicit seeds, process-group termination, source/input drift detection, and verified cached completion |
| Runtime compatibility | Biopython CE interfaces normalized for 1.84 and 1.85; 91 tests pass in both Python 3.10 docking and Python 3.13 analysis environments |
| Storage | Full raw retention accepted; original allowance ~7.24 TB for 112 pairs. A new small native pilot projects ~3.37 TB but does not justify reducing that allowance or certify capacity. Ernie /home has ~1.50 TB free |
| ML transfer | New feature schema explicitly rejects assumed old-model compatibility; experimental holdouts and complete feature validation still required |

A separate bounded local GSK215 global-recovery benchmark is in progress: 50
global poses, five local poses each, two conformers per local pose, two workers,
20,000 maximum global trials and a one-hour time limit. Its exploratory recovery
criteria were recorded before execution. Native-frame success is not substituted
for this independent global search. Active benchmark files are deferred from
installation until they stop changing.

Two study inputs remain unresolved: the original ten-compound aim versus the five
currently verified compounds, and the allocated Ernie data directory. Read-only
software probes did not discover Gromacs in PATH or the limited locations checked;
this is not proof that no installation exists. No scheduler/allocation is invented.

Further publication evidence includes independent recovery/cross-docking,
sampling adequacy, finalized cohort and chemical states, measured Ernie resources,
replicated/converged MD, appropriate experimental validation and frozen analysis.
See NATURE_READINESS.md, ML_TRANSFER_SCOPE.md and methods/04_MD_Protocol.md.
Software validation is not a guarantee of journal acceptance.

Latest software validation: 91 tests passed in each runtime. Raw failure/pass logs,
version observations and installation-specific regressions are retained in audits/.
All historical failures and raw artifacts remain preserved.
