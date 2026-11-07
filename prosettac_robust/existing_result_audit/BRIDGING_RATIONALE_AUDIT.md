# Existing-result bridging rationale audit

**Confirmed protein-role error in all 141 core systems:** every legacy group named VHL selects Elongin B/Q15370 (Protein_chain_A). True VHL/P40337 is Protein_chain_C, with 100% sequence identity and observed-sequence coverage in all cases. The POI groups match their expected targets. Thus the table below reproduces **legacy-named group proximity**, not VHL retention. See `protein_roles_core_complete/summary.json`. Corrected contact re-analysis is complete for all 141 existing trajectories; no new MD was executed.

Scope: read-only audit of the original SJF8240-labelled experiment, not a new docking/MD campaign. Preserve original structures, results, labels and predictions. Existing chemical-identity concerns remain separate provenance findings; this report does not reinterpret old trajectories as the newly verified supplier graph.

## What 5 Å does and does not establish

Five angstroms can be a useful operational contact cutoff. CAPRI/DockQ uses a 5 Å heavy-atom definition for protein-interface contacts; that is context for a geometric threshold, not validation of a PROTAC degradation rule or of the different atom selections in this project. [Basu & Wallner, 2016](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0161879).

For the existing MD metric, define B_c(t) = 1[d_ligand,VHL(t) < c AND d_ligand,POI(t) < c]. The reported occupancy is the average over sampled, aligned frames. The contact to each protein can involve different ligand atoms. This does not establish retention of both cognate head–pocket interactions, a productive ubiquitination geometry, equilibrium stability or calibrated degradation probability. Experimental work shows that ternary-complex formation alone is insufficient for ubiquitination/degradation. [Smith et al., 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6328587/).

## Actual historical measurements differ

| Component | Actual atom definition / calculation | Interpretation and action |
|---|---|---|
| PRosettaC pose screen (`01_screen_bridging_all_pois.py`) | Minimum distances from chain-X ligand atoms excluding names starting H to configured protein chain/residue subsets; protein hydrogens are not explicitly excluded. Strict `<5` for both. | Geometric pose-screen fraction among sampled outputs. Document actual selections; this is not a thermodynamic P(productive). Atom-name filtering is not an elemental definition. |
| MD (`02_run_analysis.py`, called by `08_run_random_screen.py`) | `gmx mindist` on stored Other/VHL/Kinase groups; native nm converted to Å; strict same-frame `<5` conjunction. | Existing Other groups have 144 atoms including 65 hydrogen-named atoms; ligand-heavy count is 79. Therefore do not call this a heavy-only contact metric. Named protein groups also contain hydrogens. |
| Dynamic GNN (`04c_build_md_contact_graphs.py`) | Ligand selection `resname UNL and not name H*` versus protein Cα atoms, direct Euclidean distances. Edges multiply marginal contact frequencies. | A different Cα-proximity proxy. Product P(V)P(K) need not equal observed same-frame P(V and K). Old direct subtraction does not apply minimum-image PBC. Preserve its original features/predictions and describe them honestly. |
| Repaired derivative producer | Explicit executed role mapping, periodic distances and measured joint contacts; explicit labels/provenance, 3 node features. | New feature definition requires its own version/evaluation. It cannot retroactively validate old model performance. Existing trajectories can support a separate re-extraction; MD rerunning is not intrinsically required. |

GROMACS 2025.4 `mindist` takes PBC into account by default, matching the saved XVG command/version. This code audit does not independently recompute every stored periodic distance. [Official manual](https://manual.gromacs.org/2025.4/onlinehelp/gmx-mindist.html). Direct coordinate subtraction is not generally a minimum-image displacement. [MDAnalysis documentation](https://docs.mdanalysis.org/stable/documentation_pages/lib/distances.html).

## Existing-result reconciliation and cutoff sensitivity

The initial all-manifest audit verified 140 cases and flagged one AXL time-series issue plus 33 prediction-target cases without the expected local XVG tables. Missing tables are not missing simulations or biological negatives. Keep the initial report as an explicit incomplete audit.

For the eight core MD target groups, the explicit second audit covers **141/141 cases**. Every original saved CSV and rounded bridging header reconciles to its original XVG data before any duplicate/window handling. In AXL/pose011, the available trajectory analysis extends to 100 ns and repeats the identical sample at 20 ns. Preserve that original record. The derivative audit removes only the exact repeated sample and uses a common 0–20 ns window, yielding 201 samples per case at 0.1 ns spacing. Conflicting duplicates or mismatched timestamps are rejected. No original output is changed.

Keep 5 Å as the inherited primary rule. The following checks use no biological labels and do not optimize the cutoff. Values are equal-pose mean simultaneous-proximity occupancies (%), not degradation rates:

| Target | Poses | 4 Å | 5 Å | 6 Å |
|---|---:|---:|---:|---:|
| ABL1 | 17 | 29.27 | 34.27 | 43.23 |
| AXL | 18 | 25.37 | 35.16 | 46.96 |
| DDR2 | 18 | 21.95 | 31.12 | 40.52 |
| EPHA2 | 18 | 20.67 | 27.22 | 38.11 |
| EPHB2 | 19 | 9.87 | 18.54 | 29.35 |
| MAPK14_full | 18 | 39.50 | 46.63 | 53.92 |
| MET | 16 | 24.35 | 28.67 | 35.98 |
| RIPK2_full | 17 | 34.18 | 44.69 | 55.58 |

The absolute occupancy is cutoff-sensitive. Report this sensitivity and the exact atom definition instead of presenting 5 Å as a uniquely validated physical boundary. Frames are correlated and starting structures were bridge-selected; pose counts and frame counts are not independent experimental sample sizes. All 141 poses are retained in the stated core cohort, including poor-occupancy cases.

## Claims to keep descriptive

- Replace unqualified P_accessible/P_stable/P_ubiquitinable/P_degradation language with explicitly defined sampling/proximity/reachability scores unless independently calibrated. Multiplying heuristic fractions does not by itself produce a degradation probability.
- Distinguish any-ligand proximity from correct head/pocket engagement; pocket-specific checks can reuse existing structures/trajectories in separate audit outputs.
- Record duration, sampling, target/pose selection and actual per-frame atom groups. The AXL extension shows why a hardcoded “all 20 ns” caption is inadequate.
- Preserve original GNN predictions. Re-extracted periodic/joint features constitute a separate evaluation, never a silent replacement of the reported dataset.

## Reproduction

`bridge_cutoff_v2_core_20ns/summary.json` records scope, cutoffs, duplicate/window policies, source and output hashes. `pose_provenance.json` binds each result to existing XVGs, CSV, index and GRO. `sources/` and `source_inventory.json` retain the audited source/helper bytes. No docking executable, GROMACS simulation or original script is executed by the audit.


## Direct existing-trajectory validation

A deterministic first case, DDR2/pose001, was checked on its original 201-frame, 0–20 ns trajectory. Periodic distances for the stored groups reproduce GROMACS outputs with maximum errors 0.0000137 Å (legacy-labelled VHL/actually Elongin B) and 0.00000865 Å (POI). This validates the calculation reproduction, not the legacy role label. With sequence-identified VHL, the same trajectory has 100% simultaneous 5 Å proximity occupancy versus 1.99% for the legacy-defined groups. The heavy-only true-VHL/POI definition also yields 100% at 5 Å for this case. These numbers cannot be generalized to the full cohort before the batch finishes, and none is a degradation probability.

Consequently the legacy dynamic GNN branch's chain-A/Cα features require description as Elongin-B/E3-scaffold proximity, not validated VHL-pocket stability. The original predictions may be retained as historical model outputs; their semantics and limitations must be corrected. New true-VHL/periodic/joint features, if evaluated, are a separate analysis on existing MD and cannot silently inherit old performance claims.


## Static-graph representative check and E2 scope

The read-only static audit parsed the old graph configuration as AST data (it did not execute the original script) and selected one available retained combined PDB per core target for sequence identification. All 17 representatives identify the legacy chain-A/residue-1–103 Cα selection as Elongin B/Q15370, not VHL. This verifies the sampled representatives and shared selector; it is not an exhaustive re-analysis of every pose. Evidence: `static_role_representatives_v1/`.

The existing corrected E2 script explicitly uses Elongin C/chain B (85 Cα) for frame alignment, consistent with the sequence audit. Therefore this role finding does not blanket-invalidate that transformation. Its bridging-conditioned subsets and associated interpretation still need to distinguish the legacy Elongin-B proximity criterion from true VHL contacts.


## Completed core-cohort contact re-analysis

All 141 cases completed on the original saved trajectories, with 201 common-window samples per case (28,341 sampled frame records). Every stored-group GROMACS distance check passes; maximum errors are 0.00002093 Å (legacy-labelled VHL/Elongin B) and 0.00002114 Å (POI). Source snapshots, per-case hashes, time/duplicate rules, outputs and a 2-worker plan are retained. The analysis took about 499 s; it ran no molecular simulations.

Equal-pose mean simultaneous proximity at the inherited 5 Å cutoff (%):

| Target | Poses | Legacy Elongin B + POI, all atoms | True VHL + POI, all atoms | True VHL + POI, heavy atoms |
|---|---:|---:|---:|---:|
| ABL1 | 17 | 34.27 | 85.25 | 79.81 |
| AXL | 18 | 35.16 | 89.03 | 86.15 |
| DDR2 | 18 | 31.12 | 78.03 | 74.07 |
| EPHA2 | 18 | 27.22 | 73.52 | 70.23 |
| EPHB2 | 19 | 18.54 | 82.01 | 77.90 |
| MAPK14_full | 18 | 46.63 | 96.96 | 94.50 |
| MET | 16 | 28.67 | 90.08 | 89.43 |
| RIPK2_full | 17 | 44.69 | 72.46 | 70.76 |

These corrected measurements materially change the meaning and magnitude of the old “bridging” descriptor. They describe sampled contacts in the existing simulated systems, not degradation rates or proof of native pocket engagement. High VHL/POI proximity in these initially selected complexes is compatible with further determinants of degradation; it does not validate a sole-contact classifier. Keep the legacy E3-scaffold descriptor available under its correct name rather than silently replacing its historical model inputs.

Raw per-case contact tables and definition/cutoff summaries are in `contact_core_v1/`; the review figure is exported as vector PDF/SVG and PNG under `contact_core_v1/figures/`. No uncertainty interval is claimed for the plotted descriptive means.
