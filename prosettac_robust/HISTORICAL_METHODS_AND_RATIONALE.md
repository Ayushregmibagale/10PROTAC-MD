# Retained-structure workflow: methods and rationale

## Current scope

Use the structures, parameters and completed results from the original experiment. Review and harden PRosettaC, MD, contact/E2 analysis and GNN code. Do not replace the molecular structure, rerun the study or change the labels/predictions. The historical simulated connectivity is the reference for this review. Supplier equivalence is a separate provenance question, not a prerequisite for continuing this script/rationale work. State the actual simulated identity when writing a manuscript.

The confirmed panel is 17 historical GNN evaluation targets plus SIK3 and RPS6KB2. A shared script does not mean all targets have equivalent inputs, MD coverage or experimental evidence. Coverage is recorded in `existing_result_audit/scientific_signoff_v1/coverage_19_targets.json`. Directory names such as `degraders` describe historical organization; they are not an independent source of experimental truth.

## What the pipeline actually measures

| Stage | Implemented measurement | Defensible rationale and interpretation |
|---|---|---|
| PRosettaC sampling and scoring | Finite sampled protein arrangements, linker completions and Rosetta scores | Generate geometrically feasible hypotheses. Record attempted/rejected/valid counts, head constraints and selection. Search yield depends on budget and protocol. It is neither a biological negative when zero nor a thermodynamic probability when nonzero. Rosetta score order does not by itself establish which pose is most stable in solution. |
| Historical static proximity screen | Minimum ligand-heavy/name-filtered distances to configured chain/residue groups; strict <5 A conjunction | Operational screen of retained structures. The original code's `P(productive given ternary complex)` wording overstates the measurement. The sampled fraction is conditional on that search and selection, not a calibrated productive-complex probability. In the checked static representatives, the group called VHL is actually Elongin B. |
| Historical MD distances | Stored GROMACS minimum distances, with hydrogens included in the old groups | Track proximity over saved frames. Name the exact groups, PBC, units and selected time interval. Nm from GROMACS become A by multiplying by ten. Same-frame conjunction is essential. Any two ligand atoms can supply the contacts; this is not proof of retained head/pocket binding. |
| Corrected retrospective contact analysis | Sequence-identified true VHL and POI; explicit all-atom and heavy-atom alternatives; triclinic PBC | Reuse the same simulations while making selections and definitions correct and visible. Keep old metrics/predictions and corrected measurements separately versioned. Correction changes feature meaning; it does not retroactively change the evaluated GNN. |
| E2 descriptors | Distances/reach fractions to modeled E2 positions and target lysines | Geometric reach hypotheses under an assumed E2-position prior. The Gaussian sample is not an MD equilibrium distribution. Contact-conditioned contrasts must identify their state mask and paired/unpaired estimand. Do not call these catalytic probabilities. |
| Historical GNN | Static and dynamic graph encoders fused with four E2 scalar inputs | Learn an association with recorded labels under target-held-out evaluation. Preserve architecture and original predictions. A sigmoid output is a model score, not an experimentally calibrated degradation probability. Label provenance, cohort selection, missingness and feature semantics remain part of its interpretation. |

The original dynamic graph definition uses ligand-heavy to protein-CA proximity, direct Euclidean distances and products of marginal frequencies. This differs from all-heavy contacts and from observed same-frame joint frequency. The new supplier-specific future producer has corrected periodic/joint semantics but is not used to overwrite or reinterpret the historical dataset. Read-only retrospective tools use the old structures directly.

The retained local adapter replaces the PatchDock executable with random rigid-body sampling constrained by anchor-distance and CA filters, followed by PyRosetta local docking/linker work. `Patchdock_Results` is a compatibility directory name, not proof that PatchDock ran. Describe this as an adapted PRosettaC/PyRosetta workflow, with the actual algorithm and source version. The future centered-linker controls did not generate the historical results.

## Actual compiled MD protocol, not today's template

All 174 production TPRs associated with the saved MD graph pose records were inspected with compatible GROMACS 2025.4 using `dump -om`. No `grompp` or `mdrun` was executed. The compiled input records show:

- MD integrator; 0.002 ps step; 10,000,000 nominal steps (20 ns).
- Compressed coordinate interval 50,000 steps (100 ps).
- PME electrostatics, Verlet cutoff scheme and potential-shift Coulomb/Lennard-Jones modifiers; dispersion correction off.
- V-rescale at 300 K with tau-t 0.1 ps for the three listed coupling groups.
- Isotropic Parrinello-Rahman pressure coupling at 1 bar with tau-p 2 ps.
- LINCS algorithm, order 4, one iteration; `continuation=false` in the retained input record. The template's `h-bonds` selection is not reconstructed from inputrec alone; constraint interactions are topology information.

**The cutoffs were not uniform:**

| Target | Cases with 1.0 nm Coulomb/LJ cutoffs | Cases with 1.6 nm Coulomb/LJ cutoffs |
|---|---:|---:|
| ABL1 | 17 | 0 |
| AXL | 8 | 10 |
| DDR2 | 18 | 0 |
| EPHA2 | 0 | 18 |
| EPHB2 | 19 | 0 |
| MAPK14_full | 18 | 0 |
| MET | 16 | 0 |
| RIPK2_full | 17 | 0 |
| SIK3 | 0 | 15 |
| RPS6KB2 | 0 | 18 |
| Total | 113 | 61 |

Today's `prod_20ns.mdp` specifies 1.6 nm. It cannot be cited as proof that all historical runs used that setting. Retain protocol class per case, disclose the difference and avoid attributing a cross-target difference solely to biology. This is a methods/provenance finding, not automatic grounds to discard trajectories or launch a rerun. Numerical and structural sensitivity cannot be inferred merely from the parameter difference.

TPR duration is a requested protocol length, not proof of actual completion or the origin/duration of a postprocessed trajectory. The retained AXL/pose011 analysis trajectory extends beyond 20 ns; the common-window audit explicitly uses 0–20 ns and handles its identical duplicate at 20 ns. The original preprocessing code skips existing `prod_nojump.xtc` by filename alone, so a current file is not an immutable as-run transformation receipt. Keep raw trajectory, processed trajectory, time window and execution provenance distinct.

The saved historical dynamic GNN metadata has 21 sampled frames per pose across all 174 cases, consistent with subsampling a 201-frame, 100 ps output series by stride 10. Do not repeat the source comment claiming approximately 1,000 samples. Reconstructed corrected contact summaries use 201 frames per core pose. Those sample counts describe different features.

## Why retain 5 A

Keep 5 A as the inherited primary operational cutoff; do not tune it to increase GNN performance. Explicitly state strict `<`, exact atoms, periodic boundaries, simultaneous contacts, frame interval and equal-pose aggregation. Report 4/5/6 A sensitivity as descriptive analysis. The already saved sensitivity tables use the original structures and do not change experimental labels. A heavy-atom interface convention can motivate a scale but does not validate a different CA-proximity definition or a degradation threshold. Primary references and detailed definitions are in `existing_result_audit/BRIDGING_RATIONALE_AUDIT.md`.

## Code repairs under this retained-structure scope

1. The shared executed heavy-atom mapping now rejects fractional/boolean/empty indices, overlaps, solvent indices outside the mapped solute, hydrogens and nonfinite masses. Contacts and GNN production loading use the same validation.
2. Contact analysis rejects nonfinite time controls, invalid boxes/coordinates, nonfinite or non-increasing trajectory times, and input drift. XTC reading uses the qualified sidecar-free reader, preserving source directories.
3. Thermodynamic extraction uses the existing indexed-legend parser, including named temperature/pressure/potential/density ordering and finite/physical-value checks. Requested order is not assumed to equal GROMACS output order.
4. `workflow/scripts/02_analysis/09_reanalyze_existing_contacts.py` is an explicit retained-data entrypoint. It uses the historical trajectory directly, checks both stored distance time axes before combining them, identifies proteins by observed sequence and verifies distance reproduction. It neither requires supplier-replacement parameters nor fabricates a historical production receipt.

The new retained-data CLI exactly reproduces the earlier 201-row DDR2/pose001 contact table byte for byte; its metadata additionally binds configuration and helper sources. It is qualified on the existing core layout with legacy indices/XVGs. It deliberately fails if those inputs are absent; the two prediction-target folders need separate explicit adapters rather than guessed index files.

## Safe reproduction of the retained-data contact check

Use the analysis environment with MDAnalysis/Biopython/NumPy. Select a fresh output directory outside both original projects:

```bash
python workflow/scripts/02_analysis/09_reanalyze_existing_contacts.py \
  --project /home/ayushregmibagale/protac-md \
  --poi DDR2 --pose-id pose001 --time-end-ns 20 \
  --output /tmp/my-fresh-retained-DDR2-contact-audit
```

This runs analysis, not MD. It does not launch docking, change model weights or overwrite old predictions. Compiled inspection files under `compiled_md_protocol_v1` are diagnostic exports and must not be reused as simulation MDP inputs.

## Evidence and limits

`existing_result_audit/historical_methods_v1/source_snapshot` freezes the reviewed source bytes, with original paths and hashes. These are the retained current sources, not proof of every historical execution version. `compiled_md_protocol_v1` binds all 174 settings inspections to actual TPR hashes. The initial 2024.3 GROMACS probe failed because the old TPR format required a compatible reader; that failure and the successful 2025.4 route are recorded. No original file changed.

Code/rationale correctness can be reviewed on this retained-structure basis without substituting chemistry. Biological identity, assay correspondence and predictive validity still limit publication claims. Keep the prior scientific review as evidence of those limits, not as an instruction to replace the experiment.

Current validation: 129 analysis tests pass; 126 docking tests pass with three PyG-dependent GNN replay tests skipped in that environment. The mapped-contact analyzer additionally passed a real analysis-only check on the previously completed separate 20 ps control, explicitly labelled diagnostic. No simulation or historical GNN training was repeated. Full receipts and logs are in `validation_controls/software_checkpoint_v7.json` and `existing_result_audit/historical_methods_v1`.
