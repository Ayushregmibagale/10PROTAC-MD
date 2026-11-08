# Retained-structure E2 imaging review — 2026-11-08

This is a retrospective geometric sensitivity analysis, not replacement MD, new GNN predictions or a catalytic probability estimate. Original outputs remain unchanged.

## What changed, how and why

Separate periodic molecule images can distort direct distances to modeled E2 points even after each molecule is made whole. `workflow/scripts/lib/coherent_images.py` validates triclinic cells and reconstructs a coherent assembly using exact integer lattice translations. Elongin C anchors the reference fit; true VHL uses the reference-predicted image. A complete covalently reconstructed ligand connects VHL to the target through closest heavy-atom contacts. Independent reference and protein-contact image conventions provide sensitivity controls. E2 points do not choose the images. Intramolecular conformations remain unchanged. This is an explicit geometric convention; it does not demonstrate a unique physical assembly after dissociation.

## Completed evidence

The full core cohort contains 141 poses from eight targets, each with 201 unique frames over the explicit 0–20 ns window: 28,341 frames. Hash-bound existing trajectory/metric receipts and duplicate-time rules preserve the separate extended AXL series. No cases failed. All 22,753 true-VHL/ligand/target heavy-atom bridging frames agree across the three image conventions. Ligand versus protein-contact images differ in 16 nonbridging frames: AXL pose007 (one), EPHB2 pose003 (15). Reference versus protein-contact images differ in ten frames. These findings apply to the sampled frames and conventions, not all possible assemblies.

Independent integer-box translations were tested on 715 real frames: maximum descriptor discrepancy 9.53e-13 Å. Maximum reproduced periodic heavy-contact discrepancy is 0.000020 Å. Ten poses change mean minimum E2 distance by more than 0.001 Å relative to the old representation. E3 reference-fit RMSD can reach approximately 16 Å; no rigid-E3 or conformational-convergence claim follows.

The 51 modeled E2 points are a structural/geometric prior including Gaussian perturbations, not 51 independent experimental observations or an equilibrium ensemble. The same point cloud occurs in multiple retained files. No full E2 protein steric or reaction model is certified here.

## Publication artifacts and interpretation

`existing_result_audit/coherent_e2_publication_v1` contains PDF/SVG/PNG figures, per-pose and per-target tables, paired conditional-state comparisons and minimization-exception sensitivity. The figure compares old and ligand-connected image conventions for retained poses and equal-pose target means. Lines connect the same target; they are not confidence intervals. Poses/frames are not independent biological replicas.

Bridging means simultaneous strict heavy-atom distances <5 Å from the ligand to true VHL and target. This is geometric occupancy, not thermodynamic stability or degradation probability. Historical ELOB-based and corrected true-VHL masks are retained separately; conditional contrasts use within-pose paired states, with absent states missing rather than zero. Main summaries retain every pose; omission of the EM exception is sensitivity only.

The historical main GNN consumes four E2 scalars: mean minimum distance, paired reach contrasts at 50 and 80 Å, and a static nearest lysine/E2 distance. The new image analysis does not automatically replace those inputs or enter old weights. The two extension predictions belong to the separate predictor without E2 scalars.

Reproduction tools and immutable source hashes are under `existing_result_audit/review_tools/` and `ligand_connected_e2_cohort_v1/`. The earlier pilot and failed unit-test attempt remain retained. Three added tests cover lattice invariance, invalid cells and whole-ligand paths. No biological model was fitted.
