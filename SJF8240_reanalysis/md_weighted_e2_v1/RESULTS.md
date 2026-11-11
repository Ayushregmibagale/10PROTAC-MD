# Class-weighted MD GNN with E2 reachability features

Eight-target corrected-label LOPO: six degraded binders (DDR2, EPHA2, EPHB2, MAPK14_full, MET, RIPK2_full); two non-degraded binders (ABL1, AXL). All targets have MD. The original two-layer MD GAT backbone receives two graph-level E2 geometry descriptors at its classifier head. No static encoder or conventional predictive model is used.

| Arm | AUROC | Accuracy | Balanced accuracy | Sensitivity | Specificity |
|---|---:|---:|---:|---:|---:|
| ELOB control | 0.0000 | 12.5% | 8.3% | 16.7% | 0.0% |
| True VHL | 0.0000 | 25.0% | 16.7% | 33.3% | 0.0% |

Loss uses balanced per-target weights N_train/(2*N_class), derived from the seven training targets only. E2 means and standard deviations are also fitted only to those seven targets, saved in each fold checkpoint, and applied to the held-out target without refitting. No missing target scalar is zero-filled. The epoch budget (300), seed (42), optimizer, dropout and threshold remain fixed before fitting. No outcome-based tuning was performed. An always-degraded call gives 75% accuracy in this six-versus-two cohort; only two negatives limit the interpretation of specificity and generalization.

The E2 inputs are: (1) equal-pose mean of the time-averaged minimum E2-point/lysine-NZ distance in Angstrom, across all sampled frames; (2) the equal mean of within-pose bridge-minus-nonbridge 80-Angstrom geometric reachability fractions, using only poses containing both states. They use the audited ligand-connected image convention and true-VHL heavy-contact state mask. All 141 core cases/28,341 frames and the target feature table were independently reproduced from existing packet/contact files. Paired-state coverage varies from 4 to 10 poses per target and is retained in provenance. Both graph roles receive the same E2 values, so role comparison changes the MD graph selection alone.

These quantities are modeled geometric proxies, not experimental E2 occupancy, steric accessibility, catalytic probability, or evidence that all lysines are buried. Corrected ABL1 delta80 is approximately -0.03685, not the legacy -0.0466 value. Frames with alternative image sensitivity remain documented; no pose is removed to improve a model score. E2-candidate data are used at the user's explicit request for exploratory fitting, not upgraded to scientific certification.

All 16 saved model folds passed independent CPU replay. Maximum differences: ELOB 2.98e-08; VHL 5.96e-08. Source/input hashes, class weights, exact held-out membership, training-only E2 scalers and metrics were independently checked. No new docking or MD was run. The unweighted and weighted MD-only versions remain unchanged. Adding E2 expands the head by two inputs, so any difference includes both descriptor information and the small added head capacity.

No E2-augmented SIK3/RPS6KB2 predictions are produced: audited E2 descriptors for those extension targets are absent, and no fabricated zeros are used. This retrospective, eight-protein LOPO evaluation does not establish independent biological validation or calibrated degradation probabilities. Original structures, trajectories, assay tables and earlier results are preserved.
