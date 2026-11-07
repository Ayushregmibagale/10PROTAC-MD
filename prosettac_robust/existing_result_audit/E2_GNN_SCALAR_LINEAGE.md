# Existing E2/GNN scalar lineage audit

Read-only reconciliation of current retained files, preserving every original value and prediction. A declared source path and arithmetic agreement do not prove the exact historical training inputs without an immutable execution receipt.

The saved runner `04g_e2_scalar_gnn.py` declares `analysis/features/ml_table.tsv`. Its builder reads `e2_reachability_all_pois.tsv` and `features/lysine_features.csv` for the four E2 scalars. The AXL-only `axl_mechanics/e2_reachability_corrected.tsv` is a distinct output; do not assume its existence means it fed the saved model.

| Field | Declared builder arithmetic / actual interpretation |
|---|---|
| e2_mean_min_dist_A | Equal-pose mean of each pose's frame mean of the minimum modeled E2-ensemble-point/kinase-Lys-NZ separation. This is over the producer's frames, not only bridging frames. |
| e2_delta_reach_50A / 80A | Mean of available within-pose bridge-minus-nonbridge reach differences. Both conditional values must exist; counts differ when a pose has no frames in one state. The historical state mask measures Elongin-B/POI proximity. |
| lk_lys_nearest_e2_dist_A | Rank-sorted first nonmissing static lysine distance, using a single rank-specific modeled E2 point. It is neither the nearest distance across all MD frames nor the minimum over the E2 ensemble. The source comment claiming a mean fallback does not match `.first()` behavior. |

Of 68 current scalar values across the 17 reference targets, 67 reconcile to the builder arithmetic at its retained output precision. MAPK14_full/e2_delta_reach_80A is saved as 0.022575, while the within-pose reconstruction is 0.02569375. Its difference-of-available-means alternative is 0.02256875, within the retained table precision of the saved value. The alternative uses 18 bridge-pose means and 16 nonbridge-pose means; the paired contrast uses 16 poses. This demonstrates a compatible alternative estimand, **not** proof of a manual patch or of the historical execution sequence.

The primary mismatch remains flagged. Do not overwrite it to force agreement. `gnn_scalar_lineage_v1/aggregation_sensitivity.tsv` compares both estimands for all eight MD targets without labels or model performance. A single future aggregation policy should be chosen for its scientific meaning and recorded before evaluation.

The all-target E2 producer correctly anchors to Elongin C/chain B, but truncates to the minimum trajectory/metric frame count without checking timestamps. It unwraps protein fragments and computes direct modeled-point/lysine distances. Those facts do not establish coherent imaging of a whole noncovalent complex; reference-frame and imaging checks remain a distinct task. No new MD is required to inspect them.

These scalar meanings are geometric descriptors. They are not catalytic ubiquitination probabilities, and concatenating them into a GNN does not make its scores calibrated degradation probabilities. Original model outputs remain preserved; a corrected feature evaluation would be separately versioned and assay-reviewed.

Coverage: 33 finite numeric pairs, 35 matching missing-value pairs, 32 numeric matches and 1 numeric mismatch. Matching missingness is not numerical validation of an unobserved descriptor.

## Additional source evidence: target-specific investigation script

The retained `16b_mapk14_e2_fix.py` explicitly patches MAPK14's `e2_delta_reach_80A` using the difference of available bridge/nonbridge means and rounds it to six decimals. This provides a source-level pathway consistent with the mixed estimator found above. It does not establish when that code ran or which inputs it used without execution provenance. Its header says to retain LYS211, while its executable worker still excludes that residue and writes original aggregate tables if run. Treat it as a historical investigation, not a safe routine pipeline step; this audit did not execute it or change any source table.

## Static reference-frame and ensemble checks

All 141 retained random-pose structures have the same Elongin-C Cα coordinates as the common reference (raw RMSD 0 Å), with matching ordered residue identities. The producer's eight selected reference structures are also checked explicitly. The E2 JSON has 13 keys with 51 points each, but all arrays are byte-identical: flattening yields 663 records representing only 51 unique positions, each repeated 13 times. Equal duplication leaves minima and equal-weight fractions unchanged; it must not be counted as independent structural diversity. Evidence: `e2_reference_frame_v1/summary.json`.

The static checks support a common reference frame; they do not validate dynamic imaging, a Gaussian flexibility prior or reaction probabilities. The source perturbs one structural point with 50 Gaussian samples at sigma 20 Å and a fixed seed. That is an assumed positional ensemble, not an equilibrium ensemble sampled by the MD trajectories.
