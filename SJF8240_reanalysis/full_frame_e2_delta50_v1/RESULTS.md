# Full-frame true-VHL per-lysine reach50 comparison

All141 retained cases and201 frames per case (28,341 total) were recomputed. Same50 A threshold and eight current-labelled proteins; ABL1/AXL are negatives, EPHA2 is positive. SIK3/RPS6KB2 remain unknown and excluded. No GNN fitting or new MD.

| Target | Legacy positive lysines | Corrected true-VHL positive lysines | Legacy mean delta | Corrected mean delta | Corrected paired poses |
|---|---:|---:|---:|---:|---:|
| ABL1 | 0/18 | 0/18 | -0.01965 | -0.01169 | 9 |
| AXL | 2/16 | 1/16 | -0.01676 | -0.01034 | 8 |
| DDR2 | 8/14 | 0/14 | +0.00101 | -0.00772 | 7 |
| EPHA2 | 9/20 | 0/20 | -0.00686 | -0.00715 | 10 |
| EPHB2 | 5/13 | 2/13 | +0.00178 | -0.00543 | 6 |
| MAPK14_full | 15/19 | 3/19 | +0.01615 | -0.02227 | 4 |
| MET | 2/19 | 0/19 | +0.00219 | -0.00066 | 4 |
| RIPK2_full | 1/13 | 1/13 | -0.00409 | -0.00721 | 7 |

## Protein-level evidence

- mean_delta: exact two-sided rank p=0.285714; Holm-adjusted p=0.571429; mean degrader-minus-nondegrader difference=+0.002607.
- fraction_positive: exact two-sided rank p=0.714286; Holm-adjusted p=0.714286; mean degrader-minus-nondegrader difference=+0.033527.

The independent units are8 proteins, not132 lysines or28,341 frames. Exact rank p-values enumerate all28 allocations of two negatives, with averaged ranks for ties and twice the smaller inclusive tail. Difference-statistic permutation p-values are separately labelled in the tables. These are retrospective, observational comparisons, not external validation; labels may be confounded by protein identity. A small p-value alone would not establish a degradation predictor. Positive means strictly delta>0; a1e-12 numerical-zero sensitivity count is retained.

## Meaning and controls

For each lysine/frame, take the nearest modeled E2 point distance. Reach50 is the fraction of bridge/nonbridge FRAMES with that distance<50 A. Delta is bridge minus nonbridge, calculated within each pose and then equally averaged over poses containing both states. It is not the fraction of E2 points and not evidence of solvent burial or ubiquitin-transfer chemistry. The reference file's663 concatenated rows contain51 unique coordinates; duplicate rows do not affect a nearest-point distance.

All four geometry/mask combinations use identical201-frame grids. The corrected geometry uses ElonginC fitting and ligand-connected whole-component imaging. True-VHL bridging uses periodic ligand/VHL and ligand/target heavy-atom contacts<5 A. The legacy mask uses the original metrics.csv convention, whose VHL label actually selected ELOB. The geometry factor bundles reference/image convention; the mask factor bundles identity/atom-selection policy. We do not attribute a change to one subcomponent alone.

The primary comparisons use every eligible paired-state pose for each mask. Matched-support sensitivity uses only poses eligible under both masks, preventing different pose subsets from driving that particular comparison. Coverage tables retain all poses, including those with no valid contrast. No unsupported contrast is imputed as zero. All8 targets must be present before protein-level tests are computed. Distances at the old41 sampled times are checked against the saved mechanistic_v3 packets.

The legacy per-pose and published per-lysine tables are independently checked to their rounding precision. Original files, labels, old GNN scores and older analyses remain unchanged. This result tests the specific reach50 hypothesis; it does not retrain or retrospectively validate the0.750 temporal GNN, which did not explicitly receive per-lysine delta50.

## Files

results/target_summary.tsv and legacy_vs_corrected.tsv show all targets and conventions. results/per_lysine.tsv and per_pose_per_lysine.tsv retain each lysine. results/protein_statistics.tsv and primary_statistics.tsv contain complete exact tests; matched_support_*.tsv contains pose-matched sensitivity. results/state_coverage.tsv records state support. validation/ records frozen-source, legacy/coarse reproduction and independent unwrap replay. Local analysis/poses/ contains full distance packets with source hashes.
