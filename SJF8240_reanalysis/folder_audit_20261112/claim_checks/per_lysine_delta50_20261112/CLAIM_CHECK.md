# Check of the per-lysine delta_reach_50A claim — 12 November 2026

The quoted summary does not exactly reproduce the saved e2_geometry table. ABL1 has0/18 positive deltas; AXL has2/16 (12.5%), not1/16. AXL positive residues are225 (+0.0114) and272 (+0.0014), in simulation residue numbering. Confirmed degrader positive fractions range7.7–78.9%, with many exact zeros. Positive means strictly greater than zero, not nonnegative.

| Target | Legacy positive / lysines | Legacy fraction | Corrected sampled positive / lysines |
|---|---:|---:|---:|
| ABL1 | 0/18 | 0.0% | 0/18 |
| AXL | 2/16 | 12.5% | 1/16 |
| DDR2 | 7/14 | 50.0% | 1/14 |
| EPHA2 | 9/20 | 45.0% | 0/20 |
| EPHB2 | 5/13 | 38.5% | 2/13 |
| MAPK14_full | 15/19 | 78.9% | 1/19 |
| MET | 2/19 | 10.5% | 0/19 |
| RIPK2_full | 1/13 | 7.7% | 1/13 |

## Significance check

The old lysine-level Mann–Whitney calculation reproduces p=3.030637600329807e-9, displayed as0.000000 at six decimals. It is not literally zero. This treats98 degrader lysines and34 non-degrader lysines as independent observations. They are clustered within only6 degraded and2 non-degraded proteins; this is not a valid independent-protein significance calculation.

Using one summary per protein, exact two-sided rank tests give p=0.07143 for mean delta and p=0.28571 for fraction-positive lysines. The legacy mean delta ranks all6 degrader targets above both non-degraders, an interesting descriptive pattern worth preserving. It is retrospective and not independent validation. An exact label permutation of the mean-difference statistic instead gives p=0.03571 (absolute-difference convention) for mean delta and0.21429 for fraction positive; these are different test statistics/two-sided conventions, not interchangeable p-values. All28 assignments and both investigated summaries are retained in verification.json, rather than choosing the most favorable test. Protein labels are observational and may be confounded by protein identity/family; permutation inference assumes exchangeability.

## Source and corrected-geometry check

The retired e2_geometry/03 script reads the original protac-md/analysis/per_lysine_e2_summary.tsv. Its upstream20_per_lysine_e2_distance.py uses the legacy metrics.csv bridge mask (the historically named VHL contact actually selects ELOB), crops frames with min(trajectory length,mask length), and does not apply the later ligand-connected coherent-image placement. It is not the corrected true-VHL feature table.

A supplemental calculation from all141 hash-verified mechanistic_v3 pose packets uses true-VHL heavy-contact bridging and corrected coherent E2/lysine images. It retains41 sampled frames per pose and averages within-pose bridge-minus-nonbridge contrasts only for poses containing both states. Its positive fractions are shown above: EPHA2 and MET both have0 positive lysines, as does ABL1. Thus the broad separation in the quoted legacy summary does not reproduce in these corrected sampled data. Sampling, bridge definition and imaging all differ, so this is not an isolated causal comparison; a matched full201-frame recomputation would be needed to isolate those changes. This supplemental check is not a replacement for that full-resolution validation.

Negative delta means reduced modeled50-Angstrom reach occupancy during the chosen bridge state. It does not mean a lysine is physically buried, chemically inaccessible or unable to undergo ubiquitination. No new GNN was fitted, no label changed, and no original result was overwritten.

## Reproduction

Counts, full corrected sampled per-lysine values, per-target paired-state coverage, exact statistics and source hashes accompany this report. Calculation source is preserved in the persistent workspace SJF8240_per_lysine_claim_check_20261112/check.py with its SHA256 in verification.json. The canonical audit source registry is unchanged because this delivery adds report/data artifacts only.
