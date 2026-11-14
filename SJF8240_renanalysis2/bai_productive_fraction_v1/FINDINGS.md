# Bai-style productive-fraction test — 14 November 2026

The paper's central productive/unproductive rule has now been implemented as a separate VHL adaptation. It does not improve the current eight-target ranking: MD productive-frame fraction AUROC=0.708,versus0.750 for the earlier modified CRL2 score. The closer clustered static-starting-pose fraction gives AUROC=0.583. Neither result is a validated biological classifier or a literal replication of the CRBN paper.

## Rule and denominator

A configuration is productive if at least one target lysine is exposed and≤16Å from ubiquitin's C-terminus,and the target does not clash with the chosen ligase scaffold. This implementation defines the atoms as lysineNZ and modeled ubiquitinGly76 carbonylC,exposure as NZ accessible-sphere fraction≥0.1,and a clash as maximum heavy-atom van-der-Waals overlap>0.8Å. These quantitative atom/exposure/clash definitions are explicit local choices; the available paper methods do not specify them sufficiently to claim exact reproduction.

The extra CE→NZ orientation requirement and the per-frame PROTAC bridging requirement have both been removed. A pose/frame qualifies if either declared VHL scaffold permits it. The main terminal model is fit75;fit74 is a separately reported sensitivity,not a second chance used to increase the score. Genuine clashes remain unproductive and remain in the denominator. No labels,cutoffs or scaffold weights were selected to improve the result.

Productive fraction=number of productive configurations/number of evaluated configurations. A count of lysines,contact atoms,raw docking poses or a docking energy is not this score.

## Results

| Target | User label | MD productive fraction,% | Clustered-starting-pose member-weighted fraction,% |
|---|---|---:|---:|
| ABL1 | Non-degrader | 0.000 | 0.000 |
| AXL | Non-degrader | 2.377 | 0.000 |
| DDR2 | Degrader | 3.096 | 0.000 |
| EPHA2 | Degrader | 6.357 | 5.556 |
| EPHB2 | Degrader | 11.783 | 0.000 |
| MAPK14_full | Degrader | 0.000 | 0.000 |
| MET | Degrader | 1.057 | 0.000 |
| RIPK2_full | Degrader | 1.902 | 0.000 |

At the unvalidated descriptive rule score>0,the MD score has sensitivity0.833,specificity0.500,balanced accuracy0.667. Its exact two-sided target-label permutation p=0.464286. There are eight independent labelled targets,not28,341 independent observations. AXL remains nonzero and MAPK14_full remains zero. This is a retrospective geometric ranking,not LOPO-trained GNN performance.

An explicit3Å anchor-fit QC sensitivity produces a different MD AUROC;it is saved in the full table and is not promoted as the new primary result. Closed-only versus closed-plus-modeled-open scaffolds,both terminal placements,exposure cutoffs0.05/0.1/0.2 and16/20Å distances are all retained in the tables.

## Closer correspondence to the paper

For the pose-count analysis,one production-frame0 configuration per retained MD pose was aligned on trueVHL and clustered at3Å target Cα RMSD. Complete-linkage clustering was used,requiring matched sequence/coordinate correspondence. Each cluster's geometric medoid was evaluated;its productive/unproductive classification was propagated to all retained members. Both equal-cluster and member-weighted fractions are reported. This avoids passing off an MD frame percentage as a clustered-pose percentage.

These are the141 retained starting configurations of completed MD,not the paper's complete Rosetta ternary-model population. The available set does not include a comparable complete clustered Rosetta ensemble for all eight proteins. The paper uses MMTSB clustering and a top interface-energy representative;our retained-pose diagnostic uses complete-linkage and a coordinate medoid. That difference is recorded rather than hidden. No new docking was run.

## What still differs from a literal Bai reproduction

- The published system is CRBN/CRL4A;this adaptation is VHL/CRL2 because SJF8240 recruits VHL.
- The local closed template is a trapped mutant with incomplete sidechains. The open catalytic-module placement is modeled.
- DonorGly76 is absent from the deposited catalytic structure;its endpoint remains modeled. The unpassed native-reference calibration from the earlier analysis is not repaired by simply changing the scoring rule.
- Source chemistry and actual-coordinate stereochemistry of historical trajectories remain separate from supplier-SJF8240 equivalence.
- The full original Rosetta model populations,author-selected representative energies and fully specified exposure/clash/atom criteria are not reproduced here.

Therefore this is the paper-style productive-fraction hypothesis tested on our available data,not an exact reproduction or evidence that the published method fails. A faithful VHL analogue needs validated complete catalytic reference geometry and a suitable full ternary ensemble. Experimental ubiquitination data remain the appropriate validation target.

## Saved files

- config/protocol.json:criteria and explicit reproduction gaps.
- results/primary_target_scores.tsv:MD productive fractions for the eight targets.
- results/frame0_pose_flags.tsv and cluster_target_scores.tsv:pose-count comparison.
- results/cluster_membership.tsv:every cluster,representative and member weight.
- results/metrics.tsv and target_scores_all_variants.tsv:all predefined sensitivities.
- validation/independent_verification.json:all141 pose/28,341 frame Boolean checks,pairwise AUROC verification and member-count arithmetic.
- validation/frame0_coordinate_sources.json:source hashes for all frame0 geometries.
- scripts/:reproducible scoring,clustering and verification/report generation.

Original inputs and prior analyses are unchanged. No new MD or GNN training. SIK3/RPS6KB2 remain unknown and are excluded from labelled evaluation.

Source: Bai et al.(2022), https://pmc.ncbi.nlm.nih.gov/articles/PMC9019245/ . Exact modeled-scaffold provenance is in ../crl2_ubiquitination_v1/validation/scaffold_preparation.json.
