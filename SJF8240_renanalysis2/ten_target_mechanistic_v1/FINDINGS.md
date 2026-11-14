# Ten-target mechanistic audit and frozen-model predictions

## Main findings

The matching input gap for SIK3/RPS6KB2 is now resolved. Both are scored above0.5 by every one of the60 frozen catalytic-feature models (five seeds × twelve matched folds). Mean catalytic scores: SIK3 0.809841; RPS6KB2 0.805951. The zero-padded control also predicts both above0.5 in all60 models. These are classifier outputs, not calibrated degradation probabilities or experimental labels. The model fails to reject known ABL1/AXL negatives, so the unknown calls have weak biological support.

The geometric discrimination problem is present in the measured signal itself: ABL1 can satisfy the frozen catalytic geometry while confirmed degraders such as MAPK14 and RIPK2 seldom do. Successful numerical checks do not make this a physically validated catalytic model.

- ABL1: bridge 79.8%; at least one lysine within20A 15.89%; exposed + within20A 15.83%; scaffold compatible 95.4%; all conditions including bridging 3.92%.
- AXL: bridge 86.2%; at least one lysine within20A 3.32%; exposed + within20A 3.32%; scaffold compatible 94.3%; all conditions including bridging 0.91%.
- MAPK14_full: bridge 94.5%; at least one lysine within20A 0.06%; exposed + within20A 0.06%; scaffold compatible 100.0%; all conditions including bridging 0.06%.
- RIPK2_full: bridge 70.8%; at least one lysine within20A 6.99%; exposed + within20A 6.67%; scaffold compatible 88.8%; all conditions including bridging 0.67%.

The eight preselected control poses reproduce every compared archived base/native feature array exactly (maximum error0). Independent geometry checks on12 preselected raw frames across ABL1,AXL,SIK3,RPS6KB2 passed using MDAnalysis unwrap,Bio.SVDSuperimposer and brute-force SciPy distances/overlap. These checks share the imaging convention; they do not validate the biological pose ensemble, catalytic-reference ensemble coverage, lysine rotamers or the cutoff choice.

## Step1 — Source and identity audit

Retained existing trajectories only; no docking, MD, minimization, model fitting, label changes or threshold tuning. Core141 poses reused; all33 unknown poses extracted (15 SIK3,18 RPS6KB2). Each has201frames over0–20ns. Topology hashes bind independently sequence-identified ELOB(chainA),ELOC(B),trueVHL(C),POI(D). The full name in a target directory does not imply a full-length protein: see modeled residue and lysine counts in the overview. Missing residues/lysines are outside this test, not negatives. Historical molecule provenance retained; supplier equivalence is not claimed.

### Protein coverage and numbering

| target       | uniprot   |   modeled_residues |   canonical_first |   canonical_last |   full_canonical_length |   observed_sequence_identity |   unmapped_modeled_residues |   sequence_mismatches |   ambiguous_positions |
|:-------------|:----------|-------------------:|------------------:|-----------------:|------------------------:|-----------------------------:|----------------------------:|----------------------:|----------------------:|
| ABL1         | P00519    |                263 |               235 |              498 |                    1130 |                     1        |                           0 |                     0 |                     0 |
| AXL          | P30530    |                270 |               521 |              813 |                     894 |                     1        |                           0 |                     0 |                     0 |
| DDR2         | Q16832    |                275 |               558 |              849 |                     855 |                     0.988889 |                           5 |                     3 |                     0 |
| EPHA2        | P29317    |                287 |               600 |              896 |                     976 |                     1        |                           0 |                     0 |                     0 |
| EPHB2        | P29323    |                248 |               615 |              893 |                    1055 |                     1        |                           0 |                     0 |                     0 |
| MAPK14_full  | Q16539    |                337 |                 4 |              352 |                     360 |                     0.988131 |                           0 |                     4 |                     1 |
| MET          | P08581    |                278 |              1053 |             1345 |                    1390 |                     1        |                           0 |                     0 |                     1 |
| RIPK2_full   | O43353    |                279 |                 8 |              316 |                     540 |                     1        |                           0 |                     0 |                     1 |
| RPS6KB2_full | Q9UBS0    |                241 |                92 |              332 |                     482 |                     1        |                           0 |                     0 |                     0 |
| SIK3_full    | Q9Y2K2    |                241 |                86 |              326 |                    1321 |                     1        |                           0 |                     0 |                     0 |

SIK3 is a241residue fragment of a1321residue canonical sequence; RPS6KB2 is241/482. The directory suffix_full does not establish complete modeled protein coverage. Lysine tables keep topology_resid/lysine_resid and add canonical_resid only for unambiguous sequence-matched lysines. Observed DDR2 and MAPK14 sequence mismatches and unaligned/tag residues are retained in construct_sequence_differences.tsv; they are not silently repaired. In particular a canonical DDR2 lysine maps to observed glutamate, so the modeled lysine inventory differs from the canonical sequence. Intent or experimental relevance of those sequence differences is not inferred.

## Step2 — Matching inputs

Original29 node features and contact edges follow the archived four-window implementation, including41 sampled frames at0.5ns,heavy-atom Shrake–Rupley SASA,legacy E2 ensemble and equal pose weighting. Five new lysine-node features use201frames at0.1ns and the exact archived native-reference definitions. Four graphs were prepared: native and zero-padded variants for each unknown,34features. Model architecture and trained weights were unchanged. Original dynamic graph.npz files were not rebuilt or overwritten.

## Step3 — Frozen prediction-only inference

| target       | variant             |   mean_score |   min_score |   max_score |   fixed_degrader_vote_fraction |   training_threshold_degrader_vote_fraction |   models |
|:-------------|:--------------------|-------------:|------------:|------------:|-------------------------------:|--------------------------------------------:|---------:|
| RPS6KB2_full | native_features     |     0.805951 |    0.630337 |    0.966777 |                              1 |                                           1 |       60 |
| RPS6KB2_full | zero_padded_control |     0.843072 |    0.667004 |    0.974038 |                              1 |                                           1 |       60 |
| SIK3_full    | native_features     |     0.809841 |    0.585869 |    0.980475 |                              1 |                                           1 |       60 |
| SIK3_full    | zero_padded_control |     0.839531 |    0.661495 |    0.9794   |                              1 |                                           1 |       60 |

120 existing checkpoints × two unknowns =240 predictions. Each was replayed with a fresh model instance; both original held-out scores per checkpoint were also reproduced (240 reference replays). Unknowns never entered training, threshold selection or AUROC calculations. Correlated folds and seeds do not create60 independent biological predictions. Thresholds come from the earlier training-only calibration; fixed0.5 and training-only decisions are both saved.

## Step4 — Mechanistic target comparison

| target       | group        |   poses |   modeled_POI_residues |   lysines |   bridge_fraction |   any_lysine_reachable20_fraction |   any_lysine_exposed_reachable20_fraction |   scaffold_compatible_fraction |   joint_bridge_productive_fraction |   median_VHL_anchor_RMSD_A |
|:-------------|:-------------|--------:|-----------------------:|----------:|------------------:|----------------------------------:|------------------------------------------:|-------------------------------:|-----------------------------------:|---------------------------:|
| ABL1         | non_degrader |      17 |                    263 |        18 |          0.798068 |                       0.158911    |                               0.158326    |                       0.954346 |                        0.0392157   |                    2.28663 |
| AXL          | non_degrader |      18 |                    270 |        16 |          0.861526 |                       0.0331675   |                               0.0331675   |                       0.94251  |                        0.00912106  |                    2.37959 |
| DDR2         | degrader     |      18 |                    275 |        14 |          0.740741 |                       0.0121614   |                               0.011885    |                       0.971255 |                        0.00138198  |                    2.28154 |
| EPHA2        | degrader     |      18 |                    287 |        20 |          0.702322 |                       0.0784964   |                               0.07822     |                       0.890547 |                        0.0550028   |                    2.21599 |
| EPHB2        | degrader     |      19 |                    248 |        13 |          0.779    |                       0.284629    |                               0.284368    |                       0.869076 |                        0.118356    |                    2.29674 |
| MAPK14_full  | degrader     |      18 |                    337 |        19 |          0.944997 |                       0.000552792 |                               0.000552792 |                       1        |                        0.000552792 |                    2.21025 |
| MET          | degrader     |      16 |                    278 |        19 |          0.894279 |                       0.198694    |                               0.196206    |                       0.92444  |                        0.11194     |                    2.395   |
| RIPK2_full   | degrader     |      17 |                    279 |        13 |          0.707638 |                       0.0699444   |                               0.0667252   |                       0.887621 |                        0.00673105  |                    2.19201 |
| SIK3_full    | unknown      |      15 |                    241 |        20 |          0.695522 |                       0.0361526   |                               0.0361526   |                       0.993698 |                        0.0222222   |                    2.3885  |
| RPS6KB2_full | unknown      |      18 |                    241 |        17 |          0.9356   |                       0.0257048   |                               0.0254284   |                       0.97319  |                        0.00690989  |                    2.25442 |

Fractions are means over frames within each pose, then equal means over poses. Any-lysine productive fraction is an OR over lysines per frame, not the mean lysine fraction used by the simple baseline. Productive here means distance<=20A, NZ sampled-surface exposure>=0.1 and whole-target scaffold overlap<=0.8A; joint adds simultaneous trueVHL/ligand and POI/ligand distances<5A. Deposited8R5H SY8_C3 is a crosslink surrogate, not a native Gly76 carbonyl. This is geometric access, not demonstrated ubiquitin transfer.

### Which gate removes access?

Compare reachable20 with exposed+reachable20, then scaffold compatibility and joint bridging. The scaffold-compatible fraction is high for the known targets, so blanket steric exclusion is not the sole cause of failure. ABL1 retains geometric access despite being a non-degrader; MAPK14 frequently bridges but almost never brings a modeled lysine near the chosen donor surrogate. Pose- and lysine-level tables identify the exact sites; this cannot be repaired by assuming every bridged binder degrades.

### Time and pose dependence

Four-window and per-pose results are saved separately with top five modeled lysines per target. Pose percentile ranges describe variation across retained docking starts; they are not biological confidence intervals. Consecutive-sample counts describe0.1ns sampled observations, not an experimentally validated catalytic residence time.

### Reference dependence

| target       |         min |         max |
|:-------------|------------:|------------:|
| ABL1         | 0           | 0.0392157   |
| AXL          | 0.000552792 | 0.0367606   |
| DDR2         | 0.00138198  | 0.0259812   |
| EPHA2        | 0.0163074   | 0.0892758   |
| EPHB2        | 0.0264467   | 0.12883     |
| MAPK14_full  | 0           | 0.000552792 |
| MET          | 0.00279851  | 0.11194     |
| RIPK2_full   | 0           | 0.0187299   |
| RPS6KB2_full | 0           | 0.0116086   |
| SIK3_full    | 0           | 0.0222222   |

Range spans the three archived reference constructions and VHL/ELOC alignment choices. These alternatives are sensitivity analyses, not independent physical replicates; no reference was selected against target labels. Median/max alignment RMSDs are recorded for each pose/reference.

## Step5 — Legacy50A signal checked separately

| target       |   modeled_lysines |   lysines_with_defined_delta |   positive_delta_fraction |   median_delta |
|:-------------|------------------:|-----------------------------:|--------------------------:|---------------:|
| ABL1         |                18 |                           18 |                 0         |              0 |
| AXL          |                16 |                           16 |                 0.0625    |              0 |
| DDR2         |                14 |                           14 |                 0         |              0 |
| EPHA2        |                20 |                           20 |                 0         |              0 |
| EPHB2        |                13 |                           13 |                 0.153846  |              0 |
| MAPK14_full  |                19 |                           19 |                 0.157895  |              0 |
| MET          |                19 |                           19 |                 0         |              0 |
| RIPK2_full   |                13 |                           13 |                 0.0769231 |              0 |
| RPS6KB2_full |                17 |                           17 |                 0         |              0 |
| SIK3_full    |                20 |                           20 |                 0.1       |              0 |

Source inspection of20_per_lysine_e2_distance.py corrects the earlier verbal explanation: original reach_50A is the fraction of frames with a lysine minimum distance to any E2 ensemble point<50A. It does not count the fraction of E2 points within50A. The table above reproduces that frame-access definition using corrected ligand-connected imaging and trueVHL bridging. Delta is bridge-state mean minus nonbridge-state mean within each pose, then equal averaging over poses with both states present. Undefined-state poses are excluded with counts retained. The corrected original eight-target result is reproduced; the original clean separation is not recovered.

A separate new sensitivity measure counts the fraction of663 E2 points within50A per frame; it is explicitly named ensemble_point_reach50 and saved in separate columns/tables. It is not presented as reproduction of the historical metric. Lysines/frames are not independent biological replicates; no inflated per-lysine p-value is used. Neither legacy50A measure is interchangeable with the native-reference20A catalytic feature.

## Step6 — Eight-labelled-target descriptive scores

| score                                   |   descriptive_eight_target_AUROC |   positives |   negatives | direction                       | unknowns_excluded   |
|:----------------------------------------|---------------------------------:|------------:|------------:|:--------------------------------|:--------------------|
| bridge_fraction                         |                         0.333333 |           6 |           2 | higher_is_degrader_prespecified | True                |
| any_lysine_reachable20_fraction         |                         0.5      |           6 |           2 | higher_is_degrader_prespecified | True                |
| any_lysine_exposed_reachable20_fraction |                         0.5      |           6 |           2 | higher_is_degrader_prespecified | True                |
| scaffold_compatible_fraction            |                         0.333333 |           6 |           2 | higher_is_degrader_prespecified | True                |
| productive_any_lysine_fraction          |                         0.5      |           6 |           2 | higher_is_degrader_prespecified | True                |
| joint_bridge_productive_fraction        |                         0.5      |           6 |           2 | higher_is_degrader_prespecified | True                |

SIK3/RPS6KB2 excluded from all performance metrics. EPHA2 remains label1. These are retrospective descriptive scores from eight independent proteins, not independent test-set validation. The previously completed five-seed GNN study remains unchanged (control meanAUROC0.733; catalytic0.667).

## Conclusion and remaining limits

The new audit supports numerical consistency, not a reliable degradation classifier. ABL1/AXL failures cannot be dismissed as a simple wrong-chain or serialization bug in these checked calculations. Ensemble/reference incompleteness, the historical molecule, missing protein regions, lysine sidechain sampling, coarse exposure/steric gates and biology beyond contact geometry remain possible causes. More tuning on eight labelled targets cannot resolve those limitations. Experimental degradation tests for SIK3/RPS6KB2 would provide new evidence; current predictions must not become labels.

## Saved material

Numbered scripts, source snapshots, manifests, unknown graphs, per-frame packets, pose/lysine/window tables, sensitivity tables, figures, CPU receipts and this Markdown/Word report are saved in this study folder. Original inputs and model checkpoints are referenced by hashes and remain unchanged.
