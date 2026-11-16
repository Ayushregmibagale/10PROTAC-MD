# Papers, reusable methods and campaign decisions — 16 November 2026

This consolidates today's discussion and the earlier papers brought into it. It records our source-review conclusions and proposed uses, not reproduced performance. Detailed evidence, pinned source hashes and access limitations remain in [the method review](method_selection_review/RECOMMENDATION.md), [source manifest](method_selection_review/SOURCE_MANIFEST.json), [P4ward review](method_selection_review/P4WARD/REVIEW.md) and [PotentialLinker review](method_selection_review/POTENTIALLINKER_REVIEW.md).

## Our question and fixed scope

For each purchased PROTAC, distinguish engaged kinases that degrade from engaged kinases that do not. Generate static ternary ensembles for six Tocris/Bio-Techne compounds: GSK215, SJF8240, SJF1521, JNJ1013, SJF_delta and MS143. No new MD. MS143 is an ATP-competitive AKT degrader, not an allosteric compound. MZ1 is an external structural control, not a seventh screening drug.

The immediate priority is qualifying the docking generator and preparing the Ernie pilot. Downstream full-CRL analysis and GNN development must not become additional launch prerequisites. No publication or software package automatically validates our adapted pipeline.

## What each paper can contribute

| Paper / approach | What it addresses | What we can use | Qualification or limitation |
|---|---|---|---|
| **Bai 2021**, Rosetta ternary modeling | Ternary formation, linker compatibility and fully compatible model fraction; includes kinase/foretinib examples | Original Rosetta reference implementation and formation/linker baseline | Different from the London pRosettaC implementation. Compatibility does not establish ubiquitination or degradation. |
| **Bai 2022**, modeling CRL4A | Kinase ubiquitination using productive ternary geometry in larger CRL machinery | Main scientific template: align ternary ensembles into catalytic assemblies and calculate productive fractions | CRBN/CRL4A requires checked VHL/CRL2 adaptation. Its approximately 16 Å criterion is not universal. |
| **COMPASS**, ChemMedChem 2026 | Protein docking, clustering, ubiquitination-accessibility filters and linker feasibility | Static geometry baseline; clash, orientation/accessibility and compatible-cluster logic; selected code after adaptation | Reported 93% recall accompanies 65% false-positive rate. Source has proxy/mapping assumptions requiring explicit checks. Not an established solution to our binder/nondegrader problem. |
| **PROTACMap**, J. Cheminformatics 2026 | Maps productive arrangements and compatible linker windows | Clustering, attachment-distance logic, structural benchmark and linker-compatibility comparison | Study covers VHL and CRBN. Public scripts are components, not a complete open docking replacement. Full workflow uses commercial tools and MD; we need not reproduce those stages. |
| **Kumar & Sobhia**, JPCB 2026 | Full CRL2-VHL geometry, lysine positioning, clustering and dynamics | VHL assembly/reference cross-check and mechanistic comparison | Main representative application is PTP1B; earlier BRD4-only description was incorrect. Incomplete full-method access means exact thresholds and executable implementation remain uncertified here. |
| **Crowe et al. 2024**, experimental VHL catalytic assembly | Target orientation and lysine ubiquitination in CRL2-VHL | Primary VHL structural reference and independent catalytic controls | Engineered trapped assemblies, mutations and missing atoms must be mapped explicitly. |
| **Ciulli/Crowe 2025 VHL homodimerizers** | Stable but nonproductive homo-PROTAC assemblies | Experimental rationale that ternary stability can coexist with failed ubiquitination | Different targets/system. Does not establish the mechanism for SJF8240–ABL1/AXL or prove our novelty. |
| **Integrative structural analysis**, Cell Reports Physical Science 2026 | MD, lysine topology and CDK12 degrader optimization | Motivation for accessibility and lysine-context features | Dynamics results cannot be claimed from our static data; no new MD adopted. |
| **P4ward**, JCIM 2025 | Automated ternary modeling plus larger-CRL clash and lysine-access filters | Strong candidate implementation for the post-docking full-CRL analysis layer, reusing our saved poses | Selected source and bundled models checked; not installed or validated end to end. One VHL reference, fixed Ub point and geometry assumptions need repair/checking. |
| **PotentialLinker**, J. Med. Chem., online 16 November 2026 | Reinforcement learning for linker optimization with NWHM degradation score; experimental FGFR4 D11 discovery | Candidate scoring-function comparator or later structural features, conditional on full methods/code review | Only publisher-deposited abstract verified. NWHM components, full-CRL construction, atom/cutoff definitions and author code not verified. New-linker generation is outside our fixed-compound screen. |
| **Pocket-PROTACs** | Drug–target–E3 prediction with ligand-conformer ensemble pooling | Best reviewed external reference for pooling/training organization; retain a compact model | Ligand conformers are not whole ternary-pose graphs. Explicit ternary contacts and catalytic features would be our new inputs. |
| **DegradeMaster** | 3D ligand and separate target/E3 pocket branches for high/low degradation | Pocket-graph baseline and possible encoder/pretraining comparator | Does not directly consume our full docked ternary graph. Low degradation does not necessarily mean a confirmed binder/nondegrader. |
| **SE(3)-PROTACs** | Chemical subgraphs and protein embeddings for degradation prediction | Optional chemistry/sequence baseline | Inspected implementation is not a graph of explicit whole-ternary contacts or full catalytic machinery. |
| **DeepTernary** | Pretrained ternary coordinate prediction | Possible bounded generation speed/quality comparator if later needed for launch-critical recovery | Structural generator, not degradation classifier. No new comparison project is currently authorized by our priority decision. |
| **DegradoMap** | Target/E3 degradability without compound identity | Optional target-only baseline | Cannot distinguish our six compounds against the same target/E3. Not a primary drug-specific predictor. |

## Corrections to preserve

1. Bai's own Rosetta workflow and pRosettaC are distinct. Our current engine is an adapted random-global/PyRosetta/fixed-head linker workflow, not an exact reproduction of either paper.
2. PROTACMap's 33/34 structure recovery (about 97%) and P4ward's reported 76.5% structural hit rate are not degradation-classification accuracies.
3. A 16 Å Ub-terminal criterion, a 50 Å Cα proxy and a 60 Å accessibility rule refer to different geometry. They are not interchangeable or universal. Lysine NZ is the nucleophile; ubiquitin Gly76's terminal carbonyl participates in transfer. Gly76 Cα is a positional proxy, not the attacking atom.
4. Using a scaffold/arm template from another cullin requires explicit structural validation; it is not automatically invalid merely because the template contains CUL1.
5. The P4ward VHL distance is to a stored Ub reference point, not E2 Cys. Its source computes target-alone SASA, not E3-induced lysine burial. Its terminal-point mapping, chain-aware selection, missing NZ handling, hydrogen selection and filter toggles need checking before reuse.
6. Lysine burial is a testable feature, not an established explanation for ABL1/AXL. Earlier quoted correlations or extremely small per-lysine p-values must not be treated as independent target-level validation without auditing provenance and nesting.
7. With label 1 = degrader, BCEWithLogitsLoss `pos_weight=3` increases positive-class weighting, not the negative-class penalty. Preserve historical reproduction settings; version any correction separately and derive training weights within folds.

## Retained pipeline, in order

1. **Qualify inputs and generation.** Exact supplier connectivity/stereochemistry, appropriate kinase-domain/head placement, preserved residues/PTMs, linker feasibility, independent sterics, known-structure recovery and ranking. Preserve raw failures and hashes. Sampling counts alone do not establish convergence.
2. **Run a bounded six-drug Ernie technical/resource pilot when qualified and approved.** Freeze targets, sampling, seeds, resource limits and resumability before full-scale production. Do not call a failed search a biological negative. Actual prepared coverage must be reported rather than assuming 600 usable kinases.
3. **Reuse saved ternary poses for full-CRL analysis.** Current docking contains target + PROTAC + VCB (VHL/ELOB/ELOC). Separately align through VHL to validated CRL2/E2/Ub references; assess larger-machinery clashes, exposed lysines, burial, orientation and continuous distances. P4ward and COMPASS are candidates for this layer; Bai supplies the productive-fraction framework. See [full-CRL memory](FULL_CRL2_ANALYSIS_MEMORY.md).
4. **Compare finite baselines.** Contacts/ternary feasibility; any-compatible-cluster rule; productive fraction; simple regression/tree on geometry; compact static GNN; same GNN with catalytic features. PROTACMap supplies useful structural/clustering comparisons. PotentialLinker scoring is optional only after method verification.
5. **Pool poses at compound–target level.** Begin with compact shared encoding and simple pooling, borrowing Pocket-PROTACs design where appropriate. All poses/construct variants from a target stay grouped across evaluation splits. Attention weights and docking frequencies are not equilibrium probabilities.

Productive-fraction denominators must retain chemically valid but catalytically incompatible poses as zeros. Malformed inputs and absent searches are failures/missing, not negative experimental labels. Record rejection counts; compare equal-pose and equal-cluster summaries and catalytic-reference sensitivity. Domain-only geometry does not establish full-length lysine accessibility.

## Experimental labels and publishability

SJF8240 user-designated degraders: DDR2, EPHA2, EPHB2, MAPK14_full, MET and RIPK2_full. Nondegraders: ABL1 and AXL. SIK3_full and RPS6KB2_full remain unknown, including when explored as hypothetical negatives. Do not transfer SJF8240 labels to another drug or the historical wrong-connectivity branch. Confirm binding, assay context, dose and time before combining datasets.

Eight target labels support exploratory held-out comparisons, not strong generalization claims. Thousands of poses and unmeasured kinases do not increase independent labelled sample size. Publication strength must come from qualified structures, honest baselines, independent labels, grouped evaluation and clear limitations; it cannot be inferred from choosing a particular paper or reproducing a headline structural hit rate.

## Campaign snapshot checked while saving this note

The local soft-anchor repair qualification has finished both references at 50 global × 10 local × up to 20 conformers, seed 42. GSK215: 11/500 linker constructions, 11 engine-QC retained. MZ1: 480/500 constructions, 392 retained. Independent source/input/output integrity checks PASS for both. These are not structural recovery results. Frozen reference assessment is complete: GSK215 1/11 passes (best target CA RMSD after VHL fit 3.71393 Å; top1 fails, top10 contains a pass); MZ1 0/392 passes (best RMSD 6.14650 Å). The gate remains RMSD ≤5 Å AND native residue-contact fraction ≥0.5, plus inherited chemical/gross-clash QC. Generator qualification across both references therefore still fails. See anchor_repair_results_complete/ for compact receipts. Raw coordinates currently remain at the completed local run root; durable full-raw archiving is still pending. The staged six-drug Ernie package remains on HOLD; full production has not been submitted by this work. Do not relaunch these completed qualification tasks.

## Primary sources and author code

- Bai 2021: https://pmc.ncbi.nlm.nih.gov/articles/PMC8866032/ ; https://github.com/karanicolaslab/PROTAC_ternary
- Bai 2022: https://pmc.ncbi.nlm.nih.gov/articles/PMC9019245/
- COMPASS: https://doi.org/10.1002/cmdc.70385 ; https://github.com/ssueron/compass
- PROTACMap: https://doi.org/10.1186/s13321-026-01268-5 ; https://github.com/Husam-PSE/PROTACMap
- Kumar & Sobhia: https://doi.org/10.1021/acs.jpcb.5c07053
- Experimental VHL assembly: https://pmc.ncbi.nlm.nih.gov/articles/PMC11468923/
- VHL homodimerizers: https://doi.org/10.1101/2025.11.18.689129
- Integrative study: https://www.sciencedirect.com/science/article/pii/S2666386426003267
- P4ward: https://doi.org/10.1021/acs.jcim.5c00614 ; https://github.com/SKTeamLab/P4ward ; https://skteamlab.github.io/P4ward/
- PotentialLinker: https://doi.org/10.1021/acs.jmedchem.6c02194 ; abstract metadata https://api.crossref.org/works/10.1021/acs.jmedchem.6c02194
- Pocket-PROTACs: https://pmc.ncbi.nlm.nih.gov/articles/PMC13501292/ ; https://github.com/Adochew/Pocket-PROTACs
- DegradeMaster: https://pmc.ncbi.nlm.nih.gov/articles/PMC12261415/ ; https://github.com/ABILiLab/DegradeMaster
- SE(3)-PROTACs: https://pubmed.ncbi.nlm.nih.gov/42108630/ ; https://github.com/drugparadigm/SE3-protacs
- DeepTernary: https://doi.org/10.1038/s41467-025-61272-5 ; https://github.com/youqingxiaozhua/DeepTernary
- DegradoMap: https://arxiv.org/abs/2606.04021 ; https://github.com/bryanc5864/DegradoMap

This document consolidates existing source checks; it does not claim all papers were fully accessible or all public code was executed. In particular, PotentialLinker remains abstract-only and Kumar/Sobhia access was incomplete. Ernie-first scope remains in force.
