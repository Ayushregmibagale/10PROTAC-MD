# Method selection for six PROTACs across the kinome

Evidence review: 16 November 2026. Question: among proteins engaged by a particular PROTAC, which are degraded and which escape degradation? Constraints: six supplier-verified compounds; static structures only; reuse existing preparation/execution; later GNN; no new MD. The choices below are recommendations from source review, not measured superiority on our data.

## Recommended decision

Use the Bai productive-geometry question, experimentally grounded in the VHL catalytic work of Crowe, with COMPASS-style static CRL accessibility analysis on qualified docking ensembles. Keep our repaired execution and chemical-integrity framework. Retain the original Bai lab Rosetta implementation as a reference implementation to benchmark if our adapter remains unreliable.

For the primary small-data experiment, reuse our existing compact static graph encoder and add pose-ensemble pooling. Pocket-PROTACs is the best external starting point for the ensemble-pooling/training design among the reviewed degradation models. Its conformers are ligand conformers, so feeding entire ternary graphs requires an explicit new input branch; it is not a drop-in pretrained solution. DegradeMaster is a valuable alternative pocket-graph baseline. Do not assemble all published architectures into a single large model.

DeepTernary deserves a bounded pose-generation comparison before committing months of CPU time. It is a pretrained ternary-structure GNN, not a degradation classifier. Adopt it as a generator only if exact molecular connectivity, head placement, clashes, reference recovery and diversity are acceptable on our shared benchmark. A fast generator cannot repair a wrong kinase-bound warhead input.

This decision does not require completing the downstream GNN before Ernie docking starts. Store sufficient coordinates, mappings and rejection receipts so analysis can follow generation.

## What was actually checked

Accessible primary papers/methods, authors' repository documentation, and selected implementation files were inspected. Eight repositories were retrieved at pinned commits, recorded in SOURCE_MANIFEST.json. Source review focused on inputs, atom/chain assumptions, graph meaning, aggregation, dependencies and evaluation splits. No third-party training, complete pipeline replication or accuracy reproduction was performed. One AST-isolated COMPASS selector was exercised on a synthetic multi-chain counterexample; this is not a validation of the authors' full workflow. Our ongoing docking benchmarks remain separate.

The original six supplied publications were considered, plus Bai 2021/2022, the experimental VHL catalytic study, and four additional relevant computational comparators. This is a targeted review of plausible foundations, not a claim to have audited every PROTAC method.

## Scientific and implementation comparison

| Approach | Actual endpoint | Relevance to binder versus degrader | Role in our project |
|---|---|---|---|
| Bai 2021 Rosetta / PROTAC_ternary | Ternary formation and fraction of fully compatible models | Directly relevant foretinib/kinase/VHL examples; formation alone cannot establish degradation | Original Rosetta reference and linker-compatibility baseline |
| Bai 2022 CRL4A | Kinase ubiquitination and productive ensemble fraction | Closest scientific design: same broad kinase PROTAC, different kinase outcomes; CRBN machinery requires VHL adaptation | Main mechanistic question and fraction baseline |
| COMPASS | Linker/ternary feasibility with ubiquitination filters | Useful static VHL machinery; substantial inactive compounds remain feasible | Preferred reusable geometry-analysis foundation |
| PROTACMap | Productive conformation mapping and linker design | Good clustering/benchmark examples; full study includes MD and commercial docking | Selected analysis concepts and reference benchmark |
| Kumar & Sobhia | Full CRL2-VHL geometry, lysine proximity and dynamics | Same ligase; principal application PTP1B rather than kinome discrimination | Structural cross-check; not principal software foundation |
| Crowe 2024 | Experimental VHL ubiquitination mechanism | Direct evidence linking target orientation and lysine placement | Primary VHL reference/control source |
| Crowe 2025 homodimerizers | Experimental nonproductive VHL homo-PROTAC complexes | Establishes a possible failure mode, not the ABL1/AXL mechanism | Biological rationale, not prediction software |
| Integrative structural analysis 2026 | MD-based interpretation and CDK12 optimization | Mechanistic inspiration, but its main dynamics evidence is outside our scope | Context; no MD adoption |
| SE(3)-PROTACs | Degradation classification from ligand subgraphs and protein embeddings | Drug-conditioned output, but no explicit docked whole-ternary contact graph | Optional sequence/chemistry baseline |
| DegradeMaster | High/low degradation prediction using 3D ligand and protein-pocket graphs | Useful graph encoder and public labels; endpoint differs from experimentally established binder/nondegrader | External pocket-graph baseline / possible pretraining comparator |
| Pocket-PROTACs | Drug-target-E3 degradation classification with conformer pooling | Best reviewed external ensemble-model scaffold, but conformers describe the drug rather than entire ternary poses | Preferred ensemble-pooling/training reference |
| DeepTernary | Ternary coordinates | Addresses expensive pose generation, not degradation conditional on binding | Bounded acceleration benchmark |
| DegradoMap | Target/E3 degradability without drug identity | Cannot distinguish different PROTACs for the same target/E3 | Target-only baseline; not primary foundation |

### Bai 2021 and 2022: three different fractions must not be confused

Bai 2021's public Rosetta code includes a fully compatible fraction calculator; its examples include foretinib-bound kinases paired with VHL or CRBN. Bai 2022 instead classifies ternary ensembles using a CRL4A assembly and surface-lysine proximity to ubiquitin, then evaluates productive fractions against kinase ubiquitination assays. This is much closer to our target-selectivity question than comparing linkers for one target. The 2022 staged Rosetta sampling is substantially larger than our diagnostic runs. Its CRBN geometry and 16 Å criterion cannot simply be declared universal for VHL. These are related Rosetta workflows, not interchangeable with the London lab pRosettaC implementation.

Sources: [Bai 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8866032/), [code](https://github.com/karanicolaslab/PROTAC_ternary), [Bai 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9019245/).

### COMPASS: strongest directly reusable static geometry foundation

The authors combine RosettaDock, clustering, CRL clash/accessibility checks and fragment-based linker reconstruction. Degradation feasibility means at least one productive cluster supports reconstruction. The paper reports 93% recall with 65% false-positive rate, not strong general discrimination of all nondegraders. Crystal recovery was tested without the ubiquitination filter: reproducing a crystal arrangement and predicting catalytic competence are separate tasks.

Code inspection of ubiq_processor.py finds curated chain/residue inputs, a SER-based COM proxy, a fixed-radius line-of-sight test and a coarse minimum-distance clash gate. The example proxy must be traced to the supplied reference before treating it as a catalytic cysteine atom. A synthetic selector check shows equal lysine residue numbers on different chains overwrite each other. This is an adaptation requirement for our multichain inputs, not evidence their curated examples fail. Preserve outputs: several scripts remove intermediates.

Sources: [paper](https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cmdc.70385?af=R), [code](https://github.com/ssueron/compass).

### PROTACMap: useful components, not a complete open docking replacement

The paper benchmarks both VHL and CRBN, with 33/34 experimental structures recovered; this is not degradation classification accuracy. Its public files provide attachment filtering, clustering and distance analyses; actual docking/induced-fit stages use external commercial tools, with MD in the full workflow. Code has editable selections, fixed filenames/frame ranges and limited rejection logging. Clustering fits the combined proteins; this differs from target RMSD after fitting only VHL. Equal atom counts are checked, but our adaptation must also check residue correspondence. The scripts are research components requiring automation and data-contract checks for a large campaign.

Sources: [paper](https://link.springer.com/article/10.1186/s13321-026-01268-5), [code](https://github.com/Husam-PSE/PROTACMap).

### VHL mechanism and the two experimental references

Crowe 2024 directly studies the complete VHL catalytic machinery and BRD4 lysine ubiquitination, making it a stronger VHL reference than assuming a transplanted CRBN cutoff is correct. Its deposited trapped assembly contains engineered mutations; atom identity and modeled chemistry must be explicit. Crowe's 2025 homo-PROTAC work supports stable but nonproductive assemblies in a different system. Neither proves our ABL1/AXL explanation. Static ensembles can test that hypothesis, with uncertainty from catalytic reference conformations and missing full-protein context.

Sources: [Crowe 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11468923/), [authors' publication listing for 2025](https://sites.dundee.ac.uk/alessio-ciulli/publications/).

### Kumar & Sobhia and the MD framework

The Kumar/Sobhia primary abstract establishes PTP1B as the representative application, 6900-plus poses, CRL2-VHL assembly, clustering and dynamics. An accessible supporting-information excerpt identifies an E2 Cys93 reference, underlining that catalytic numbering is not interchangeable across E2 proteins. Full-text/supplement access was insufficient here to certify all detailed thresholds or released executable scripts; earlier BRD4-only descriptions were incorrect. The separate 2026 integrative study uses MD and experimental CDK12 optimization; its surface/lysine ideas are relevant but the full method conflicts with our static-only constraint.

Sources: [Kumar/Sobhia](https://pubmed.ncbi.nlm.nih.gov/41577441/), [supporting-information landing page](https://pubs.acs.org/doi/abs/10.1021/acs.jpcb.5c07053), [integrative study](https://www.sciencedirect.com/science/article/pii/S2666386426003267).

### GNN foundations: what each graph actually represents

**Pocket-PROTACs:** source model_protac.py encodes each ligand conformer, pools conformer embeddings by attention, and combines them with protein embeddings/pocket masks. It provides molecule-disjoint and additional cold-start evaluation designs. Reusing this pooling design is attractive, but whole-ternary pose graphs, explicit target–VHL contacts and E2 geometry are new inputs. Protein pocket attention is not proof of modeled interprotein contact or ubiquitination. Data preprocessing, sample identity and checkpoint compatibility still need independent reproduction.

**DegradeMaster:** source model.py encodes PROTAC, target pocket and E3 pocket in separate branches before fusion. It is useful for 3D pocket graphs, but does not directly consume our complete docked assembly. Its main reported evaluation is random splitting with duplicate target/SMILES exclusion, not a guarantee of unseen-kinase generalization. Its high/low activity labels and pseudolabeling objective do not establish that every negative is a confirmed binder. Reproduce a supervised baseline before considering pseudolabeling.

**SE(3)-PROTACs:** model.py uses ligand subgraphs and ESM protein embeddings, without explicit ternary pose contacts or CRL geometry. It is a reasonable chemistry/sequence comparator, but adapting it to our proposed mechanistic input is substantial. A repository checkpoint existing does not establish applicability to our six exact structures or held-out proteins.

**DegradoMap:** source and author README describe a drug-independent target/E3 model with lysine-weighted pooling. It cannot distinguish SJF8240 from another compound against the same kinase/VHL. The README also reports weaker broader-seed and VHL target-unseen results than headline best-seed metrics; use it as a context-only comparator, not our core predictor. External ubiquitination annotations must not encode our held-out experimental outcome.

Sources: [Pocket-PROTACs paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13501292/) and [code](https://github.com/Adochew/Pocket-PROTACs), [DegradeMaster paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12261415/) and [code](https://github.com/ABILiLab/DegradeMaster), [SE(3)-PROTACs](https://pubmed.ncbi.nlm.nih.gov/42108630/) and [code](https://github.com/drugparadigm/SE3-protacs), [DegradoMap](https://arxiv.org/abs/2606.04021) and [author code/results](https://github.com/bryanc5864/DegradoMap).

### DeepTernary: plausible acceleration, different learning objective

The pretrained model predicts ternary structures and supplies downloadable weights; its public inference interface still needs two proteins, exact ligand/warhead mappings and multiple conformer seeds. The quoted per-prediction speed is not end-to-end campaign timing. It omits our full CRL context unless separately evaluated. We should compare it on identical inputs and controls, keeping geometry/chemical checks independent of its confidence. Do not label its structures degradation predictions or assume it solves the existing pocket clashes.

Sources: [paper](https://www.nature.com/articles/s41467-025-61272-5), [code/inference instructions](https://github.com/youqingxiaozhua/DeepTernary).

## Concrete primary model and baselines

The prediction unit is one compound–target–assay context, with a bag of qualified ternary poses. Do not train each pose as an independent experiment. A shared compact encoder produces a representation of each pose; an invariant mean/max or attention pool combines pose representations into one prediction. Start with mean pooling and a small encoder before adding attention. Use the existing static EnsembleEncoder as an implementation starting point, checking its graph/schema compatibility; remove the MD branch explicitly rather than feeding fabricated dynamics.

Nodes should retain protein/component identity, amino-acid identity for residue nodes, exact ligand atom chemistry, and explicit lysine flags. Edges should distinguish intraprotein proximity, ligand bonds and measured intermolecular contacts. Include continuous lysine exposure/burial and distances/orientation to a checked catalytic reference as declared features or a separate feature branch. A global distance or clash score need not be attached to every node. ELOB/ELOC can provide context, but VHL identity determines the biological interaction arm. Avoid hardcoded chain letters.

For each target, pool a standardized set of diverse clusters, recording pose sampling/rejection rates. Multiple poses are plausible configurations, not independent labels or equilibrium samples. If many duplicate poses inflate a cluster, compare equal-pose and equal-cluster summaries. A learned attention weight is not a thermodynamic probability or a causal mechanism.

Required comparisons, using identical outer splits:

1. Binding/contact or ternary-feasibility baseline.
2. COMPASS-style any-compatible-cluster rule.
3. Bai-style productive fraction, with an explicit denominator and catalytic-reference ensemble.
4. Logistic regression or a small tree model on geometry summaries.
5. Existing compact static GNN without catalytic features.
6. Same GNN with catalytic features and ensemble pooling.
7. Optional frozen/public-model sequence/chemistry baseline, screened for train/test overlap.

Keep the analysis finite: these ablations answer whether the graph and catalytic features add information. Running a large model zoo is unnecessary for the initial decision.

## Productive fraction and filtering contract

Ternary formation compatibility, catalytic productivity and successful cellular degradation are distinct endpoints. Bai 2021 fully compatible fraction, Bai 2022 productive fraction and COMPASS binary feasibility must have separate column names. Report all search, chemical-QC and geometric failure counts.

A useful proposed score is the weighted mean over chemically valid ternary poses of the fraction of predefined valid catalytic references allowing lysine access without prohibited clashes. Explicitly retain catalytically incompatible poses as zeros in that denominator; deleting them first would inflate the productive fraction. Treat absent/malformed inputs or failed searches as missing/failed, not experimental nondegraders. Define reference validity independently of our eight labels. Report equal-pose and equal-cluster sensitivity because computational sampling is not an equilibrium distribution.

Bai's Ub-terminal distance, PROTACMap's C-alpha to C-alpha distance and COMPASS's lysine-NZ/proxy geometry measure different quantities. Do not interchange their 16, 50 and 60 Å thresholds. Save E2 identity, residue mapping, atom names and reference provenance. Continuous values plus prespecified sensitivity analyses are safer than tuning one cutoff to separate ABL1/AXL.

Full-chain context matters for clashes and lysines. A kinase-domain screen may generate useful candidates but cannot establish full-length ubiquitination access. Report construct coverage and perform static full-protein context checks on prioritized candidates when reliable structures exist.

## Labels and evaluation determine whether the result answers our question

Our intended endpoint is degradation conditional on engagement: distinguish binders that degrade from binders that do not under specified experimental conditions. Many public PROTAC datasets instead label high versus low degradation without establishing binding in every negative. Public pretraining is potentially useful, but does not replace this conditional validation cohort.

The present SJF8240 set has six designated degraders and two nondegraders. SIK3/RPS6KB2 remain unknown. Do not create confirmed labels from hypotheses, docking failures, predicted productive fractions or GNN outputs. Audit exact compound connectivity, binding evidence, cell line, concentration and time for all six compounds before pooling labels. An unmeasured kinase is not a negative. Thousands of poses or thousands of unlabelled drug–kinase combinations do not increase independent labelled sample size.

Use target-grouped outer evaluation, grouping variants/poses across compounds. Where enough data exist, test kinase-family and held-out-compound generalization separately; six compounds alone provide limited compound-level uncertainty. Fit normalization, thresholds, class weights and hyperparameters only within training/inner folds. Eight labels support exploratory LOPO comparisons, not strong generalization claims. AUROC cannot be calculated for a one-target held-out fold alone: report pooled held-out predictions and per-target scores, with sensitivity/specificity and uncertainty. A GNN improvement must be measured against simpler models on the same held-out samples.

Class-weight note: with label1=degrader, BCEWithLogitsLoss pos_weight weights positives. A value3 does not increase the penalty on negative examples. A conventional training-fold balance choice is n_negative/n_positive, not its reciprocal; use the current frozen experiment's declared loss when reproducing historical results, and version any corrected-loss experiment separately.

## Work order and decision gates

1. Complete the already-running centroid-repair reference benchmark; independently assess both reference recovery and ranking, chemistry, head engagement and sterics. Preserve failed baseline and repaired results.
2. Check original Bai/standard Rosetta docking as a reference implementation if the adapter remains weak. Do not claim that adopting a published package transfers its validation to our modified inputs.
3. Run a bounded pretrained DeepTernary comparison if runtime/checkpoint/input compatibility can be established without delaying the existing path. Select a generator by common-control quality, diversity and measured end-to-end cost, not paper headline speed.
4. Map reusable COMPASS CRL/accessibility components to our verified VHL reference; add chain-aware residues, explicit catalytic atoms, immutable inputs, failure receipts and independent geometry checks. Use PROTACMap analysis as a comparison, not a second mandatory full campaign.
5. Qualify six compound inputs and a small Ernie technical/resource pilot, then freeze full-campaign sampling and targets. Obtain final approval for the concrete remote submission. No GNN-training dependency is required for docking submission.
6. Build pose graphs from saved outputs, run simple geometry baselines and the compact ensemble GNN. Borrow Pocket-PROTACs pooling/training organization where useful; use DegradeMaster/SE3 as optional external baselines rather than wholesale replacement.

No production switch, new remote launch, third-party installation or MD was performed for this review. No method is presently certified to separate our binders/nondegraders. The recommended route maximizes reuse while preserving a testable comparison and our existing campaign work.

## Latest local benchmark evidence during this review

The explicit-centroid repair completed both references and passed independent source/input/output integrity checks. MZ1:1/129retained poses passes the frozen coarse gate, best4.4000A, passing rank16. GSK215:0/3passes, best35.1094A, zero native contacts. Thus the repaired adapter is NOT qualified as a general production generator. Do not deploy it because one reference improved, and do not revert to the known-defective representation merely to recover a better score. This strengthens the recommendation to compare a standard/original Rosetta implementation and a bounded pretrained DeepTernary generator on common inputs before choosing the campaign generator. Full findings are in ../centroid_repair_results_complete/. This review has not launched those new comparisons.

Selected source snapshots are preserved in a private archive outside Git; SOURCE_ARCHIVE_RECEIPT.json records its path/hash. Only original review text and source provenance are committed publicly.
