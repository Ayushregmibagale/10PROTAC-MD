# Six-PROTAC static ensemble study: detailed plan and current position

Prepared 15 November 2026 from a read-only review of the local 10PROTAC-MD project.

## 1. Current decision and scope

Use chemically verified pRosettaC ternary ensembles as the main structural dataset. New MD is not a required stage. Existing MD remains historical/supporting evidence and is not mixed into the primary static feature schema. Develop a static ensemble GNN only after structural quality and experimental labels are established. This document records a plan; it does not launch docking, MD, training or remote jobs. The old supplier-correct MD pilot remains paused.

Primary question: among proteins with evidence of binding to the relevant compound, do static ternary ensemble geometry and lysine accessibility help distinguish experimentally degraded proteins from experimentally binding non-degraders? For unmeasured targets, outputs are structural predictions, not confirmed binding or degradation.

Six compounds are now registered with supplier-verified molecular identities. The user confirmed all six were purchased from Tocris/Bio-Techne on 15 November 2026. Bound-head preparation and docking qualification remain incomplete for MS143. Approximately 600 kinase candidates is an intended scope, not a verified current number of usable proteins. Count UniProt proteins, domains, structure variants and compound-target pairs separately.

## 2. What we actually have

| Compound | Verified supplier identity in registry | Local structural starting point | Current interpretation |
|---|---|---|---|
| GSK215 | Tocris 7818 | Native ternary reference 7PI4; PTK2/FAK control | Best existing native control; supplier graph verified, campaign recovery still needs evaluation |
| SJF8240 | Tocris 7266 | MET–foretinib parent reference 3LQ8/88Z; supplier-correct preparation branch | Exact complete molecule verified; a parent ligand structure is not an observed full-PROTAC ternary structure |
| SJF1521 | Tocris 7261 | EGFR–lapatinib parent reference 1XKK/FMM; modeled head candidates | Supplier graph verified; modified head placement needs retained validation |
| JNJ1013 | Tocris 8029 | IRAK1 parent reference 6BFN/DL1; modeled head candidates | Supplier graph verified; head placement is not a measured ternary pose |
| SJF_delta / SJFδ | Tocris 7267 | Foretinib-parent references and separately recorded MAPK13/p38δ overlay control | Distinct chemistry/exit vector must be preserved; do not copy SJF8240 parameters |
| MS143 | Tocris 8081; CAS 2376137-41-4 | Supplier graph registered; AKT1–AZD5363 4GV1 is a parent-reference candidate | Molecular identity verified; bound-head maps, parameters, pocket transfer and docking remain pending |

Important identity correction: the original MS143 publication identifies MS143 (compound 20) as a VHL-recruiting AKT degrader derived from the ATP-competitive inhibitor AZD5363. It is not described there as the allosteric AKT degrader. AKT1–AZD5363 structure 4GV1 is a parent-reference candidate, not an MS143 ternary structure. If an allosteric compound was intended, resolve its exact identity before substituting a structure. Do not silently rename MS15 or another compound MS143.

Evidence: https://pmc.ncbi.nlm.nih.gov/articles/PMC8900464/ . Supplier identities and reference IDs above come from config/molecules.json and molecules/verification.json; they have not been freshly re-downloaded in this planning review.

### Existing preparation and execution evidence

- Six supplier molecular graphs are now verified. Existing compound-specific head/anchor preparations cover the original five; MS143 preparation is pending. Molecular identity verification is distinct from reliable bound placement and successful docking.
- Saved preparation summaries report 770 domains across 435 UniProt proteins and 496 exact KLIFS fingerprints. The original CE preparation has 112 compound-target candidates; a separate explicit-pocket comparison has 161. These are alternative preparation cohorts, not 273 unique pairs or complete kinome coverage. Recompute current inventories before freezing a campaign.
- The frozen small five-compound pilot produced two QC-valid SJF1521 poses and three SJFδ poses; GSK215, SJF8240 and JNJ1013 had no valid pose at that small search budget. These are sampling outcomes, not biological negatives. A native-frame GSK215 control separately produced one valid pose at a 100-conformer budget; it is not independent global recovery.
- The later paused supplier pilot reports 30 Boltz reference models with inconsistent recovery, a failed native MZ1 single-pose QC attempt, MET domain placement passing while EPHA2/ABL1/AXL placements failed retained checks, and incomplete global-recovery qualification.
- Raw saved global-search summaries now contain 123 finite-score MZ1 models in 114 clusters and three GSK215 models in two clusters. These raw counts are not QC-certified native recovery. Both report zero clusters with at least five members. Assess controller receipts, full-complex QC and native alignment before interpreting them. Do not automatically restart completed raw searches because the stopped STATUS.json contains stale 'active' stage descriptions.
- No completed six-compound production screen is documented by the reviewed receipts. MS143 is integrated at the molecular-registry level only; docking inputs are not yet prepared. Ernie campaign release remains pending in the saved root documentation.
- Historical SJF8240 chemistry differs in connectivity from the supplier molecule. Archived kinome poses also have placement/geometry failures. The bounded repair report passes only FER, TGFBR1 and PRKG1 under archived-model operational checks; it does not certify them as supplier-SJF8240 binders or validate a top ten.
- Corrected native-VHL mechanistic/GNN analyses are useful historical diagnostics. Their eight-target results are not a validated classifier for these six supplier-correct compounds. Do not use the best seed AUROC as campaign performance.

## 3. Phase A — freeze compounds, endpoints and target inventory

Deliverables: compounds.tsv, labels.tsv, target_structures.tsv, coverage.tsv and a machine-readable source/hash manifest.

1. Reuse the five verified supplier graphs without altering their identities. MS143 catalog/CAS/isomeric SMILES/InChIKey are now recorded from Tocris and pass the existing identity verifier. Prepare bound heads from its exact graph and compare attachment/stereochemistry to primary material. Record the warhead, linker, VHL head and attachment atoms explicitly.
2. Separate binary reference crystals, observed ternary crystals, transferred heads and modeled candidates. Store atom correspondence by chemical graph/CCD atom names, never by coordinate order alone. Validate modified conjugated heads rather than assuming a parent inhibitor is chemically identical.
3. Build compound-specific labels with source, dose, time, cell line, binding evidence, degradation endpoint and uncertainty. Confirmed binding non-degraders require appropriate binding and negative degradation evidence in a stated context. A missing assay, failed docking or below-detection protein is not a negative.
4. Preserve the current SJF8240 panel: DDR2, EPHA2, EPHB2, MAPK14_full, MET and RIPK2_full positive; ABL1 and AXL negative, subject to recording the underlying experimental evidence. SIK3_full and RPS6KB2_full remain unknown. Do not carry these labels onto another PROTAC.
5. Resolve structure/domain boundaries, sequence identities, pocket coverage, missing residues and mutations. A kinase-domain graph cannot certify access to lysines elsewhere in the full protein. Exclude missing lysines from the measured denominator and explicitly flag incomplete constructs.
6. Replace the '600' shorthand with verified totals at each level. Record why every candidate is included, deferred or excluded. Binding evidence and structural feasibility remain separate fields.

Exit condition: the six compound identities are unambiguous; every structure has a provenance/selection record; labels are compound- and assay-specific. Five-compound work can progress independently if the sixth remains unresolved.

## 4. Phase B — qualify pRosettaC before scaling

Deliverables: frozen_benchmark.json, recovery.tsv, rejection_counts.tsv and runtime_storage.tsv.

1. Audit the stopped MZ1/5T35 and GSK215/7PI4 runs before rerunning anything. Reuse valid outputs with matching source/input hashes. Record interrupted, failed, no-valid-pose and complete outcomes separately.
2. Separate native-frame integration controls from global recovery using independent bound binaries. Evaluate chemistry, head retention, all-component clashes, interface/ligand RMSD and score-ranked top-1/top-10 recovery plus best sampled recovery. Define criteria before evaluating the next benchmark.
3. Use at least three declared docking seeds (proposed 42, 7, 123) on the small reference panel. Track attempted/accepted global poses, local models, linker closures, scores and final QC counts. Do not silently drop failed attempts.
4. Use PROTACMap as a structural reference: assess orientation diversity and recovery against experimental structures. Its reported 33/34 structural recovery is not a degradation-accuracy target. Compare matched inputs/metrics where possible; use label-blind held-out references when selecting thresholds. Do not assume its commercial MOE steps are available.
5. Resolve failed head transfers with suitable structures or separately validated placement methods. Do not loosen pocket checks until all targets appear usable. For a genuinely allosteric compound, use the correct allosteric pocket and required domain context instead of an ATP-pocket transfer.
6. Measure CPU-hours, peak RAM and raw storage on representative easy/difficult pairs. Choose production sampling budgets from measured diversity/recovery and cost. No fixed large pose count is assumed sufficient, and finite failure is not evidence of non-binding.

Exit condition: documented recovery and reproducibility are adequate for the stated structural claim, and the complete QC/rejection path works. If recovery fails, finish the failure diagnosis before a kinome campaign.

## 5. Phase C — bounded six-compound pilot, then staged kinome docking

Pilot: one experimentally supported target/control per registered compound plus selected binding non-degraders where evidence exists. Include the SJF8240 four-target subset MET/EPHA2/ABL1/AXL after their structural placement issues are resolved. Select controls by evidence, not by whether a preliminary productive score looks favorable.

Use identical documented search policies across classes. Introduce additional seeds/budgets according to predefined convergence rules. Cluster poses using a frozen ligand-independent protein-orientation/interface metric with correct chain mapping. Do not select graph inputs using the desired label or only a favorable productive-geometry score.

Inspect feature/rank stability over nested sample sizes and independent seeds. Preserve per-seed denominators. A computational pose count or cluster size is not a thermodynamic population. Report sampling-weighted and equal-cluster summaries separately.

After the pilot passes, process validated compound-target pairs in resumable batches. Store raw poses, scores, cluster assignments, reject reasons, seeds, versions, hashes and terminal receipts. Budget storage from the measured revised campaign; old MD-inclusive storage projections do not define the new static-only requirement. Failure/timeout/missing pocket remain explicit statuses, never biological zeroes.

## 6. Phase D — frozen static feature extraction

Deliverables: per_pose_features.parquet, per_lysine_features.parquet, ensemble_features.tsv, selections.json and schema.json.

- Dual-head contacts: verify the chemically defined target head engages the intended pocket and VHL head engages actual VHL. Include contact counts/chemistry, not only any ligand atom within 5 Å. ELOB/ELOC are accessory proteins, not interchangeable VHL selections.
- Interface descriptors: protein-protein contacts, interface area, polar/nonpolar interactions, clashes and available Rosetta energies. Keep energy terms and normalization explicit; docking energies are not measured affinities.
- Lysine exposure: whole-residue relative SASA and separately NZ SASA. Calculate the target alone and in the same fixed conformation within the ternary complex, with identical radii/probe settings. Their difference measures static occlusion upon adding partners, not MD-induced conformational burial. State whether ligand/accessory chains are included.
- Productive accessibility: align a validated native CRL2–VHL–E2–Ub reference ensemble through the correct VHL frame; assess machinery clashes, exposed target lysines and declared distances. Keep E2-Cys, Ub-Gly76-Cα and Ub carbonyl-carbon metrics distinct. Reconstructed missing donor atoms are modeled, not experimental coordinates. Test reference uncertainty independently of class labels.
- Ensemble summaries: productive fraction among eligible pose/reference combinations, cluster-balanced productive fraction, accessible lysine counts, feature quantiles and variability across seeds/reference assemblies. Keep all denominators and missing-reference counts.
- Sensitivity: freeze a small physically justified range of contact/exposure/distance definitions using references before outcome testing. Do not choose the cutoff that best separates eight targets and present it as an independent result.

Bai motivates productive fractions; PROTACMap motivates structural benchmarking and productive-conformation assessment; COMPASS motivates linker/machinery compatibility checks. They are adapted concepts, not reproduced implementations. No new MD or static 'contact persistence' claim is required.

Exit condition: independent recalculation on a small selected set matches; atom selections and quantities have physical definitions; results are stable enough to report or their uncertainty is explicit.

## 7. Phase E — static ensemble GNN, after data freeze

Deliverables: graph_schema.json, split_manifest.tsv, model/config versions, per-fold weights, replay receipts, CPU verification and held-out_predictions.tsv.

1. Make one graph per pose with target/VHL residues and explicit ligand representation; distinguish target, VHL, ELOB/ELOC and ligand roles. Preserve ligand chemistry and pose geometry. Pool across poses/clusters into one prediction per compound-target pair. Record residue mapping, construct coverage and pose masks.
2. Compare an interface-only graph with a lysine/catalytic-feature graph under the same folds, seeds and tuning budget. Include contact/productive-fraction baselines; any statistical baseline is a declared comparator, not a replacement objective. Do not feed unavailable MD features as zero or reuse old weights under an incompatible schema.
3. Keep all related poses, structure variants and replicates from a compound-target pair together. For unseen-target claims hold out the target across ALL compounds, with family/sequence clustering sensitivity; for unseen-compound claims leave one compound out, with related-warhead/linker sensitivity. If label coverage is inadequate, report that evaluation as not estimable.
4. Fit preprocessing, feature selection, class weights, calibration and thresholds on training folds only. With standard BCEWithLogitsLoss and positives encoded 1, pos_weight weights POSITIVE errors; setting it to 3 does not upweight the negative class. If used, compute n_negative/n_positive in each training fold, or use explicitly defined class-balanced sample loss. This corrects the earlier inverted weighting suggestion.
5. Freeze multiple seeds and a small model/hyperparameter budget; use nested/grouped validation when feasible. With eight SJF8240 labels, a high-capacity GNN remains exploratory. Pretraining requires audited compatible graphs and removal of held-out target/compound overlap; unlabelled kinases do not become degradation labels.
6. Report AUROC, AUPRC with prevalence, balanced accuracy, sensitivity, specificity, confusion matrices, all held-out scores and per-seed variation. Estimate uncertainty over independent target/assay units, not poses; explain instability from two negatives. Compare models on identical held-out examples and retain all outcomes.

SIK3/RPS6KB2 and unassayed kinases receive prediction-only outputs. Separate structural compatibility from predicted degradation; do not claim a calibrated probability without suitable independent calibration data.

## 8. Immediate order of work

1. MS143 identity registration is complete (Tocris 8081); prepare its exact heads and reference mapping. It is the AZD5363-derived ATP-competitive VHL degrader, not an allosteric degrader.
2. Reconcile current raw files and receipts with stale pilot statuses; produce a factual docking inventory. Preserve stopped work and historical molecule branches.
3. Rebuild the current target/compound coverage table and assay-linked label table from existing records. These are small read-only/data-curation steps.
4. Prepare a frozen static-only reference benchmark and pilot budget, then repair the demonstrated structure-placement failures.
5. Run the bounded docking pilot after qualification, then static feature extraction and independent checks.
6. Expand docking in measured batches; develop and evaluate the static ensemble GNN only when the labelled dataset and split policy are ready.

No MD launch is part of this sequence. This planning request itself does not restart the paused pilot or submit remote jobs.

## 9. Proposed ordered scripts (planned, not implemented here)

1_inventory_inputs.py → source/identity/structure/receipt inventory.
2_curate_compounds_labels.py → six compound graphs and assay-specific labels.
3_prepare_validate_references.py → bound heads, pocket maps and VHL machinery references.
4_benchmark_docking.py → recovery, seed diversity and resource accounting.
5_plan_run_static_screen.py → resumable, receipt-bound pRosettaC batches.
6_validate_cluster_poses.py → complete structural QC and pose clustering.
7_extract_static_features.py → contacts, SASA and productive geometry.
8_build_ensemble_graphs.py → versioned static graphs and ensemble membership.
9_train_evaluate_gnn.py → grouped evaluation, training-only transforms, weights and receipts.
10_verify_report.py → independent metric/CPU replay checks and prediction-only reporting.

Reuse or wrap audited existing code where appropriate instead of copying its assumptions blindly. Verify changes with meaningful identity, chain-role, missing-atom, alignment, denominator and grouped-split regression checks. The existing 91-test report is historical software evidence, not a new folder-wide validation certificate.

## 10. Source records and boundaries

Local paths relative to /home/ayushregmibagale/10PROTAC-MD:
- config/molecules.json; molecules/verification.json
- README.md; docs/REMEDIATION_STATUS.md; docs/ML_TRANSFER_SCOPE.md
- validation_evidence/verified_checkpoint_summary.json
- SJF8240_kinome_screen_analysis/repair_v3/FINDINGS.md
- SJF8240_renanalysis2/ACTIVE_WORK.md
- SJF8240_renanalysis2/supplier_correct_pilot_v1/{STATUS.json,TODO_AND_RESUME.md,WORK_LOG.md}

Raw stopped benchmark summaries: /tmp/sjf_supplier_pilot_runs/global_5T35_MZ1_seed42/engine_raw/result_summary.txt and corresponding global_7PI4_GSK215_seed42 path. /tmp is not the permanent raw-data archive; inventory and archive relevant outputs before cleanup.

Method references:
- Bai 2022: https://pmc.ncbi.nlm.nih.gov/articles/PMC9019245/
- PROTACMap: https://pmc.ncbi.nlm.nih.gov/articles/PMC13404913/
- COMPASS: https://pmc.ncbi.nlm.nih.gov/articles/PMC13372713/
- Native VHL mechanism: https://pmc.ncbi.nlm.nih.gov/articles/PMC11468923/
- MS143 primary identity/binding mode: https://pmc.ncbi.nlm.nih.gov/articles/PMC8900464/

This review inspected selected registry, status, receipt and report records. It did not independently rerun the test suite, authenticate every structure, quantify every pose file or certify a completed campaign. The source snapshot records exactly what was used.

## 11. Resumption instruction — 15 November 2026

The user requested continuation and regular GitHub updates, prioritizing pRosettaC on Ernie as soon as practical. This authorizes progressing toward a bounded remote docking pilot and then staged release under measured resource and structural checks; it does not revive the old MD plan. Current next actions: MS143 bound-head preparation; reconcile stopped reference outputs; check remote access/runtime/storage and prepare an isolated repaired engine. Do not launch the legacy invalid inputs or certify failed recovery as successful.

## 12. Revised sampling strategy after scale review

The user challenged the cost of200global ×50local ×100conformers across6compounds ×600kinases. These are inherited defaults, not an authorized production budget. No full-scale campaign has launched. No MD will be run for these six compounds.

Published methods are not equivalent counting units:
- Bai2022:20initial global decoys; top10seed a10×1000global search; top5000clustered at3Å;200cluster representatives seed200×1000local decoys; top20000considered for full ternary generation. These are docking decoys, not200000valid full-PROTAC poses. https://pmc.ncbi.nlm.nih.gov/articles/PMC9019245/
- PROTACMap: approximately400–1400generated protein-pair conformations per system, followed by filtering/clustering and subsequent linker assessment. This is not the cost of an identical number of full pRosettaC models. https://pmc.ncbi.nlm.nih.gov/articles/PMC13404913/
- COMPASS:20000RosettaDock candidate protein arrangements per system, then filtering/clustering and separate linker screening. https://pmc.ncbi.nlm.nih.gov/articles/PMC13372713/

Proposed budget to TEST on references and a representative subset:
1. Broad first pass: N=50accepted global orientations ×10local poses ×up to20linker conformers. This is500local arrangements and at most10000conformer proposals per compound-target pair, one seed initially. At3600eligible pairs this is1.8million local arrangements and up to36million conformer proposals, excluding global rejected trials/internal embedding retries. This is100-fold fewer proposals than the inherited defaults, not a measured100-fold speedup.
2. Selective expansion: N=100 ×20local ×up to50conformers =2000local arrangements/up to100000proposals per pair. Use for sampling-unstable cases, benchmark controls, and a prespecified representative sample of apparently negative/no-hit cases. Do not deepen only promising predictions or select sampling by degradation labels.
3. Deep diagnostics only: inherited200×50×100, if needed to investigate convergence on a small subset and if measured cost permits. Not the default for3600pairs.

The first-pass budget is provisional: validate recovery, valid-pose diversity and runtime before release. Use three seeds on selected controls to assess reproducibility, not automatically every pair. Keep primary comparisons at a common budget, and report expanded analyses separately to avoid sampling-effort bias in productive fractions and GNN inputs. If the small budget misses known structures or changes rankings substantially, revise the policy before scaling rather than declaring unproductive targets negative.

Keep a proposed20–50diverse cluster representatives for downstream GNN input, testing score stability to that cap. Do not fabricate additional representatives when few valid poses exist. All raw counters/statuses remain available; sampled conformers are not independent biological observations or equilibrium populations.
