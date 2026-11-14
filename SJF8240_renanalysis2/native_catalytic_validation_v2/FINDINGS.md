# Native catalytic-reference validation and follow-through — 14 November 2026

All four requested workstreams were assessed: native control/reference correction,reference-ensemble checks,the full eight-target MD comparison,and GNN-input/wider-kinome assessment. The native structural check improves,but the independently sourced primary descriptor does not separate the eight targets: AUROC0.500. No new biological classifier is released.

## 1. Correction to the earlier native-control interpretation

The earlier control superposed crystallographic5T35 BRD4 through VHL and carried over the crystal distal backbone/sidechain placement. That is not the target geometry present in the catalytic cryo-EM model. It made K456 appear23–25Å from the modeled donor endpoint. The corrected control keeps the deposited catalytic BRD4 backbone and transfers only missing sidechain heavy atoms through residue-local N/CA/C fits,with matching amino acids and RMSD≤0.25Å. Existing observed sidechain coordinates remain untouched. Crystal copiesA/E are retained as a sensitivity. No rotamers were optimized toward the donor or the eight labels.

After that correction,K456 in8RX0 is approximately19.14–19.18Å from the prior modeled Gly76 endpoint,with exposed NZ. This is approximate consistency with the published18.4Å mechanistic illustration,not an exact numerical reproduction: our missing-tail and sidechain models remain distinct. The original native-control failure must therefore be qualified as a representation/correspondence problem,not evidence that the native VHL system is unproductive.

The old hemisphere orientation rule still rejects these frozen K456 rotamers,so it is not reinstated as a mandatory reaction criterion. A static lysine rotamer is not an experimentally established transition-state attack geometry.

## 2. Independently sourced catalytic-reference test

We downloaded8R5H from the Li et al.(2024) VHL–MZ1–BRD4–UBE2R2 study and its SY8 chemical-component definition. It contains a chemically trapped three-way junction. PDB connectivity explicitly links E2Cys93SG to SY8C3 and UbGly75C to SY8N. CCD bonds give the path N–C4–C3. Thus C3 is a deposited,protein-tethered reaction-center surrogate;it is not native ubiquitinGly76 carbonylC. SY8's oxygen is attached elsewhere in the artificial linker. We do not silently relabel this as a resolved physiological thioester.

The C93SG–SY8C3 distance is approximately1.705Å. With the native observed scaffold,target and missing-sidechain procedure,the known K445 site is16.46Å from that surrogate,exposed,and within the coarse20Å geometry criterion. A500Å-translated target fails reach. The trapped K368C site is retained asCYS and never counted as a wild-type lysine.

Three declared references were constructed before eight-target rescoring: deposited8R5H/SY8,8RX0 with the prior modeled Gly76 fit75,and8RX0 with SY8 transferred through the common E2 active-loop backbone85–100. The active-loop fit is0.558Å. No endpoints were fitted to target degradation outcomes. Two anchors,VHL and ELOC,are compared. Missing atoms/components remain explicit:8R5H has no resolvedNEDD8 and8RX0 has multiple backbone-only catalytic components. These are not complete atomic reaction simulations.

## 3. Full eight-target MD result

All141 retained cases and201 saved frames each were processed:28,341 frames. The primary reference/anchor was frozen as8R5H/SY8 with VHL anchoring. The primary rule is: at least one target lysineNZ within20Å of C3,accessible-sphere fraction≥0.10,and maximum target/new-scaffold heavy-atom overlap≤0.8Å. No orientation or bridging condition is imposed on that primary productive fraction. Joint-bridge and fit3Å variants are separate sensitivities. Target-induced clashes remain in the denominator;frames and lysines are not independent labelled observations.

| Target | User label | Primary productive fraction,% | Joint bridged/productive fraction,% |
|---|---|---:|---:|
| ABL1 | Non-degrader | 11.384 | 3.922 |
| AXL | Non-degrader | 0.967 | 0.912 |
| DDR2 | Degrader | 0.221 | 0.138 |
| EPHA2 | Degrader | 5.777 | 5.500 |
| EPHB2 | Degrader | 17.334 | 11.836 |
| MAPK14_full | Degrader | 0.055 | 0.055 |
| MET | Degrader | 12.220 | 11.194 |
| RIPK2_full | Degrader | 1.112 | 0.673 |

Primary target-ranking AUROC=0.500;exact two-sided permutation p=1.000 across all28 allocations of two negative labels. All eight targets have nonzero scores. The unvalidated score>0 classifier therefore has sensitivity1,specificity0,balanced accuracy0.5. ABL1 is geometrically compatible despite being a user-confirmed non-degrader. No threshold was tuned to remove it.

| Reference | Anchor | AUROC of the same productive-fraction rule |
|---|---|---:|
| 8R5H_observed_SY8 | ELOC | 0.458 |
| 8R5H_observed_SY8 | VHL | 0.500 |
| 8RX0_E2_aligned_SY8 | ELOC | 0.750 |
| 8RX0_E2_aligned_SY8 | VHL | 0.625 |
| 8RX0_G76_fit75 | ELOC | 0.833 |
| 8RX0_G76_fit75 | VHL | 0.625 |

The earlier0.750 result remains a different exploratory analysis with different endpoint/scaffold choices and additional rules;it was not validated by the new controls. Some alternative reference/anchor combinations score higher,but selecting those after examining eight labels would overstate performance. The stable conclusion is reference dependence and inadequate current evidence for a useful degradation classifier. This does not prove that full-machinery physics cannot help,or that a GNN cannot use finer patterns.

The exposure is the NZ accessible-sphere fraction from180 Fibonacci points and1.4Å probe,not an experimental SASA measurement. Overlaps use declared heavy-atom radii and omit missing sidechain atoms;they are not MolProbity scores. Distances to the SY8 surrogate are not catalytic probabilities.

## 4. GNN input assessment

Prepared16 matched four-window graphs:eight with five provisional native-reference lysine-node descriptors and eight zero-padded controls. Existing29 node features and every edge are byte-identical between controls and augmented inputs;both have34 features. Added columns are mean NZ–SY8C3 distance/100Å,reachable20Å fraction,exposure fraction,exposed/reachable/scaffold-compatible fraction,and joint bridge/productive fraction. These are explicitly engineered descriptors,not new full-scaffold nodes or observed ubiquitination labels. Each new descriptor uses equal pose weighting over the full201-frame data,while the original base features remain unchanged.

The existing TemporalEncoder was software-tested with only the first GAT input width changed29→34,using identical34-input architecture forboth variants. All16 graphs pass finite untrained CPU forward checks;there were zero optimizer steps. An initial unguardedPyG import hit an incompatible optional native extension;the existing process-local native-PyG guard resolves it without modifying installed packages. No trained predictions or replacement weights were generated.

This completes input/runtime readiness assessment,not GNN validation. A new training campaign is not automatically justified by a descriptor whose independent primary ranking is chance-level and whose historical chemistry remains unresolved. If a controlled exploratory GNN ablation is subsequently pursued,it must hold out entire proteins,retain the same architecture/splits/seeds/epochs forboth variants,and report all folds and CPU replay. Eight labelled targets cannot establish generalization regardless of graph/frame count.

## 5. Wider kinome assessment

Annotated all723 original QC-shortlist poses across196 accessions with the deposited native-reference surrogate. 89 pass the raw exposed/reachable/no-added-scaffold-clash test. They are hypotheses about geometry only,not predicted binding or degradation. No unknown target was assigned a negative experimental label.

There are0 scientifically eligible new-MD candidates in this unrepaired archive. All723 inherit the historical molecule-connectivity mismatch and shared internal linker-overlap problem;original head-pocket placement validation is unresolved. New compatibility annotation cannot cure those defects. The existing ten-pose local repairs were a separate bounded historical-model diagnostic and are not substituted into the723-pose population. The top-ten table is labelled DIAGNOSTIC_ONLY and must not become a production shortlist.

## Verification and saved evidence

- Native PDB/CCD connectivity,local sidechain-transfer receipts,and every native lysine distance/exposure/clash are saved.
- All141 MD source trajectories/topologies/reference PDBs were rehashed against the earlier source-bound audit. Contact and frame-grid replay covers all28,341 frames.
- Independent MDAnalysis-unwrapping/Bio.SVDSuperimposer/dense-distance replay checks six raw frames from ABL1/DDR2 for distances,overlaps and exposure. Four target AUROCs are independently verified by pair counting.
-16 matched GNN graphs retain original base tensors and edges;finite untrained forward checks pass.
-723 archived kinome source PDBs were hash-checked and sequence-role checked by the existing parser/context contracts.
- Original structures,trajectories,labels,models and previous outputs remain unchanged. No new docking,MD production or GNN training was run.

## Decision

The native-control correspondence issue is corrected and qualitative20Å reference access is recovered. The full native-reference rescoring does not separate the eight labels. The GNN data path is prepared but untrained,and the wider archive remains unsuitable for supplier-SJF8240 selection. The justified next scientific work is an independent assay/chemical-provenance review and a validated candidate ensemble,not further tuning of geometric thresholds on these same eight labels. The chemistry-correct and historical-model branches remain separate. No general claim about binder-versus-degrader prediction is supported.

## Sources

- Crowe et al.(2024),VHL catalytic architecture: https://pmc.ncbi.nlm.nih.gov/articles/PMC11468923/
- Li et al.(2024),VHL/UBE2R2 substrate priming: https://pmc.ncbi.nlm.nih.gov/articles/PMC10997478/
- Bai et al.(2022),productive-fraction concept: https://pmc.ncbi.nlm.nih.gov/articles/PMC9019245/
- Deposited8R5H andCCD: https://www.rcsb.org/structure/8R5H ; https://files.rcsb.org/ligands/view/SY8.cif

## Additional chemical-tail reconstruction: saved after the primary comparison

A chemically explicit missingGly76 tail was reconstructed using only the ubiquitinGly75 backbone andE2Cys93 sulfur placement. No substrate lysine or degradation label enters these fits. Early fixed-atom trials were retained:they could satisfy some bonds while violating neighboring peptide planarity or creating clashes,so they were not released. The peptide carbonyl planes,trans-amide geometry,declared ideal bond lengths/angles and bonded1–2/1–3 exclusions were checked explicitly.

The final bounded attempt permitted only the observedGly75 carbonyl oxygen plus missingGly76 atoms to move,with all other protein backbone atoms fixed. It achieved the declared geometry/clash gates for8RX0 after movingGly75O by0.944Å. Independent dense-distance replay over all deposited E2/Ub atoms confirms maximum checked nonbonded overlap below0.8Å and C–S distance1.800Å. Cys93SG remains transferred/modelled because it is absent from8RX0;unresolved sidechains remain absent. This is a geometrically consistent local model,not a physiological reaction simulation or an energy-minimized complete catalytic complex. The corresponding8R5H constrained reconstruction failed the unchanged gates and is not promoted.

The new8RX0 modeled endpoint places native-backbone K456 at approximately18.7Å,close to the paper's18.4Å illustration. This is an independent geometric consistency check,not exact reproduction and not a distance fitted to K456.

The corrected endpoint was then assessed on all141 poses/28,341 frames as a separately frozen secondary sensitivity,not a replacement for the initial8R5H primary. Both raw-frame geometry and AUROC arithmetic have an independent replay.

Secondary productive-fraction AUROC=0.625,exact target permutation p=0.678571. The original primary remains0.500. Neither supports a reliable degradation classifier. All alternative thresholds and orientations remain in chemical_tail_sensitivity/results/;no best-case variant is selected as validation.

The matched GNN inputs still encode the declared primary SY8 surrogate;they were not silently replaced by this secondary endpoint. The full kinome annotation is likewise the declared primary diagnostic. Keeping these separate prevents a succession of physical reference changes from becoming undocumented outcome-driven feature selection.

## Completion and remaining scientific limits

Completed:native-backbone correspondence correction;three-reference preparation andcontrols;two full eight-target MD comparisons including the locally repaired modeled tail;independent coordinate/metric checks;16 matched GNN input/runtime checks;723-pose kinome diagnostic;all scripts,receipts and findings. No new training or new simulations were performed.

Still unresolved:full physiological catalytic-reference validation,the archived-versus-supplier molecular identity,correct pocket placement throughout the kinome,and independent biological validation with more than two non-degrader targets. These are explicit limitations,not PASS statuses. The saved provisional GNN data do not constitute a trained or validated GNN result.
