# Active campaign checkpoint — 15 November 2026

Six supplier-verified compounds: GSK215, SJF8240, SJF1521, JNJ1013, SJF_delta, MS143. No new MD. MS143 is ATP-competitive, not allosteric. MZ1 is an external method control, not a seventh screening compound.

User asks to prioritize submission before the weekly allowance runs out, while retaining publishable-quality checks. Prepare concrete inputs, quality reports and resource commands before the final submission approval. Do not launch the full kinome campaign while qualification fails. Existing remote reference pilot is complete; no additional remote job has been submitted.

## Completed checks

Ernie 50 global x10 local x20 conformer pilot: GSK215 3/8 retained poses pass the frozen coarse structural gate; top1 passes; best target RMSD after VHL alignment 1.5628 A. MZ1 0/187 pass, best 8.3451 A. Integrity receipts pass. These are reference redocking results, not degradation validation. Native bound binary structures/heads were used; native ternary arrangement did not guide global generation.

Full lossless raw archive: /mnt/c/Users/Ayush Regmi Bagale/10PROTAC_static_campaign_raw/ernie_reference_pilot_20261115/full_pilot_raw.tar.gz . SHA256 cc4ae4ec61ab7c4de4c87675cce1a74e1104ae5702956bd0aea740a468114554 . Remote originals retained.

Native fixed-head linker construction passes for both controls at seeds42/7/123; incompatible head placements reject. This is construction, not recovery or steric validation. Explicit source chain mapping confirms B/C ELOB/ELOC swap between references. MS143 strict seed4 remains the candidate; Rosetta retains phosphorylated Thr308. Modification presence/missing-phosphate/undeclared-residue checks pass. Existing ATOM-only geometry filters omit HETATM TPO, so explicit PTM-aware downstream selection is required.

Prepared cross-kinase structures are not production-qualified: many transferred heads clash, including backbone clashes. SJF_delta MAPK13 overlay is unsuitable. Its MET parent-head control must not be described as a confirmed degradation positive. Source calibrated steric reports are in pocket_obstacle_review/. Adapter is seeded random-global/PyRosetta/fixed-head linker workflow, not exact published PatchDock pRosettaC; see ENGINE_METHOD_REVIEW.md.

## Running local comparison

/tmp/six_protac_matched_budget_comparison_v1/STATUS.json ; controller unified exec session41483. Four cores,24GiB,10GiB output cap. Budget500 global x1 local x20 conformers. MZ1 complete:0/201 pass, best13.2798 A. More global sampling alone did not fix recovery. GSK215 was still running at this checkpoint. Check fresh status; do not relaunch completed tasks. STOP sentinel at that root stops the controller; namespace PIDs are not host PIDs.

When complete, use scripts07 and01 to independently verify and evaluate both tasks, archive raw records outside Git, then checkpoint compact reports. Global coverage diagnostic /tmp/six_protac_global_coverage_diagnosis_v1 has exact first50 agreement with original pilot; no native-based orientation selection permitted.

## Next work

Resolve search/refinement and pocket fit before full production. Do not keep increasing counts blindly or relax thresholds to force passing. If preparing a bounded six-compound technical/resource pilot, explicitly label its scope and unresolved biological limitations; assemble all inputs/checks before requesting approval. No MD. Keep unknown SIK3/RPS6KB2 labels unknown. Preserve original Linux working tree and its unrelated changes. Work and GitHub checkpoint from this isolated writable checkout.

16 November update: local comparison completed both tasks. GSK215 assessment pending; do not relaunch. Read RESUME_HERE.md for durable handoff.

16 November assessment completed: both local tasks integrity PASS. GSK0/5pass (best5.6865A), MZ1 0/201pass (best13.2798A). Full raw archive receipt and compact reports saved. All500paired orientations show the closest initial orientations are lost in local refinement. See LOCAL_COMPARISON_FINDINGS.md. Representation diagnostic pending. No new Ernie jobs.

Confirmed implementation defect: inherited DockingLowRes receives full-atom pose with no centroid conversion; diagnostic before/after fullatom=true and0CENproteinresidues. Separate explicit-conversion/restore smoke test passes and heads remain within0.001A in firstGSKpose. Frozen local repair benchmark now RUNNING: /tmp/six_protac_centroid_repair_pilot_v1/STATUS.json , controller unified session92185. Budget50x10x20, two references, fourcores/24GiB/10GiB/4h-per-reference. Do not relaunch; inspect status/logs. Original baseline unchanged. Next independently verify/recover/steric-review completed repair outputs. No Ernie submission.

Latest partial repair result16November: MZ1 integrityPASS,1/129coarsegatepass,best4.4000A,rank16,60%nativecontacts. All50global inputs byte-identical to originalbaseline. Passingpose graph/stereo matches, ligand-proteinmaxoverlap0.3795A. GSK215 stillRUNNING. See centroid_repair_results/PARTIAL_FINDINGS.md; not production qualification.

MS143 centroid/PTM continuation: uncorrected sidechain restoration drops phosphate atoms; pose-levelfullatomflag falselylooksadequate. A separate explicit-original-fullatomTPO replacement passes no-motionroundtrip,allatomsrestored,0.0Adisplacement,finiteenergy. Check moved-frame/full-docking compatibility before integration; current running reference candidate unchanged. Script16 and failed/pass receipts under MS143/centroid_modification_review/.

MS143 follow-up: explicit originalTPO restoration now also passes a90degree rotation+(13,-7,5)Atranslation static control,allTPOatoms match expectedmovedframe,finitefullatomscore. See explicit_TPO_restoration_moved_frame.json. Full-docking integration remains pending; running reference snapshot unchanged.

Method-selection review completed: read method_selection_review/RECOMMENDATION.md, METHOD_MATRIX.tsv and DECISION.json. Eight pinned repositories inspected at selected source files, no third-party training/replication performed. Closest scientific template=Bai2022; VHL reference=Crowe2024; reuse COMPASS static accessibility ideas/code with explicit mapping/QC; keep existing compact static GNN, borrow Pocket-PROTACs ensemble pooling/training design; DegradeMaster baseline; optionally benchmark pretrained DeepTernary for generation speed/quality. Original Bai2021 Rosetta code supports foretinib/kinase/VHL and is an alternative reference implementation. No full method adoption or remote submission yet. No MD. GNN work need not precede docking launch.

LATEST16November: centroid-repair reference benchmark FULLY COMPLETE, both integrityPASS. MZ1 1/129gatepass best4.4000A rank16; GSK215 0/3gatepass best35.1094A. Generator NOT qualified; no local benchmark still running. Complete reports and durable rawarchive receipt in centroid_repair_results_complete/. Method review recommends common-control comparison of original/standard Rosetta and pretrained DeepTernary before choosing generator. No new comparisons or Ernie jobs yet launched. Source review and recommended GNN approach documented in method_selection_review/.

LATEST USER STEERING: stop broader method-comparison work and focus solely on the six-PROTAC Ernie campaign. Read ERNIE_CAMPAIGN_PRIORITY.md. No new alternative-model research project or GNN work in immediate scope. Resolve only launch-critical qualification and assemble concrete inputs/scripts/resources/commands for final approval. No MD. No full campaign submitted.

LATEST16Novemberactualwork: anchor drift audit430/500MZ and452/500GSKoutside sampled interval. Separate soft-anchor repair smoke finishedboth,integrityPASS;5staticcontractcontrolsPASSincludingactualTPOrestoration. FULL50x10x20qualification RUNNING at /tmp/six_protac_anchor_repair_qualification_v1,controller session11223,GSKthenMZ,4cores24GiB10GiB4h/reference. Do not relaunch/mutate;STOPsentinel available. Three in-range failedGSKgeometries:20vs100embeddingretries recovers0/3both; more requested conformers maynotincreasesearchdueinitial-failureearlystop. Six technicalinputs nowstagedandhash/steric-audited in ernie_submission_v1/; private six-taskdraftbuiltat /tmp/six_protac_ernie_pilot_draft_v1,executiongateverifiedHOLD. Read-onlyErnieconnectivityPASS,96visibleCPUs~1.49TBfree,no new remotejobs. See anchor_repair_candidate_v1/FINDINGS_AND_NEXT.md and ernie_submission_v1/README.md. Only launch-critical campaign work;noMD/noGNN/no widercomparison.

USER-REQUESTED MEMORY: retain Bai2022-style full-machinery analysis after small-complex docking. Read FULL_CRL2_ANALYSIS_MEMORY.md. Current docking=target+PROTAC+VCB; planned downstream=align saved poses throughVHL onto checked fullCRL2/E2/Ub references, assess clashes/exposedlysinegeometry/productivefractions. Not yet integrated for sixdrugs. Validate VHL-specific reference/atoms/cutoffs; do not blindly useCRBN16A. No newMD or automatic redocking requirement. Ernie launch-critical work remains first.

P4WARD CHECK16November: author code/documentation and bundledCRLmodels reviewed at2d5cf1c0dc88083995b707371199b56444513fe6; see method_selection_review/P4WARD/REVIEW.md. Strong candidate for downstream static fullCRLclash/lysine-access analysis without replacing ongoing docking. Current source uses1VHLreference and fixedUbproxy;VHL60Acriterion is not an E2-Cys distance. Validate point/atom/chain mappings, SASA context, heavyatom selection and switch behavior before reuse. Not installed/integrated or end-to-endvalidated. NoMD, no pipeline switch or newErniejobs.


LATEST SAVED CHECK 16 November 2026: soft-anchor qualification COMPLETE, both integrity PASS. Frozen gate: GSK215 1/11 retained poses passes, best RMSD 3.71393 A, top1 fails/top10 contains pass; MZ1 0/392 passes, best RMSD 6.14650 A. Generator still NOT qualified across both controls. No relaunch; six-drug Ernie package remains HOLD. Compact receipts in anchor_repair_results_complete/. Consolidated papers, reusable components, corrections, source-access limitations and campaign decisions saved in PAPERS_AND_REUSE_NOTES_2026-11-16.md. No MD or new GNN work.
