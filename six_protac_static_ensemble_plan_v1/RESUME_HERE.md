# Resume the six-PROTAC static campaign

Saved 16 November 2026. Read this file, ACTIVE_WORK.md, STATUS.json, RESEARCH_CONTRACT.md and PLAN.md before continuing.

## Objective and instructions

Prepare publishable-quality static pRosettaC ensembles for six Tocris/Bio-Techne compounds across the kinase panel, followed by static geometry features and an ensemble GNN. NO NEW MD. Compounds: GSK215, SJF8240, SJF1521, JNJ1013, SJF_delta and MS143. MS143 is ATP-competitive AKT-targeting, not allosteric. MZ1 is only an external reference benchmark, not part of the six-drug screen.

Work autonomously, conserve allowance, checkpoint scripts/results on GitHub, and ask the user only when a concrete Ernie submission package is ready. Do not sacrifice qualification to meet a quota deadline. No full campaign has been submitted. Prior reference pilot authorization does not authorize a new campaign submission.

## Locations

Writable checkout: /mnt/c/Users/Ayush Regmi Bagale/10PROTAC_static_campaign_work
Original project, preserve unrelated changes: /home/ayushregmibagale/10PROTAC-MD
Study folder: six_protac_static_ensemble_plan_v1
GitHub: https://github.com/Ayushregmibagale/10PROTAC-MD ; latest prior checkpoint912b731 onmaster.

## Completed / blockers

Ernie reference pilot finished. At50global x10local x20conformers, seed42, GSK215:3/8retained poses pass frozen coarse gate, top1passes, best1.5628A target CA RMSD after VHL fit. MZ1:0/187pass, best8.3451A. Gate: RMSD<=5A and native residue-contact fraction>=0.5 plus inherited chemical/gross-clash QC. This is redocking, not biological validation. Bound binary/head coordinates derive from crystal references; native ternary arrangement did not guide generation. Integrity checks passed.

Native fixed-head linker construction works for both controls at3seeds, incompatible controls reject. Many kinase-head transfers have serious steric/backbone clashes. SJF_delta MAPK13 overlay is unsuitable; MET can only be a clearly labeled parent-head technical control. MS143 strict seed4 is a modeled candidate, Rosetta typing/load succeeds and phosphorylated Thr308 is retained. Explicit PTM-aware downstream geometry selection is needed because ATOM-only selection omits HETATM TPO. See reports and ENGINE_METHOD_REVIEW.md: current engine is adapted seeded-random-global/PyRosetta/fixed-head linker sampling, not exact published PatchDock pRosettaC. Do not invent qualification or loosen gates to force success.

## Immediate next action

The local matched-nominal-budget comparison has NOW FINISHED, as checked16November. Both tasks statusCOMPLETE, overallREFERENCE_PILOT_COMPLETE_NEEDS_RECOVERY_ASSESSMENT. Root /tmp/six_protac_matched_budget_comparison_v1 . Budget500global x1local x20conformers. MZ1 already assessed:0/201pass, best13.2798A. GSK215 completed but has NOT yet been independently assessed in saved reports. Run scripts07 and01 for the complete pilot using their --help arguments; verify all receipts, assess GSK215, preserve full raw records outside Git and update compact reports. Do not relaunch these completed tasks. Existing local_budget_comparison_results/ is partial MZ1-only and must not be mistaken for full assessment.

Next resolve search/refinement recovery and pocket fit, then assemble actual six-compound technical/resource pilot inputs, manifest, QC, runtime/storage estimates and commands for final Ernie approval. Full3600pair production is not ready. Failed docking is not a biological negative; poses are not independent biological samples. Keep SIK3/RPS6KB2 labels unknown. Never transfer SJF8240 experimental labels to other compounds.

## Raw preservation

Ernie originals /home/s1638308/ernie_static_reference_pilot_v1 retained. Local full archive /mnt/c/Users/Ayush Regmi Bagale/10PROTAC_static_campaign_raw/ernie_reference_pilot_20261115/full_pilot_raw.tar.gz ;409094335bytes; SHA256cc4ae4ec61ab7c4de4c87675cce1a74e1104ae5702956bd0aea740a468114554 . Extracted review /tmp/ernie_static_reference_review_20261115/extracted . Treat /tmp as temporary, preserve needed local comparison raw output durably.

Latest targeted tests:11passed. No new Ernie campaign running or submitted by this work. Read actual current files/status before reporting new progress.

## Latest continuation16November — takes precedence over earlier pending-status text

The previous matched-budget comparison is fully assessed and archived: both integrityPASS; GSK0/5pass,best5.6865A;MZ1 0/201pass,best13.2798A. Reports in local_budget_comparison_results/. Raw archive receipt there points to durable542145922bytearchive with SHA256a61310bf80a2da75e203793c99c9fdb13d588e4f0f9acd3d73a05272ffc65a56. Paired500orientation traces implicate local refinement. Representation diagnostic confirms missing centroid conversion in inherited low-res adapter. Explicit transition/restore smoke passes; see REFINEMENT_REPAIR_NOTES.md.

A NEW SEPARATE LOCAL bounded repair benchmark is RUNNING at /tmp/six_protac_centroid_repair_pilot_v1 . Controller session92185; read STATUS.json before any action. FirsttaskMZ1, thenGSK;50global x10local x20conformers,seed42. Candidate explicit centroid ligand typing/conversion/restoration, otherwise existing fullatom and linker pipeline; source/plan/manifests saved in refinement_candidate_v1/.4cores,24GiB,10GiB,4h/reference. STOP sentinel stops it. No new Ernie job. Next assess integrity/recovery/head engagement/independent sterics and preserve raw results. Passing execution is not structural qualification.

Latest partial repair result16November: MZ1 integrityPASS,1/129coarsegatepass,best4.4000A,rank16,60%nativecontacts. All50global inputs byte-identical to originalbaseline. Passingpose graph/stereo matches, ligand-proteinmaxoverlap0.3795A. GSK215 stillRUNNING. See centroid_repair_results/PARTIAL_FINDINGS.md; not production qualification.

MS143 centroid/PTM continuation: uncorrected sidechain restoration drops phosphate atoms; pose-levelfullatomflag falselylooksadequate. A separate explicit-original-fullatomTPO replacement passes no-motionroundtrip,allatomsrestored,0.0Adisplacement,finiteenergy. Check moved-frame/full-docking compatibility before integration; current running reference candidate unchanged. Script16 and failed/pass receipts under MS143/centroid_modification_review/.

MS143 follow-up: explicit originalTPO restoration now also passes a90degree rotation+(13,-7,5)Atranslation static control,allTPOatoms match expectedmovedframe,finitefullatomscore. See explicit_TPO_restoration_moved_frame.json. Full-docking integration remains pending; running reference snapshot unchanged.

Method-selection review completed: read method_selection_review/RECOMMENDATION.md, METHOD_MATRIX.tsv and DECISION.json. Eight pinned repositories inspected at selected source files, no third-party training/replication performed. Closest scientific template=Bai2022; VHL reference=Crowe2024; reuse COMPASS static accessibility ideas/code with explicit mapping/QC; keep existing compact static GNN, borrow Pocket-PROTACs ensemble pooling/training design; DegradeMaster baseline; optionally benchmark pretrained DeepTernary for generation speed/quality. Original Bai2021 Rosetta code supports foretinib/kinase/VHL and is an alternative reference implementation. No full method adoption or remote submission yet. No MD. GNN work need not precede docking launch.

LATEST16November: centroid-repair reference benchmark FULLY COMPLETE, both integrityPASS. MZ1 1/129gatepass best4.4000A rank16; GSK215 0/3gatepass best35.1094A. Generator NOT qualified; no local benchmark still running. Complete reports and durable rawarchive receipt in centroid_repair_results_complete/. Method review recommends common-control comparison of original/standard Rosetta and pretrained DeepTernary before choosing generator. No new comparisons or Ernie jobs yet launched. Source review and recommended GNN approach documented in method_selection_review/.

LATEST USER STEERING: stop broader method-comparison work and focus solely on the six-PROTAC Ernie campaign. Read ERNIE_CAMPAIGN_PRIORITY.md. No new alternative-model research project or GNN work in immediate scope. Resolve only launch-critical qualification and assemble concrete inputs/scripts/resources/commands for final approval. No MD. No full campaign submitted.

LATEST16Novemberactualwork: anchor drift audit430/500MZ and452/500GSKoutside sampled interval. Separate soft-anchor repair smoke finishedboth,integrityPASS;5staticcontractcontrolsPASSincludingactualTPOrestoration. FULL50x10x20qualification RUNNING at /tmp/six_protac_anchor_repair_qualification_v1,controller session11223,GSKthenMZ,4cores24GiB10GiB4h/reference. Do not relaunch/mutate;STOPsentinel available. Three in-range failedGSKgeometries:20vs100embeddingretries recovers0/3both; more requested conformers maynotincreasesearchdueinitial-failureearlystop. Six technicalinputs nowstagedandhash/steric-audited in ernie_submission_v1/; private six-taskdraftbuiltat /tmp/six_protac_ernie_pilot_draft_v1,executiongateverifiedHOLD. Read-onlyErnieconnectivityPASS,96visibleCPUs~1.49TBfree,no new remotejobs. See anchor_repair_candidate_v1/FINDINGS_AND_NEXT.md and ernie_submission_v1/README.md. Only launch-critical campaign work;noMD/noGNN/no widercomparison.
