# Six pilot controls and paper-backed targets

The original six-control run used the exact compounds listed in its frozen manifest. Scientific interpretation distinguishes the compound's reported cellular targets from the target and conformation chosen for a docking control.

| Compound | Primary evidence | Project control and current limitation |
|---|---|---|
| GSK215 | FAK/PTK2 is the target. 7PI4 is a deposited FAK–GSK215–VCB structure. [Law et al.](https://onlinelibrary.wiley.com/doi/10.1002/anie.202109237) · [PDB 7PI4](https://www.rcsb.org/structure/7PI4) | PTK2 is the native structural benchmark. Independent review verifies 37/37 retained ligand identities. The matching 50 × 2 benchmark yielded three coarse native-recovery passes outside the top 10 by raw score. |
| SJF8240 | Burslem's compound 7 induced MET degradation. Smith et al. also reported submicromolar degradation of p38α and p38δ in their assay. [Burslem et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC5831399/) · [Smith et al.](https://www.nature.com/articles/s41467-018-08027-7) | MET is a biologically relevant target. The MET pilot produced zero feasible linkers. A synthetic same-molecule head geometry closed successfully, while three archived MET refinements with the smallest anchor-distance shifts remained infeasible with 100 internal embedding iterations per outer attempt. The existing p38δ placement has 1.43 Å of maximum backbone overlap. |
| SJF1521 | Burslem's lapatinib-based compound 5 favors EGFR degradation in OVCAR8. Burslem's gefitinib-based compound 3 is a distinct PROTAC. [Burslem et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC5831399/) · [SJF1521 follow-up study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11688428/) | EGFR matches the reported target. The head is modeled and its original binary maximum overlap is 0.25 Å. The pilot anchor indices 1 and 1 map to full-molecule atoms N24 and O37 in the stored SMILES. Independent review verifies 55/55 retained ligand identities. |
| JNJ1013 | Degrader-3 targets IRAK1 and achieved a DC₅₀ of 3 nM in HBL-1. [JNJ1013 discovery paper](https://pubmed.ncbi.nlm.nih.gov/34279092/) | IRAK1 matches the reported target. The complete run with repaired SDF serialization scored three poses. The head audit initially rejected them because two files had the same hash. The repaired audit resolves the run-declared head filenames and verifies their hashes. All three scored poses pass on re-audit. |
| SJFδ | p38δ/MAPK13 is the selective degradation target. Smith et al. measured about 46 nM DC₅₀ for p38δ and no detectable degradation of p38α/MAPK14, β or γ under their assay conditions. [Smith et al.](https://www.nature.com/articles/s41467-018-08027-7) | MET served as a technical control. All 50 retained MET poses have the correct ligand identity, but their median maximum ligand–protein overlap is 1.52 Å. The modeled MAPK13 input has improved after a rigid refit and protein repack. Its maximum overlap fell from 1.61 to 0.56 Å, with 1.13 Å RMS head displacement from the starting overlay. It remains a modeled binary with no ternary validation. |
| MS143 | Compound 20 is a VHL-recruiting AKT degrader derived from AZD5363. The paper establishes cellular AKT degradation and does not establish AKT1-only selectivity. [MS143 discovery paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8900464/) | AKT1 is a reasonable structural control. The modeled head uses 4GV1 as a capivasertib-bound template and retains phospho-Thr308. The repaired archived worker produced one mechanical-QC passing pose and retained the declared phosphate atoms. |

## Independent steric readouts

| Pilot control | Identity pass | Maximum ligand–protein overlap, min / median / max (Å) | Maximum VBC–target overlap, min / median / max (Å) |
|---|---:|---:|---:|
| GSK215 / PTK2 | 37/37 | 0.693 / 0.693 / 1.185 | −0.225 / 0.259 / 1.863 |
| SJF1521 / EGFR | 55/55 | 0.715 / 0.715 / 1.515 | −0.404 / 0.393 / 1.003 |
| SJFδ / MET | 50/50 | 0.849 / 1.516 / 1.882 | 0.016 / 0.432 / 0.759 |

These are element-radius geometric readouts, alongside the pilot's minimum-distance screen. The reference-like posture, cellular target evidence and score rank are separate observations. Bulk coverage and submission were not exercised.

## Additional checks

All 12 declared pilot anchor assignments have a unique full-ligand attachment-atom identity across their chirality-aware matches. The detailed one-based mappings are in `ANCHOR_MAPPING_REVIEW.json`.

The 37 retained GSK215 pilot PDBs are byte-identical to the earlier 50 × 2 reference run. Their scores match exactly by pose ID. Only the score row ordering differs.

The provenance repair reproduced the original bound-head QC on representative GSK215, SJF1521 and SJFδ retained poses. Altering the bytes of a declared head caused a hash rejection as intended.

The modeled SJFδ/MAPK13 run used four global poses and one refinement each, with five conformers per feasible geometry. It completed with three scored and two retained ternaries in 146 seconds. Both retained ligands passed independent identity review. Their maximum ligand–protein vdW overlaps were 1.22 Å and 1.80 Å. The prepared binary's improvement therefore did not carry through to contact geometry in these two recovered ternaries.

## Terminal MS143 validation

The full repaired MS143 run finished all 100 refinements, scored 27 conformer-built ternaries and exited zero. The earlier hash-uniqueness QC classified all 27 as rejected. Re-audit with the declared-head path repair accepts all 27 under the existing mechanical rules. All 27 retain the complete declared phospho-Thr308 heavy-atom set and pass the 0.05 Å bound-head core displacement check. The underlying scored PDBs and scores were preserved. These are technical checks on a modeled AKT1 control.
