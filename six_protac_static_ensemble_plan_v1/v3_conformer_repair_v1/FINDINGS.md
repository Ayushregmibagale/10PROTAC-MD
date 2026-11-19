# Six-control v3 follow-up

The failed archived JNJ1013 and MS143 linker SDFs included 18–20 records per failing worker. RDKit rejected all records tested with aromatic-bond kekulization errors. The linker writer and Rosetta parameter input writer had each explicitly disabled kekulization. The separate repair runtime uses the normal RDKit writer bond-order serialization at both boundaries.

The archived JNJ1013/IRAK1 31_1 and MS143/AKT1 8_1 workers were replayed on Ernie. Both used their original docked PDB, head templates, parameter provenance, deterministic worker key and conformer budget. Both wrote 20 conformers, completed scoring, and produced one combined pose passing the existing mechanical QC and ligand checks. MS143's final PDB retained the full declared TPO308 heavy-atom set. These are two repaired worker replays. Their results do not retroactively qualify the original whole-compound attempts.

SJF1521/EGFR uses one-based anchor 1 for each head. A chirality-aware full-ligand graph check identifies the VHL nitrogen as one-based supplier-SMILES atom 24, attached to linker atom 25. The POI oxygen is atom 37, attached to linker atom 36. All observed stereochemical graph mappings give those same atom identities. This establishes the attachment identities in the frozen files. The source head is still modeled and the preparation is flagged production_ready false.

All 100 SJF8240/MET workers completed linker sampling. All 2,000 recorded embedding attempts failed. Median absolute anchor-distance shift during refinement was 0.26 Å. This specific MET technical input and budget did not produce a full ligand. The rejection ledger supports a finite linker-feasibility outcome.

The held bulk draft remains separate. It covers 433 distinct compound–kinase pairs. The current wider inventory has candidates for 490 of 2,610 pairs, including unreviewed candidates. Releasing a subset as the first tranche would change the requested all-kinase scope and requires an explicit release decision.
