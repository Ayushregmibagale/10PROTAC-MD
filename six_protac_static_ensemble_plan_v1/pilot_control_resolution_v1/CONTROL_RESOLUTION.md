# Remaining control issues

## GSK215 / PTK2 ranking

Compared archived REF2015 totals with two fixed-coordinate score differences. The protein-interface diagnostic subtracts isolated VBC and PTK2 scores from the combined proteins. The ligand–protein diagnostic subtracts the combined proteins and the isolated PROTAC from the full complex. These differences are model scores. They do not estimate binding free energy. Rescoring the reloaded PDBs differed from recorded totals by up to 1.85 score units in the newer ensemble.

| Ensemble | Raw score passing ranks | Protein-interface passing ranks | Ligand–protein passing ranks |
|---|---|---|---|
| Repaired 50 × 2 pilot, 37 retained | 15, 16, 25 | 9, 10, 11 | 1, 4, 6 |
| Earlier 50 × 10 pilot, 8 retained | 1, 4, 7 | 2, 4, 8 | 1, 2, 3 |

The frozen geometry labels were joined after scoring. The same native structure underlies both ensembles. The candidate ligand–protein score is exploratory and the production score has not been changed.

## SJFδ / MAPK13 sterics

The original saved conformers all retained substantial overlaps. A protein-aware UFF minimization held both complete heads fixed and retained ligand identity and bond lengths. The best worst overlap improved from 1.22 to 1.07 Å. Its coordinates are diagnostic. In the original two retained poses the worst overlap involved a linker atom and VHL Ile109 or His110.

A wider 20-global × 2-local test of the repaired binary on Ernie generated 40 refined poses, scored 36, and retained 24. All 24 retained ligands passed the independent identity check. The best worst overlap was 0.85 Å. The limiting contact in two of the better poses was a fixed VHL-head carbon against VHL Tyr112 OH. This documents a reduction in clashes. It leaves the modeled ternary control unqualified for production.

## SJF8240 / MET linker feasibility

A ligand-driven test rigidly aligned the MET and VHL binaries to each archived free SJF8240 conformer. All 189 identity-verified conformers failed to preserve both bound head shapes within the original 0.35 Å RMSD threshold. The best maximum-of-two RMSD was 1.50 Å.

A second diagnostic sampled 40 full ligands with intra-head distance bounds imposed independently on the two heads. Seven proposals satisfied both head shape checks. All seven failed the existing gross-clash screen, with minimum VBC–MET heavy-atom distances from 0.12 to 0.35 Å. These results describe the tested proposals. Additional unconstrained conformer sampling has no demonstrated route to a viable MET geometry.

Bulk preparation and release are unchanged.
