# Completed soft-anchor repair diagnosis — 16 November 2026

Both completed runs pass independent source/input/output integrity. Full raw archive was read back (15,218 members); receipt records its durable path, size and SHA256.

## Where recovery is lost

All 50 global and 500 local orientations per reference were evaluated retrospectively. MZ1 starts with no global orientation within 5 A (best 13.0828 A); refinement improves the best to 6.1465 A, still outside the frozen recovery gate. Linker construction succeeds in 480/500 local arrangements, so its dominant remaining failure in this run is protein arrangement recovery, not inability to construct the linker. GSK215 starts at best 8.1067 A; two local arrangements reach <=5 A, and one of 11 retained final poses passes the full RMSD/contact gate. Only 11/500 local arrangements yield feasible linkers, leaving a distinct short-linker feasibility bottleneck. No native geometry was used to generate or select docking starts.

The restraint reduces extreme anchor drift but is soft: 235/500 MZ1 and 364/500 GSK215 local outputs remain outside the sampled heuristic interval. Distance within that interval does not establish compatible exit vectors. These results do not prove either convergence or a fundamental engine limitation.

## Independent molecular readouts

All 403 retained poses match the reference ligand graph/stereochemistry under the existing RDKit review. MZ1 has 107 poses with ligand/protein maximum overlap more than 0.5 A above its native baseline, and 7 with excess protein/protein overlap. GSK215 has 0 and 2 respectively. These are continuous diagnostic summaries, not validated physical acceptance thresholds; more linker solutions are not automatically more useful structures.

## Next bounded decision

Do not change the model and sampling simultaneously again. The repaired source can now support a declared larger-global-coverage diagnostic with the same local protocol and independent gate, retaining all generated starts without native-based selection. This would test a genuine increase in total sampling, unlike the earlier 50x10 versus 500x1 comparison. Before launching, freeze the finite budget/resource caps and report its exploratory scope; do not treat any success as automatic six-drug qualification. In parallel, remaining MS143 PTM-aware final-QC and per-input geometry issues still require resolution. No full Ernie production release is authorized or ready.
