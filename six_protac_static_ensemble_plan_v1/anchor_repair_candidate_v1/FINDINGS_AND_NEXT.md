# Launch-critical docking repair — 16 November 2026

A read-only audit of all500local poses per reference shows430/500MZ1 and452/500GSK215anchor distances outside the range sampled from each molecule's free conformers. GSKmedian10.1166A, sampled1.345–5.345A. MZmedian17.5624A, sampled7.61–11.61A. These ranges are heuristic, not proofs of physical impossibility. Local refinement does not retain the global linker's distance prior in the previous protocol.

The new isolated candidate adds a generation-only flat-bottom harmonic restraint from this existing interval, sigma1A, weight1 in low/high docking. It recreates anchor AtomIDs after residue-type conversion. Final independent raw scoring and the frozen recovery gate are unchanged. It also checks restored residue atom inventories, restoring supported phosphorylated Thr explicitly and rejecting unknown mismatches. Original/failed source and results remain untouched.

First3global x3local execution smoke completed both references with independent integrityPASS. All18local restraint receipts are present. MZ1generated9/9linkers;GSKgenerated0/9. This is execution evidence only; the firstthree global orientations are not expected to establish reference recovery. Five static contract controls pass: flat interval has zero penalty, outside interval positive penalty, missing anchor rejected, actual candidate restores phosphoThr, and unexpected protein residue substitution rejected.

Full common-budget qualification is RUNNING at /tmp/six_protac_anchor_repair_qualification_v1/STATUS.json ; unifiedcontroller session11223. GSK215 first, thenMZ1;50global x10local x20conformers,seed42;4cores24GiB RSS; output cap remains10GiB in this local qualification configuration,4h/reference. Read the frozen PILOT_CONFIG.json for authoritative caps. STOP sentinel stops it. Do not relaunch or mutate the running snapshot. No native-based selection/restraints or degradation labels are used.

## Is this a budget or engine limitation?

Current evidence favors protocol/geometry failures over a demonstrated inability to construct these molecules. Native-head linker construction succeeds forbothreferences atthree seeds. Earlier50x10 versus500x1redistributed the same500local attempts and did not improve recovery; this does not test larger total budgets. The known centroid conversion defect is separate from linker-complexity limitations.

A new bounded diagnostic held3failedGSKgeometries fixed, selected by sorted index among failures inside the original sampled distance interval. Increasing outer embedding retries20to100recovered0/3in both conditions. This small selected test is not a convergence study; it suggests that repeatedly embedding these fixed arrangements is not the immediate remedy. In-range anchor distance does not ensure compatible exit vectors or full-head geometry.

There is also an early-stop rule: when the first fixed-head feasibility search fails, remaining requested conformers are skipped. Increasing the requested conformer count alone may therefore leave actual search unchanged. Search attempts, internal iterations, generated conformers and accepted geometries must be reported separately. Do not call failed finite searches biological negatives or a fundamental limitation of all pRosettaC.

Next assess the restrained full qualification, including chemistry, binding-head preservation, structural recovery and independent sterics. Only then decide whether additional geometry/orientation constraints or deeper sampling are justified. No broader literature or GNN work is in current scope.
