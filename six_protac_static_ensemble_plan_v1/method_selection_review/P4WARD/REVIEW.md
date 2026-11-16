# P4ward: focused check for our full-CRL post-docking layer

Checked16November2026. DOI:10.1021/acs.jcim.5c00614. Sources: [publisher record](https://pubs.acs.org/doi/10.1021/acs.jcim.5c00614), [official documentation](https://skteamlab.github.io/P4ward/), [configuration reference](https://skteamlab.github.io/P4ward/config_reference.html), [author code](https://github.com/SKTeamLab/P4ward). Selected source was inspected at commit2d5cf1c0dc88083995b707371199b56444513fe6. Full publisher text was blocked; accessible abstract/supporting-information excerpt, documentation, actual code and bundled coordinates were reviewed. No installation or end-to-end run was performed.

## What it does

P4ward is an automated ternary-modeling workflow built around MEGADOCK, RDKit and RxDock, with free/open-source structural tools. Its separate CRL filter aligns docked target/E3 arrangements into bundled larger ubiquitination assemblies, checks machinery clashes and tests surface lysine access to a ubiquitin reference point. Its reported76.5%hit rate concerns structural modeling, not classification of cellular degraders versus confirmed binding nondegraders.

The accessible supporting-information excerpt describes a VHL–ElonginB/C–CUL2–NEDD8–RBX1–E2–Ub model assembled using structural templates including6TTU,5N4W and1LQB. Thus it uses a constructed structural reference, rather than docking the entire CRL independently for each PROTAC. A transplanted arm is a declared modeling assumption requiring validation, not automatically an invalid reference merely because6TTUcontainsCUL1.

## Verified implementation

The principal files are protein_filter.py, structure_tools.py, structures/model_info.py and config/default.ini. SOURCE_REVIEW.json records hashes and default settings.

1. The source lists one VHL model and nine CRBN models. The VHL model is not a sampled VHL catalytic-reference ensemble.
2. Accessible-lysine evaluation uses target LYS NZ coordinates, residue SASA from BioPython ShrakeRupley and a straight-line occlusion test. The SASA calculation uses the target alone; it does not directly calculate burial caused by binding the E3/PROTAC.
3. A fixed stored coordinate supplies the ubiquitin reference point. Code/documentation must not be relabeled an E2-Cys distance. Inspecting the bundled VHL PDB, the nearest atom to this point is GLY75 C in chainU, about0.953Aaway; the supplied point does not coincide with a named terminal atom. The coordinate model's chainU contains75observed residues. Resolve terminal-atom numbering/missing atoms before using it as an atomistic transfer-site reference.
4. Defaults are: VHLdistance60A, CRBNdistance16A; lysine residue SASA threshold2.5; at least1accessiblelysine; occlusion distance5A; clash threshold1A and tolerance10clashingreceptoratoms. These are configurable geometric heuristics, not universal catalytic distances or atom-type-aware force-field thresholds.
5. The current clash function examines all parsed atoms without explicitly removing hydrogens; heavy-atom selection must be explicit in our adaptation. The lysine occlusion exclusion uses residue number without full chain/insertion identity. Check missing NZ atoms explicitly.
6. Accessibility execution is nested under the crl_model_clash conditional in this snapshot. Independent switch combinations require verification before modular reuse. These are source-review findings, not an assertion that the authors' standard curated examples fail.

## Recommendation for our campaign

P4ward is a strong candidate reference implementation for the planned static full-CRL analysis and may reduce custom development. Reuse the relevant model/filter layer after saved ternary docking; no wholesale migration or new MD is required by that layer. This request does not authorize switching or interrupting the running docking qualification.

Before integration, map our VHL/Elongin chain identities and residue correspondence to its references, verify the exact E2/Ub point/atom definitions, retain engineered mutations and missing atoms explicitly, and compare against our checked VHL experimental reference. Preserve pose coordinates and source hashes. Use full target context where reliable, recording missing domains/lysines. Adapt the input/output layer rather than forcing pRosettaC outputs into MEGADOCK rotation-object assumptions without verification.

Record continuous closest-lysine distances, exposure, occlusion and clash values for all poses; retain failures and clear denominators. Compare its default pass/fail rule with a productive fraction under the same sampled ensemble. For later GNN use, these can be declared node/graph features and a geometric baseline. Do not use geometric pass/fail as experimental degradation labels or train/test on independent poses from the same target.

The full-CRL analysis remains planned, not integrated. The immediate six-PROTAC Ernie preparation and launch-critical docking tests continue. The selected P4ward code/model snapshot is privately archived; only our review and provenance are committed here.
