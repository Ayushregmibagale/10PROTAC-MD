# Saved workflow decision: full CRL2 ubiquitination assessment

User explicitly asked to retain this on16November2026. Immediate priority remains the six-PROTAC Ernie campaign. No new MD, broader literature project or GNN work is required before docking submission.

## Bai2022 workflow to learn from

Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC9019245/ . Their own Rosetta-based protocol generated target–PROTAC–CRBN ensembles with linker-derived constraints and clustering. The full CRL4A machinery was modeled separately: CRBN/DDB1/CUL4A/RBX1/NEDD8/E2/ubiquitin. They built12machineryconformations and retained5for their main analysis. Ternary-cluster representatives were superimposed through CRBN onto the larger assemblies. Exposed target lysines near the ubiquitin C-terminus, with no incompatible clashes, defined productive geometries; their approximately16Acriterion was a study-specific CRBN criterion. Productive-model fractions were compared with experimental ubiquitination, including selected lysine-mutant tests.

Do not describe this as docking the whole CRL4A from scratch for each pose, an exact pRosettaC implementation, a universal distance law or direct proof of cellular degradation.

## Our retained plan

1. Generate qualified target–PROTAC–VCB ensembles. VCB=VHL+ElonginB+ElonginC. This is the current six-drug docking scope; full catalytic machinery is not included in docking.
2. Separately prepare and validate full CRL2-VHL catalytic references containing VHL/ElonginB/ElonginC/CUL2/RBX1/NEDD8/E2/ubiquitin. Verify chain identities, atom/residue mappings, engineered mutations, missing atoms and reference-conformation assumptions.
3. Align each retained ternary pose through VHL to the checked full-machine reference ensemble. Reuse saved poses; a new whole-machine docking search or new MD is not automatically necessary.
4. Assess clashes against the larger machinery, target lysine exposure/burial, orientation and distances to E2/ubiquitin. Report construct coverage and missing full-protein context. Keep catalytic-reference sensitivity explicit.
5. Compute clearly defined productive fractions with rejection accounting and a fixed denominator. Do not transplant Bai's16Acutoff or interchange Ub-terminal distances with E2-Cys/CA proxies without VHL-specific checks. Use continuous measurements and prespecified sensitivity analyses.
6. Later compare these mechanistic baselines with a compact static pose-ensemble GNN using compound-specific experimental labels. Do not turn predictions or docking failures into experimental labels; SIK3/RPS6KB2 remain unknown.

Status distinction: this is the retained downstream analysis plan, NOT a claim the full CRL2 assessment is already integrated or validated for the six compounds. Current launch-critical work is generator/input qualification and a concrete Ernie pilot package. No remote submission is authorized by this memory request. Ask for final submission approval when the concrete package is qualified and ready.
