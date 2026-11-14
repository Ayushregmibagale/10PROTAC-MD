# Rosetta implementation and sampling comparison

## Main conclusion

The hardened local adapter has stronger explicit error/provenance/chemistry safeguards than the public PROTAC_ternary helper. The published Bai2022 protocol has substantially broader and differently organized protein docking sampling than our typical retained campaign summaries. Neither source review nor pose count establishes that ours is scientifically better. A same-input native-recovery benchmark remains required.

## What their public code actually provides

Public repository pinned to a174357135d599c7274d24779e0662e15461b0ef. It provides the ternary linker-placement helper, model-modification/FFC helper scripts, examples and Rosetta command templates. The2022 article describes a modified three-stage Rosetta search, RDKit linker generation and CRL4A filtering. The checked repository does not contain a complete runnable driver implementing all those2022 stages or CRL4A filtering. The public helper is not PRosettaC.

## Pose budgets: compare the same stage

Bai2022:20 initial Rosetta global decoys → top10 ×1000 global decoys → top5000 clustered at3A →200 representatives ×1000 local decoys → top20000 for ternary construction,then minimization and clustering. Total scheduled docking decoys210020; this is the published schedule,not an independently verified completed-job ledger. The200 representatives are starting points for200000 local attempts,not their entire docking output.

Our current receipt wrapper defaults to200 accepted global arrangements and50 local refinements each,or10000 local attempts if all stages complete. The underlying PRosettaC script has500/1000 global targets in different modes,but the wrapper caps the effective count. Our global stage is random rigid-body sampling with anchor-distance and CA-clash filters,not Rosetta global docking or PatchDock. A high global-trial ceiling counts rejected proposals,not energy-refined poses. Recent50×5 diagnostic budgets must not be mistaken for the old full campaign.

Four actual historical summary files report:

- ABL1: 4080 local solutions; 2500 final models; 200 clustered. Raw docking directory present: False. These are summary-reported counts,not verified distinct pose counts.
- AXL: 2500 local solutions; 2500 final models; 200 clustered. Raw docking directory present: False. These are summary-reported counts,not verified distinct pose counts.
- DDR2: 2500 local solutions; 2500 final models; 200 clustered. Raw docking directory present: False. These are summary-reported counts,not verified distinct pose counts.
- EPHA2: 2500 local solutions; 2500 final models; 200 clustered. Raw docking directory present: False. These are summary-reported counts,not verified distinct pose counts.

Thus the published200000 local attempts are about80times a2500-local historical summary count. This is not80times more accuracy. Local/refined/final/clustered/MD-selected poses are different denominators. Our200 clustered models and~20 MD starts are downstream selections,not all sampled models. Historical score labels below-zero are retained as historical text,not evidence that current repaired raw energies use the same convention.

## Direct source findings

- Public helper alignment uses selected atom coordinates or group centroids with default RMSD cutoff0.4A. It does not itself prove complete ligand graph/stereochemistry preservation or complete bound-head agreement. PDB merging intentionally omits adding the incoming bonds; the separate modifier/parameterization stage is required. This design choice is not automatically an invalid molecule.
- Public helper has no explicit duplicate unqualified atom-label guard; the dictionary silently picks the later repeated label. Chain/residue-qualified supported formats can avoid this ambiguity.
- Public helper ignores nonexistent atom-deletion labels. It prints an error for negative cutoff but still exits0; NaN cutoff also exits0 with no models. Both can be confused with an empty geometric result unless guarded by a controller. Zero alignment iterations are not rejected as invalid. Malformed atom-format raises a Python TypeError because it raises a string.
- Our repaired linker explicitly checks full stereochemical disjoint head matching,both-head displacement,finite results andoutput stereo. Our parameterization preserves named atoms,explicit H andconformer frame. Its finite ETKDG retry budget still cannot prove that rejected poses are chemically impossible.
- Our receipt wrapper checks nonzero exits,missing outputs,stage accounting,hash drift andreplay provenance. Existing141 regression tests passed in the analysis Python environment; this does not establish scientific docking accuracy.
- Our method named relax currently performs PackRotamersMover sidechain repacking. It does not call FastRelax or Bai minimize_ppi with jump_all. That is an important scientific-method difference despite the method name.
- Our local refinement calls DockingLowRes plus DockMCMProtocol through PyRosetta. This is not established as identical to the paper's Rosetta3.12 docking_protocol configuration. Current validated local runtime is a different Rosetta release.
- Our final ranking selects up to1000 by final raw score,then200 by a local total ref2015 score,and clusters at4A using receptor-frame coordinates. Public Bai methods use3A clustering and interface-energy ranking; these are not interchangeable score definitions.
- The local adapter top-level source andruntime copy differ in resource controls: execution uses the runtime copy with an environment-enforced trial cap and ledger. Auditing an unused duplicate would give the wrong effective behavior.

## Runtime and verification

Unmodified upstream example1 generated one model with alignment RMSD0.234371A in the preserved Python3.10 docking overlay. Edge-case probes are saved in code_audit/probes. PyRosetta import succeeded; querying a nonexistent module __version__ attribute failed and was not a docking failure. The overlay lacks pytest,so the existing local regression suite was run using analysis Python3.13:141passed in10.64seconds,with dependency deprecation warnings. No production docking or model training was run for this source audit.

## Benchmark status

The earlier authorized VHL–MZ1–BRD4 benchmark is still pending. Source inspection and a published example are not native recovery. The complete2022 driver/constraint inputs have not been located in the public repository; an exact reproduction cannot be claimed. Any reconstructed PyRosetta comparison must be declared an adaptation,with a frozen equal budget and independent evaluation. No SJF8240 panel rerun is justified by this audit alone.

## Sources

-2022 methods: https://pmc.ncbi.nlm.nih.gov/articles/PMC9019245/
-2021 foundation and code availability: https://pmc.ncbi.nlm.nih.gov/articles/PMC8866032/
-Public code: https://github.com/karanicolaslab/PROTAC_ternary
-Our exact reviewed files and hashes: code_audit/SOURCE_MANIFEST.json

Original trajectories,poses,labels andmodels are unchanged. Only isolated audit outputs are new.
