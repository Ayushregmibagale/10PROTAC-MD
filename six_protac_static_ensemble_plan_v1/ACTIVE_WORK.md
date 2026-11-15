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
