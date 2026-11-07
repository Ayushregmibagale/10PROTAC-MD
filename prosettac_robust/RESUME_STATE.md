# Resume checkpoint — 2026-11-07

Read `AGENTS.md`, `USER_REQUIREMENTS.md`, this file and `SJF8240_STATUS.md` first.

## Scope that must persist

Preserve/reuse the existing SJF8240 structures and completed PRosettaC/MD/GNN results. Harden scripts and audit rationales, explicitly the <5 Å bridging definition. Helpful bounded local diagnostics are authorized; a wholesale experiment rerun or Ernie production submission is not. Original `/home/ayushregmibagale/protac-md`, `/mnt/d/PROTAC degradation Prediction` and prior top-level work remain untouched. All derivative work belongs in `prosettac_robust/`. Keep committing source/evidence, and record what/how/why. No guarantee of Nature acceptance.

## Most important completed findings

- Compiled-TPR sequence audit covers all 141 core MD cases: historical group named VHL is actually Elongin B/Q15370, chain A. True VHL/P40337 is chain C, with 100% sequence identity and observed-sequence coverage in all cases. All POI groups match expected targets.
- One retained static structure per core target (17/17 representatives) has the same Elongin-B selection under the old VHL label. This is not an exhaustive static-pose audit.
- All 141 existing trajectories were re-analysed, not re-simulated, using frozen sources, two workers, common 0–20 ns windows and 201 frames/case. Stored GROMACS distances reproduce within 0.000022 Å. Original metrics/models remain preserved.
- At 5 Å, true-VHL contact occupancies differ substantially from the legacy Elongin-B descriptor. Tables and PDF/SVG/PNG figures are under `existing_result_audit/contact_core_v1/`. These are geometric occupancies, not degradation probabilities or proof of native pocket engagement.
- AXL/pose011 extends to 100 ns and has an identical duplicate at 20 ns. The derivative common-window audit explicitly handles that; originals are unchanged. Initial incomplete audits are retained.
- The corrected E2 script anchors to sequence-verified Elongin C (chain B, 85 Cα), consistent with topology identities. Do not blanket-call that transformation wrong. Bridging-conditioned subsets still use the legacy descriptor and need semantic/lineage review.
- Historical simulated ligand connectivity differs from the verified supplier reference; experimental labels also have missing/conflicting evidence. These remain separate scientific findings, not permission to silently relabel data or launch a replacement campaign.

## Software and bounded controls

- Isolated future linker method: centered ETKDG, retain experimental torsion preferences, disable basic-knowledge preferences; fixed heads/stereochemistry/QC preserved. Earlier controls/sources/failures are retained separately. No original study output was regenerated with this method.
- A broader-DG local integration control completed EM, 200 ps NVT, two 500 ps NPT stages and a 20 ps diagnostic segment. This is not scientific production/convergence validation.
- A selected-method minimization retry completed in 5,468 steps, Fmax 455.493 <500 kJ/mol/nm. Interrupted attempt retained. No following MD stages were started for that retry.
- Read-only XTC reader prevents offset/lock sidecars beside original files; tests cover normal and refreshed-offset modes.
- Future GNN evaluations now retain tensor-only per-fold weights, exact training-only scalar reference rows, held-out inputs/scores, hashes and full score precision. No historical SJF8240 GNN was retrained. Synthetic replay/one-epoch tests validate saving/restoration only.
- Analysis Python 3.13 supports the complete 122-test suite. The Python 3.10 docking overlay lacks PyG: 119 tests apply there, and 3 GNN tests are explicitly skipped. Final logs verify 122 passed in analysis and 119 passed/3 skipped in docking. Never claim 122 passed in both environments.

## Durable artifacts

`USER_REQUIREMENTS.md`: full user requests and preservation rules.
`SJF8240_REMEDIATION_LOG.md`: detailed what/how/why through section 29.
`existing_result_audit/BRIDGING_RATIONALE_AUDIT.md` and `METRIC_GLOSSARY.md`: rationale, primary sources, limits.
`existing_result_audit/protein_roles_core_complete/`: combined sequence proof for all 141.
`existing_result_audit/contact_core_v1/`: completed contact analysis, source snapshot, plan, per-case hashes/tables, figures.
`existing_result_audit/static_role_representatives_v1/`: 17 static representatives.
`sjf8240_review/`: original chemistry, assay, graph/model lineage audits and readiness metadata.
`validation_controls/`: software checks and installation receipts.

Persistent private docking overlay is `.runtime/python310-packages`; launch via `workflow/scripts/00_runtime/docking_python.py`. Base conda interpreter/fallback packages remain external dependencies. Licensed binaries stay out of Git.

## Staging / installation / Git

Staging root: `/tmp/10protac_finish/prosettac_robust`.
Installed derivative: `/home/ayushregmibagale/10PROTAC-MD/prosettac_robust`.
Guarded updater: `/tmp/update_isolated_prosettac.py`, plan `/tmp/prosettac_robust_update_plan.json`.
Generate plan, then apply with approved escalation. Source/destination hashes are checked; backups and private history retained; installed tests run. Do not bypass drift checks.
Repo: https://github.com/Ayushregmibagale/10PROTAC-MD, master. Read current `git log`/status for the actual checkpoint. Before the November 7 save, latest pushed commit was 3dfbdae and later work was partially staged; do not assume the prior promised final save had completed.

## Next work (no new MD)

1. Complete/verify/install/commit this checkpoint and preserve actual test logs and resume instructions.
2. Trace which E2 tables/scalars the saved GNN actually used. `09_build_ml_table.py` reads `e2_reachability_all_pois.tsv`, `e2_reachability_poi_summary.tsv` and other feature tables; `15_e2_reachability_corrected.py` writes an AXL-specific corrected table. Do not assume that corrected table fed the saved model.
3. Review bridging-conditioned E2 summaries and original graph semantics against sequence-identified roles. Preserve old predictions; any corrected feature/evaluation must be separately versioned and assay-reviewed.
4. Keep chemical identity, assay-label evidence, sampling bias, uncertainty and calibrated-probability claims explicit. Do not choose cutoffs or relabel targets to improve performance.

No project simulation/analysis batch remains running from the completed core re-analysis. Only checkpoint verification tests/installations may be active during this save. Do not automatically restart completed diagnostics on resume.

Final source verification logs and hashes are recorded in `validation_controls/software_checkpoint_v5.json`. The installation and Git-only export checks are additional saved evidence; inspect the latest installation receipt and Git history rather than assuming an earlier partial save was current.

## Additional completed E2/GNN audits on 2026-11-07

`existing_result_audit/gnn_scalar_lineage_v1` checks 68 scalar entries: 33 finite pairs (32 match/one mismatch) and 35 matching missing-value pairs. MAPK14_full/e2_delta_reach_80A is 0.022575 saved versus 0.02569375 by the declared paired estimator; unpaired available-means gives 0.02256875. This is compatible with an alternative estimator, not proof of a patch/history. Original values were not changed. See `E2_GNN_SCALAR_LINEAGE.md` for actual unconditional/static scalar definitions and limitations.

`gnn_availability_audit_v1` fits no model: saved GNN 15/17 correct, AUROC .8; MD availability flag 14/17, AUROC .875. These are retrospective diagnostics on recorded historical labels, not independent validation or proof of leakage/uselessness. See `GNN_AVAILABILITY_AUDIT.md`.

Next substantive work is E2 reference-frame/coherent-imaging and exact frame correspondence review on retained data, followed by final semantic/integration review. Do not repeat the completed 141-case contact batch. Current source tests after adding scalar-audit regressions: 124 passed in analysis; 121 passed with 3 GNN-only skips in docking. Latest Git history, software checkpoint and installation receipts identify which changes are durably saved.

The E2 static-frame audit is complete under `existing_result_audit/e2_reference_frame_v1`: all 141 starting Elongin-C frames agree exactly; the 13 identical 51-point arrays represent 51 unique coordinates, not 663 independent samples. A retained `16b_mapk14_e2_fix.py` explicitly patches the scalar with the alternative estimator, but no historical execution receipt proves when it ran. Do not execute that old script: it can upsert original tables and its residue-exclusion body contradicts its retention header. Dynamic imaging remains to be checked separately.
