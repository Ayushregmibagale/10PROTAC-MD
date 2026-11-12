# SJF8240 retained-structure GNN experiments

## Folder-wide audit and supported reuse path

The 12 November audit found critical regressions in exploratory chain selection, contact definitions and validation. Twelve exploratory entry points are retired with original sources archived. [Audit report](folder_audit_20261112/AUDIT.md) · [Word copy](folder_audit_20261112/AUDIT.docx) · [Reusable workflow](folder_audit_20261112/robust_analysis_v1/README.md).

Use `RUN_ANALYSIS.py` and per-compound manifests for new work. All148 reviewed saved models passed fresh CPU replay; the configurable contact extractor reproduced141 SJF8240 cases and16 role graphs. The documented verification route has76 passing tests, including a portable amendment for a historical test-path failure. Other compounds still need their own verified labels, protein identities, chemistry/MD provenance and feature calibration; similar results are not assumed. Frozen SJF8240 studies below remain preserved.

Latest completed analysis: [matched holdouts and lysine-aware GNN](mechanistic_v3/RESULTS.md). All 48 weights independently replayed on CPU. Original versus enriched matched-pair ranking: VHL 7/12 to 9/12 (0.583 to 0.750); ELOB 6/12 to 1/12 (0.500 to 0.083). Specificity remains zero at 0.5 in both enriched arms; VHL balanced accuracy is 0.458. This is a limited ranking gain, not useful binary discrimination or external validation.

Per-lysine SASA analysis covers 141 existing poses and 5,781 sampled frames. ABL1 bridge-frame exposure is approximately 98.4%; the earlier claim that all 18 lysines were buried was unsupported. SIK3/RPS6KB2 remain experimentally unknown and excluded from these new training/evaluation folds.

Complete step-by-step history: [Word document](mechanistic_v3/docs/STEP_BY_STEP.docx) or [Markdown](mechanistic_v3/docs/STEP_BY_STEP.md). [All held-out scores](mechanistic_v3/results/target_held_out_scores.tsv). No new docking/MD; previous studies preserved.

## Previous hypothesis cohort

Previous: [corrected_md_only_v2 hypothesis test](corrected_md_only_v2/RESULTS.md), using the exact v1 MD-only architecture, hyperparameters, 300 epochs, base seed 42 and unweighted loss. The six degraded labels are unchanged. SIK3_full and RPS6KB2_full are added as assumed negative targets alongside ABL1/AXL, giving 6:4. Their real experimental degradation labels remain unknown. No graphs were geometrically rebuilt and no new docking/MD was run.

| v2 arm | AUROC | Balanced accuracy | Sensitivity | Specificity |
|---|---:|---:|---:|---:|
| ELOB | 0.0417 | 41.7% | 83.3% | 0.0% |
| True VHL | 0.1667 | 12.5% | 0.0% | 25.0% |

Both SIK3/RPS6KB2 were genuinely excluded from their own training folds, yet both arms predict both as degraded at threshold 0.5. This disagrees with their assumed negative labels and does not establish actual degradation. [Held-out score table](corrected_md_only_v2/results/hypothetical_target_held_out_scores.tsv).

Both v1 arms had AUROC 0.0 and balanced accuracy 50%. V2 does not rescue classification. The report includes full ten-target metrics, the original-eight subset and pair-count decomposition. Adding targets changes class balance, geometry/composition and updates per epoch together; this test cannot prove either general geometric discrimination or its universal impossibility. Independent CPU replay verified all 20 saved models.

The confirmed/user-labelled eight-target map remains config/current_md_labels.json: six positives and only ABL1/AXL as confirmed negatives. The two hypothetical labels are scoped exclusively to corrected_md_only_v2/config/user_labels.json and are not new experimental validation.

Earlier corrected_md_only_v1, md_weighted_v1, md_weighted_e2_v1 and fix_comparison_20261111 remain unchanged. The original 16-target dual-encoder history is also preserved. Current status: STATUS.json.
