# SJF8240 retained-structure GNN experiments

The latest requested interventions are complete: [class weighting and E2 comparison](fix_comparison_20261111/COMPARISON.md). All experiments use the user-corrected eight-target MD cohort: six degraded binders including EPHA2; ABL1/AXL are the two non-degraded binders. Every target has MD.

| Primary true-VHL variant | AUROC | Balanced accuracy | Specificity |
|---|---:|---:|---:|
| Corrected-label MD-only baseline | 0.0000 | 50.0% | 0.0% |
| Balanced class loss | 0.0833 | 25.0% | 0.0% |
| Balanced loss + audited E2 descriptors | 0.0000 | 16.7% | 0.0% |

Neither intervention provides useful discrimination in this small retrospective cohort. The protocol, labels, graph inputs, epochs (300), seed (42), optimizer and threshold were fixed; E2 adds two normalized graph scalars at the classifier head. No outcome-based tuning, score inversion, conventional predictive model, new docking or MD was run. All 32 newly saved models pass independent CPU replay, including class-weight and E2 normalization provenance checks.

Versioned reports: corrected_md_only_v1/RESULTS.md, md_weighted_v1/RESULTS.md and md_weighted_e2_v1/RESULTS.md. ELOB role controls, full target scores, metrics/confusion matrices and independent receipts are retained. The comparison reports identical primary raw graph hashes across versions.

Current class definitions are in [CURRENT_LABELS.md](CURRENT_LABELS.md). The immutable label-map JSON records provenance when created; STATUS.json gives current completion. The original 16-target dual encoder, old EPHA2=0 labels and all earlier results remain frozen under the top-level config/protocol.json and historical report folders. They are not current corrected-label performance results.
