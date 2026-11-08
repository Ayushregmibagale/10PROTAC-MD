# Historical 19-protein panel

These are the confirmed saved-project groups: 17 historical GNN evaluation targets plus two prediction-only targets. Group names come from the retained configuration, not a newly validated quantitative affinity assay. High/low refers to the historical affinity grouping, not degradation magnitude. Directory placement and the `_full` suffix do not establish measured affinity.

| Protein / saved identifier | Retained project group | Saved experimental-class label | Qualification |
|---|---|---|---|
| ABL1 | High-affinity non-degrader folder | Non-degrader (0) | Reviewed proteomics and western endpoints disagree; saved label remains unchanged. |
| AXL | High-affinity non-degrader folder | Non-degrader (0) |  |
| CDK17 | Low-affinity degrader folder (legacy) | Non-degrader (0) | Legacy folder says degrader, saved label is 0; retained source is CDK14/O94921, not CDK17/Q00537. |
| CDK2 | Low-affinity non-degrader folder | Non-degrader (0) |  |
| CDK4 | Low-affinity non-degrader folder | Non-degrader (0) |  |
| CDK6 | Low-affinity non-degrader folder | Non-degrader (0) |  |
| DAPK1 | Low-affinity non-degrader folder | Non-degrader (0) |  |
| DDR2 | High-affinity degrader folder | Degrader (1) |  |
| EPHA2 | High-affinity non-degrader folder | Non-degrader (0) |  |
| EPHB2 | High-affinity degrader folder | Degrader (1) |  |
| IRAK1 | Low-affinity non-degrader folder | Non-degrader (0) |  |
| MAP4K5 | High-affinity non-degrader folder | Non-degrader (0) | Saved label conflicts with the particular reviewed proteomics endpoint; do not resolve from predictions. |
| MAPK14_full | High-affinity degrader folder | Degrader (1) | _full identifies the retained computational construct; not an extra biological target or affinity measurement. |
| MET | High-affinity degrader folder | Degrader (1) |  |
| RIPK2_full | High-affinity degrader folder | Degrader (1) | _full identifies the retained computational construct; not an extra biological target or affinity measurement. |
| RPS6KB2 | Prediction extension; affinity group unassigned | Unknown / not in evaluation | Prediction only; a NON call is not experimental negative evidence. |
| SIK3 | Prediction extension; affinity group unassigned | Unknown / not in evaluation | Prediction only; a NON call is not experimental negative evidence. |
| SLK | High-affinity non-degrader folder | Non-degrader (0) |  |
| TLK2 | Low-affinity non-degrader folder | Non-degrader (0) |  |

Totals: five saved degrader labels, twelve saved non-degrader labels and two unknown prediction targets. The retained evaluation therefore has no independently confirmed low-affinity-degrader subgroup simply because CDK17 remains in that legacy folder.

TNIK, PIP4K2C and alternate MAPK14/RIPK2 constructs appear in older project material but are not extra members of this user-confirmed 19-target set. Historical duplicate/artifact constructs must not be counted as independent targets.

For biological interpretation, freeze compound, binding/degradation assay, cell line, concentration, exposure and threshold independently of the model. The saved directory categories do not settle missing/conflicting labels or compound correspondence. See `INTEGRATION_REVIEW_RESULTS.md`, `sjf8240_review/primary_assay_label_audit.json` and `existing_result_audit/integration_19_v2`.
