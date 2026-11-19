# Six-compound kinase release checklist

The request covers six compounds against 435 kinases, for 2,610 distinct pairs. The checklist has one row per compound and UniProt ID. Every row is held. Availability and the release decision are recorded in separate fields.

| Input state | Pairs | Meaning |
|---|---:|---|
| Staged candidate held | 433 | At least one candidate from the earlier frozen bulk draft with input hashes verified |
| Prepared candidate under review | 75 | Includes 57 earlier recovery candidates and 18 newly materialized DFG candidates |
| Modeled control verified | 1 | SJF8240/MAPK14, with 12 independently identity-verified retained poses in the technical pilot |
| Missing prepared candidate | 2,101 | No candidate recorded from the available preparation routes |

The 18 DFG candidates reproduced the benchmark placement and passed independent supplier-graph, head-attachment, chain, hash and geometry checks. The geometry screen was 1.5–4.5 Å nearest head contact and at most 0.6 Å worst element-radius van der Waals overlap. DFG residue-offset fits can still misalign the binding site, so these are structural-review candidates.

The held bundle selects one input per available pair. For multiple candidates it keeps the earlier staged selection first, with the verified MAPK14 control taking precedence for SJF8240/MAPK14. Files are hash-pinned. The bundled runtime contains the pilot-validated conformer serialization change and declared-head path audit. The previous bounded 20-global × two-local settings are carried over as provisional. The queue has 509 tasks and its release flag is false. The full 2,610-pair campaign would require completing the missing inputs and reviewing the available candidates before approving a new full release.

The primary SJF8240 technical control is MAPK14. SJF8240/MET retains its published target status and the zero-linker result from the tested docking geometries. The GSK215 alternate ranking is exploratory. SJFδ/MAPK13 ternaries remain under steric review. SJF1521's anchor mapping is documented and its modeled EGFR head remains pending production qualification. JNJ1013 and MS143 have full-budget repair evidence.

Rebuild the checklist with `build_release_checklist.py`. Input archives, raw molecular artifacts and licensed runtime files stay in the workspace and on the private host. The GitHub checkpoint contains the checklist, summary, preparation scripts, and compact receipts.

## Ernie dry run

The held package was staged in `/home/s1638308/ernie_held_509_candidate_bundle_v1`. Its release file SHA256 is `e35e2c3192c3c44892671ad5855f1425524865501a7b301bec804bc4ffe702ce`. Ernie verified all 509 tasks and their bundled runtime files. The queue rejected an attempted dispatch while held. No attempt directory was created.
