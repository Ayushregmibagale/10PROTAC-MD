# ML transfer scope and validation requirements

The original 10PROTAC-MD README proposed ten VHL PROTACs, proteome-wide screening
and transfer of the protac-md XGBoost/GNN models. The current verified registry
contains five compounds, and the current structural inventory is a kinase-domain
inventory, not a complete proteome. The remaining compound list and the scientific
scope must be resolved before a ten-compound/proteome-wide claim is made.

The original protac-md project was inspected read-only. No model, input or result
there was modified. Reference files and their hashes are recorded locally in
engine_audit/protac_md_ml_readonly_audit.json. Historical LOPO results describe a
small target panel; they are not evidence of accuracy on the five new compound
ensembles. The currently inspected graph-metadata file contains unknown labels
(-1), which must not be converted into experimentally negative examples.

The existing MGTbind table has 320 structural rows, including degradative and
non-degradative mechanism labels. A structural row is not automatically an
independent compound/assay measurement. Mechanism categories, unreported
responses, censored DC50/Dmax values and confirmed assay negatives must remain
distinct. Parent-inhibitor binding is not a cellular degradation label. No such
reference data were silently imported as validation ground truth here.

The repaired generator changes stereochemical validity, head placement, linker
feasibility and rejection handling. Receipt-derived geometry features therefore
use an explicitly new schema with old_model_compatible=False. Conditional
contacts among both-head-constrained accepted poses cannot be interpreted as a
calibrated degradation probability. The full search/rejection denominator is
retained instead of relabeling failed preparations or searches as biological zeroes.

Before any model transfer: identify the exact frozen training dataset and feature
schema; reproduce its preprocessing; document compound, target and assay labels;
validate complete required MD/E2/graph features; assess applicability to the new
chemotypes and generator; and evaluate genuine experimental holdouts. Do not
fill unavailable features with zeros or report docking scores as measured affinity.

For new cross-compound claims, split by the relevant independent biological unit.
Related poses, structure variants, AF2/crystal alternatives and repeated records
of the same compound/target must not leak across training and validation folds.
Protein-held-out performance alone does not establish compound-held-out transfer.
Calibration and uncertainty must be evaluated on data appropriate to the claimed
use case. These are study requirements, not results already achieved by software QC.

No new ML training, prospective prediction table or validation accuracy has been
fabricated during remediation. The current software produces validated structural
and diagnostic outputs; the scientific model-transfer release remains pending.
