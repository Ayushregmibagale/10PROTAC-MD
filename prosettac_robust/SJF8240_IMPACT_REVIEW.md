# Preserve completed SJF8240 results; trace impact before rerunning

**Latest scope and evidence:** retain the exact old structures; focus on scripts/rationale. Supplier replacement is not the active task. The early next-step statements below are historical; actual connectivity has since been traced in 174 compiled production TPRs, and all 174 protein sets have been identified. Start with `HISTORICAL_METHODS_AND_RATIONALE.md`, `SJF8240_STATUS.md` and the saved evidence inventories.

## Established existing work

The read-only original `protac-md/analysis/gnn/e2gnn_lopo_summary.tsv` records
17 folds, 15 correct, accuracy 0.8824 and overall AUROC 0.8 for the E2-scalar GNN.
This is an existing cross-validation result, not a newly independently verified
performance claim. The corresponding script consumes `ensemble_graphs.pt`,
`md_contact_graphs.pt`, and `analysis/features/ml_table.tsv`. None was changed.

## Concrete input finding

The D: project's SJF8240.smi parses to InChIKey
`MCJFROJTDORRBB-ZWUDBOSPSA-N`; the supplier-verified registry graph is
`SCUZPTNDZHEHLB-ZWUDBOSPSA-N`. Both have C58H65F2N7O11S and 79 heavy atoms,
but different connectivity. Formula or atom count checks cannot detect this.
Some original setup scripts reference the historical Ernie SJF8240.smi path,
but that alone does not establish the actual graph in every completed pose,
parameter set, MD system or GNN feature row. No blanket rerun is prescribed.

A further read-only check found the same historical connectivity in both standard
Mol2 source copies under `protac-md/ligand_params/SJF8240_3d*` (identical SHA-256
4d9ceb1754042a711dbe29bb1fdfe4abf8a9c7aa331d417075b85dd75780a539).
The GAFF-typed Mol2 could not be parsed directly by RDKit because its atom types
use GAFF labels; that parsing limitation is not evidence of a corrupt file.
The original `08_run_random_screen.py` points to this ACPYPE parameter directory.
Per-system topology inclusion and graph-feature lineage still require verification.
See `EXISTING_MD_PARAMETER_CHEMISTRY.json`. This strengthens the reason for tracing
impact, but does not by itself establish which completed runs used these files.

## Focused next steps

1. Trace each completed SJF8240 dataset's full-ligand graph, parameters and pose provenance. Compare actual artifacts, rather than infer chemistry from names or scripts alone.
2. Check retained final poses for exact stereo/connectivity, bound-head displacement and all-chain overlap. Distinguish hard failures from absent provenance and finite sampling outcomes.
3. Trace affected pose/MD IDs into static graphs, contact graphs, feature tables, labels, CV folds and saved predictions. Preserve old evaluations as versioned historical results.
4. Specify the minimum affected upstream reruns and dependent feature/GNN recalculations. If inputs and features are unaffected, retain the completed computation. Evaluate potential leakage and uncertainty separately from execution success.

## Repaired derivative already available

The separate workflow includes checked subprocess execution, process-group
termination, deterministic seeds, immutable attempt/source/input manifests,
strict chemical and full-complex geometry checks, constrained full-head linker
construction, invariant ligand rotamer orientation, preservation of actual raw
scores, and explicit zero-pose outcomes. It also checks MD stage completion and
provenance. These protections apply to new runs; they cannot retroactively prove
that every old result passed the same checks.

## Focused sampler and GNN repairs completed

The initial MCS anchor path was replaced by exact disjoint stereochemical matching,
a unique anchor-pair contract and one RDKit thread, with an immutable per-run
histogram record. The actual 200-conformer SJF8240 paired control produced identical
conformer bytes, range and histogram. The separate SJF8240/MET 5×5×2 pilot completed
without execution failure but yielded no valid pose; it remains a sampling result.

The isolated graph/GNN scripts now enforce verified compiled chemistry, production
receipts, explicit atom mappings and cohorts, periodic distances, observed joint
contacts, deterministic CPU folds and source/input hashes. Full notes and current
scientific requirements are in SJF8240_STATUS.md and SJF8240_REMEDIATION_LOG.md.
The historical MD/GNN results remain unchanged; a corrected training release has
not been issued. The completed GSK215 recovery diagnostics are supplementary
engine evidence, not a replacement for SJF8240 validation.
