# Isolated SJF8240 publication robustness work

**Latest experimental evidence:** user-supplied RPE1 LF/PS data have been inspected (500 nM, 12 h; heavy=newly synthesized pool). LF `_new` provides matched results for all19 targets; five historical positive targets show reduced abundance in that saved version. Multiple result versions and weak/ambiguous evidence remain explicit; no label release or model refit. See [RPE1_EXPERIMENTAL_EVIDENCE_REVIEW.md](RPE1_EXPERIMENTAL_EVIDENCE_REVIEW.md).


**Latest 2026-11-08 reporting check:** 141 analysis tests pass; docking 137 passed / four PyG-only skips. Subgroup AUROC now evaluates independently and records undefined single-class results. Retained E2 scores yield overall AUROC 0.800, MD-only 0.600; no-MD AUROC is undefined. The current dual-GNN TSV does not reconstruct the older 0.850 baseline. Accession checks confirm the CDK17 abundance row cannot be assigned to the retained CDK14 structure. See [GNN_REPORTING_AND_ASSAY_RECONCILIATION.md](GNN_REPORTING_AND_ASSAY_RECONCILIATION.md). No original data, labels or predictions changed.


**2026-11-08 update:** completed the full 141-pose coherent E2 imaging sensitivity and retained stage/checkpoint audits. All 174 production checkpoints reach 20 ns. All 22,753 true-heavy bridging frames agree across three imaging conventions; 16 differences are nonbridging. One EM force-tolerance exception is disclosed. Analysis tests: 139 passed; docking: 135 passed / four PyG-only skips. See [COHERENT_E2_REANALYSIS.md](COHERENT_E2_REANALYSIS.md), [METHODS_DRAFT_RECONCILIATION.md](METHODS_DRAFT_RECONCILIATION.md), [RETAINED_GNN_MODEL_PROVENANCE.md](RETAINED_GNN_MODEL_PROVENANCE.md) and [PROTEIN_PANEL_19.md](PROTEIN_PANEL_19.md). Original results unchanged; unconditional biological/predictive sign-off remains withheld. Earlier dated counts below are historical checkpoints.


Latest completed review: [INTEGRATION_REVIEW_RESULTS.md](INTEGRATION_REVIEW_RESULTS.md), with evidence and software checks; scientific release remains withheld for documented findings.

**Active scope:** use the exact historical structures/results and review scripts plus rationale. Start with [HISTORICAL_METHODS_AND_RATIONALE.md](HISTORICAL_METHODS_AND_RATIONALE.md). Supplier replacement and new production are not the current task.


Start with [SJF8240_STATUS.md](SJF8240_STATUS.md). The other compounds are deferred.

This folder preserves the existing SJF8240 campaign and GNN work. The original
`D:\PROTAC degradation Prediction`, `/home/ayushregmibagale/protac-md`, and the
previous 10PROTAC-MD workflow are not edited by this isolated work.

- `legacy_source/`: exact, hash-recorded code/document snapshots from the D: project, for inspection only. These archived scripts contain original paths and must not be used as launchers.
- `workflow/`: separate repaired workflow based on 10PROTAC-MD commit `9d66c47`, with its own tests, configuration, references and validated supplier inputs. Results belong under this copy's own output directories. Use the controller with an explicit `--root` pointing here.
- `SOURCE_MANIFEST.json`: original file paths and SHA-256 fingerprints.
- `LEGACY_SJF8240_CHEMISTRY_CHECK.json`: direct comparison of the D: input SMILES and the verified supplier graph; historical downstream use remains to be traced.
- `SJF8240_IMPACT_REVIEW.md`: focused assessment and minimum-rerun decision process.

This is a repaired derivative of the working PyRosetta adapter, not unmodified
upstream pRosettaC. Its inherited global search samples random rigid-body
arrangements subject to anchor-distance and Cα filters; it does not run the
PatchDock executable despite compatibility folder names. Recovery adequacy is
still under assessment. Existing MD, model weights and predictions are not
silently replaced, relabeled, retrained or declared invalid.

Run regression checks from `workflow/`:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m pytest -p no:cacheprovider tests -q
```

The tested docking environment and separately licensed dependencies are described
in `workflow/baseline_engine/README.md`. Passing these tests is software evidence,
not a guarantee of publication acceptance. No Ernie campaign is launched here.
