# Isolated pRosettaC robustness work

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
