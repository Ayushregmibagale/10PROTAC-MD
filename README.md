# 10PROTAC-MD

Pre-Ernie remediation is in progress. **Production submission is blocked.**
The legacy `prosettac_inputs/` are preserved for provenance and are not validated
inputs. They must not be submitted using an old copy of the launch script.

## Verified chemical identities

Tocris catalog structures are the reference molecular graphs:
GSK215 (7818), SJF8240 (7266), SJF1521 (7261), JNJ1013 (8029), SJF_delta (7267).
`config/molecules.json` records supplier isomeric SMILES, formula, InChIKey,
CAS and source URLs. `molecules/verification.json` records machine checks.
SJF_delta is Tocris catalog 7267; its distinct VHL attachment is retained in the registry.
These five identified compounds are not a finalized ten-compound study.

## Repaired preparation and analysis

* `00c_verify_supplier_registry.py`: validates chemistry and writes separate full SMILES.
* `00d_prepare_bound_references.py`: extracts deposited ligands using CCD bond graphs
  and native coordinates; normalizes label chain IDs explicitly.
* `00e_prepare_gsk215_heads.py` and `00f_prepare_foretinib_heads.py`: derive exact bound head fragments and anchors.
* `00g_model_modified_heads.py`: models changed head atoms with a fixed bound scaffold; requires explicit review.
* `01_place_warhead_batch.py`: requires an explicit chain/domain table, maps protein
  sequence or CE structural correspondence, checks pocket fit and writes separate preparation candidates.
* `01b_preflight.py`: checks the runner contract, molecular identity, head substructures,
  chain selections, contact geometry and preparation hashes.
* `00j_prepare_delta_p38_control.py`: prepares a separately labeled p38δ control
  using the published 5IA4/4EYJ template pair and explicit KLIFS pocket mapping.
* `02_batch_prosettac.py`: read-only planning, atomic manifests and a bounded
  foreground pilot controller for the repaired original baseline. Full campaign
  mode requires versioned release evidence; remote submission remains disabled.
* `03_record_attempt.py`: isolated local benchmark harness, timeout handling and terminal
  receipts after structural validation. Not a production submitter.
* `01_extract_bridging_features.py`: computes geometry from receipt-validated PDBs.
  It does not infer distances from Rosetta score columns or claim ML compatibility.

## Local validation

```bash
python3 -m pytest -q tests
python3 scripts/01_setup/00c_verify_supplier_registry.py
python3 scripts/01_setup/02_batch_prosettac.py --plan
python3 scripts/01_setup/01b_preflight.py --output audits/legacy_preflight.json
```

The final command is expected to exit 1 while legacy inputs remain invalid.
Use `--inputs prepared_inputs` to inspect new candidates; these remain blocked by
`BENCHMARK_RELEASE_PENDING` even if their mechanical checks pass.

## Remaining release gates

Representative global recovery and cross-docking controls; sampling convergence;
justification or improvement of restricted kinome transfer coverage; assay-linked
protonation/state selection; Ernie pilot resource measurements and writable storage
allocation; converged MD equilibration and frozen statistical analysis. No production jobs have been launched by this repair.
See `docs/REMEDIATION_STATUS.md` for current scope and limitations.

## Reproduction and data retention

Code, exact supplier identities, compact reference/control inputs, ligand parameter
sets and selected validation receipts are version controlled. Large diagnostics,
legacy inputs, complete structure inventories and simulation outputs remain in the
local project; their exclusion from Git does not delete them. Full raw retention is
the agreed planning default (approximately 7.24 TB for the measured 112-pair docking
example, excluding MD, retries and backups). See docs/DATA_RETENTION.md.

Install Python dependencies from requirements-qc.txt for the software checks.
Licensed PyRosetta and AmberTools/Gromacs are separate runtime requirements.
The pinned Rosetta helpers can be obtained with baseline_engine/fetch_rosetta_helpers.py.
The current checks and open scientific requirements are recorded in
docs/REMEDIATION_STATUS.md; complete changes and failed attempts are recorded in
docs/REMEDIATION_LOG.md.

The separate explicit-pocket comparison (`01d_compare_pocket_screen.py`) records
161 candidates and all five known-target controls. Its source correspondences
are built by `00k_cache_klifs_pockets.py` and `00l_build_pocket_correspondence.py`.
It remains comparison-only; the controller does not silently replace the original
112-pair cohort. See docs/NATURE_READINESS.md for the requested publication standard.
