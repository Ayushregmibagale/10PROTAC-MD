# Check of proposed historical-adapter route — 16 November 2026

## Decision

The suggestion is right about the historical kinome launcher using the PyRosetta adapter. Reusing that installed runtime is a practical route and does not require installing standalone PatchDock/Rosetta or passing the replacement-engine MZ1 benchmark. However, the proposed symlink is not a demonstrated missing dependency, and the current multi-compound input files cannot be launched unchanged. Prepare a wrapper around the existing runner using compound-correct inputs; do not blindly expand the old shell script or modify historical results.

## Directly verified on Ernie

- `/home/s1638308/protac/kinome_scan/scripts/08_run_kinome_ernie.sh` calls `/home/s1638308/protac/scripts/02_rosetta_utils/03_run_local_prosettac.py` with `--n_dock` and `--n_workers`.
- The runner applies `pyrosetta_adapter.apply_patches()`. That function replaces PatchDock with its random-orientation search, standalone Rosetta with PyRosetta calls, and OpenBabel CLI functions with Python implementations. This establishes the inspected launch route; it is not a per-output immutable historical provenance receipt.
- Actual imports of RDKit, NumPy, SciPy, PyRosetta, sklearn and Python OpenBabel PASS in `/home/s1638308/miniconda3/envs/protac/bin/python3`. The existing runner's `--help` also executes successfully. Neither command launches docking.
- The historical source directory `/home/s1638308/protac/kinome_scan/data/kinase_pRosettaC` contains 461 parameter files. These are target inputs, not a claim there are exactly 461 poses or that all 461 finished successfully.
- `/home/s1638308/10PROTAC-MD/prosettac_inputs` is absent. The new local input tree is not already present there on Ernie.

Source SHA256 observed remotely:

| File | SHA256 |
|---|---|
| `08_run_kinome_ernie.sh` | `445329a04c64e1bfc9a83f1b55e66d5f5b3582402e3da01352942a5dcbde9aa8` |
| `03_run_local_prosettac.py` | `d4c7228472e9ce7516cd3ea496cea7a50115a0b9bdf985e872d6de112f650c20` |
| `pyrosetta_adapter.py` | `df7bf0f642fe517df2b1f57e6443d0b97ab9615c7a75733c1e1d269a789c9808` |

## Evaluation of the specific suggestions

**Missing installation claim:** Earlier searches covered selected paths and a depth-seven search of `/home/s1638308`, `/opt` and `/usr/local`, not the entire filesystem. Standalone binaries were not found there. They are unnecessary for this adapter route. Installation does not inherently require administrator rights; no setup-duration estimate is established by those searches.

**Suggested `babel` symlink:** Not needed by the normal patched adapter calls inspected here. Python OpenBabel is installed and imports. Bare `which obabel` outside the intended environment can select the wrong executable, and creating a symlink in the shared environment without checking consumers is unnecessary. No symlink was created. If an actual residual CLI call later fails, trace that call and address its specific compatibility needs.

**One test pair:** A sound operational next step, using a new isolated directory with verified molecule-specific inputs and the unchanged historical runner. Do not use the current GSK215 `prosettac_inputs` parameter file as-is. `N_DOCK=50` sets global results only: `Full: True` still requests 50 local poses per global result (up to 2500 local attempts), whereas `Full: False` requests 10. Record the real sampling configuration. A successful process and some PDBs establish execution, not automatically correct compound identity, valid geometry or degradation accuracy. Check identity, stage failures and readable outputs without imposing a new MZ1 benchmark on this operational pilot.

**Expand old batch script:** Reuse the scheduling organization and runner calls, but replace unsafe housekeeping and completion logic. The old script deletes Patchdock_Results and other intermediates, accepts any existing result PDB as completion, and suppresses failed background waits with `|| true`. These behaviors should not be copied into a campaign that must preserve raw scores, failures and restart provenance. Preserve existing runs; write new task attempts separately.

## Concrete input findings

Read-only audit of `/home/ayushregmibagale/10PROTAC-MD/prosettac_inputs` is reproduced by `check_existing_inputs.py`; results are in `INPUT_AUDIT.json`.

| Compound | Parameter files | Finding |
|---|---:|---|
| GSK215 | 650 | Every file points to SJF8240 SMILES and lacks `Protac`, `Anchor atoms`, `Full`. |
| JNJ1013 | 644 | Same wrong-compound SMILES reference and missing required keys. |
| SJF1521 | 644 | Same wrong-compound SMILES reference and missing required keys. |
| SJF8240 | 461 | SJF8240 path, but same missing required keys in this new tree. The historical Ernie input files have the actual required fields. |
| SJF_delta | 0 | Directory exists, no per-target parameter files found. |
| MS143 | 0 | Absent from this input tree; remains part of the requested six-drug campaign. |

All populated files inspected declare `Chains: ABC A`; review/normalize component identity using the already staged corrected inputs rather than silently accepting overlapping chain identifiers. This audit does not certify the target-bound warhead structures or shared molecule files.

These are real input problems, not new scientific recovery requirements. Renaming folders or fixing `babel` would not correct the drug identity. Existing six-task staged inputs in `../ernie_submission_v1/inputs` already provide molecule-specific files and explicit parameters; reuse those as the starting point for the adapter pilot. The full kinase collection must subsequently be generated from each supplier-correct drug/head template, with actual target placement and manifest checks. Do not silently drop MS143 because the legacy tree lists five directories.

## Recommended immediate action

Keep the installed historical runner, freeze source hashes, stage one supplier-correct GSK215/PTK2 pair in a separate remote task directory, and perform a bounded operational run using the user's existing pilot authorization only where applicable; no full campaign submission follows automatically. After the run, check molecule identity and actual successful outputs, preserve raw/shifted scores distinctly and measure cost. Then finish the resumable multi-compound wrapper, populate only verified task inputs and present the concrete campaign for final submission approval. No new MD, standalone-engine installation or further replacement-engine benchmarking is required by this route.

This turn performed read-only inspection and local documentation/audit generation only. No remote files, environment symlinks or molecular jobs were created.
