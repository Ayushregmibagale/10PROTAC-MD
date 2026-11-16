# Original-engine six-PROTAC Ernie submission package

Latest user decision, 16 November 2026: use the original PatchDock-based pRosettaC engine. Replacement-engine structural recovery is no longer the launch gate. Stop replacement-engine development and further sampling benchmarks. Preserve that work separately.

## What is prepared

`bundle/` contains an executable six-task input bundle: GSK215/PTK2, SJF8240/MET, SJF1521/EGFR, JNJ1013/IRAK1, SJF_delta/MET and MS143/AKT1. All six supplier input identities and source hashes are reused from the existing staged package. SJF_delta/MET is a parent-head technical control, not a newly declared biological positive. This is six technical tasks, not 600 prepared kinases per drug.

The wrapper calls `/home/s1638308/PRosettaC/main.py` directly. It does not import the PyRosetta adapter or replace PatchDock, local docking, linker modeling, restraints or scoring. A small LOCAL scheduling override bounds concurrent child commands and propagates subprocess failures previously swallowed by the existing LOCAL class. Original inner `os.system` calls remain inherited; output checks cannot prove all intermediate commands succeeded. No whole-pipeline robustness or scientific accuracy certification is claimed.

Tasks run sequentially, each with four local workers and one BLAS/OpenMP thread per child. A per-task 12-hour timeout and 50-GiB free-storage floor are operational limits, not runtime estimates or a resource allocation. No memory reservation is established. Existing attempts are preserved. Completed task scores and clustered PDB hashes are checked before skipping; failed tasks require explicit `--retry-failed` and get a fresh attempt directory. No within-task resume is claimed. `STOP` terminates the active process group and prevents the next task. Logs and task receipts stay in each attempt directory.

The bundle deliberately uses the **original short protocol**, `Full: False`: up to 500 PatchDock results × 10 local poses, with the original linker-generation default (nominal up to 100 conformers). This is not the adapter's 50×10×20 protocol and is not a frozen final kinome sampling choice. Original `Full: True` uses 1000×50. `build_bundle.py --full` is required so that this choice is explicit. Do not silently substitute a custom N through monkey-patching.

## Live Ernie path result

See `ERNIE_PATH_AUDIT.json`, checked 16 November 2026. SSH and the Python executable work. The original source files, LOCAL backend and Rosetta XML/flag files are present at `/home/s1638308/PRosettaC`.

Missing or unsuccessful at the configured locations:

- `/home/s1638308/PatchDock/{patch_dock.Linux,buildParams.pl,transOutput.pl}`.
- `/home/s1638308/rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease` and the expected clean/parameter helper scripts.
- `babel` in the protac environment (only `obabel` was found), Python 2.7 on the checked PATH, and `pymol` import in the configured Python environment.

A broader read-only file search under `/home/s1638308`, `/opt` and `/usr/local` to depth seven found PatchDock benchmark scripts but no matching PatchDock/Rosetta executables. This is a bounded search, not proof no shared installation exists. No sbatch/qsub/srun was found on the checked noninteractive PATH. This package uses the observed LOCAL route; site resource policy is still applicable.

The missing paths are actual execution blockers. They are distinct from replacement-engine scientific benchmarks. Provide actual installation roots through `--engine-config` when rebuilding; do not supply the PyRosetta package path as a standalone Rosetta root. Original Python-2-era helper compatibility and PyMOL availability require the real runtime. No automatic installation, license acceptance or download-form submission was performed.

## Historical engine provenance caveat

The current historical Ernie launcher `/home/s1638308/protac/scripts/02_rosetta_utils/03_run_local_prosettac.py` describes itself as “no Rosetta binaries, no PatchDock” and invokes `pyrosetta_adapter.apply_patches()`. Therefore these inspected files do not establish that the earlier eight-target SJF8240 outputs came from unmodified PatchDock/Rosetta. Existing results are preserved; no engine provenance or biological validation is rewritten. This caveat does not impose another replacement-engine benchmark on the requested original-engine submission.

## Commands

Build a new bundle locally (output must not already exist):

```bash
python3 six_protac_static_ensemble_plan_v1/ernie_original_engine_v1/build_bundle.py \
  --staged six_protac_static_ensemble_plan_v1/ernie_submission_v1 \
  --output /tmp/six_protac_original_release_v2 --full False \
  --engine-config /path/to/verified_ernie_engine_config.json
```

The optional engine configuration overrides installation roots and operational limits. Rebuild after changes so the manifest reflects them. Source and input hashes are recorded; no production transfer/submission has occurred. For a wider cohort, supply a staged INPUT_MANIFEST with all intended compound/target rows and actual input hashes; do not fabricate absent kinase pairs or reuse a head without its target-specific placement.

After the package is ready and final submission approval is obtained, transfer the bundle to the whitespace-free intended root `/home/s1638308/six_protac_original_campaign_v1`. On Ernie:

```bash
cd /home/s1638308/six_protac_original_campaign_v1
/home/s1638308/miniconda3/envs/protac/bin/python3 check_paths.py --config engine_config.json
bash submit_ernie.sh
# Only after successful path checks and final submission approval:
bash submit_ernie.sh --execute
# To select one task:
bash submit_ernie.sh --execute --task GSK215__PTK2_Q05397
```

Default invocation prints the plan and launches nothing. `--execute` checks original dependencies, frozen bundle hashes and inputs before starting any molecular task. Do not use the previous replacement-engine HOLD release gate for this distinct package. Original-engine runs still need molecular identity, binding-head/PTM preservation and output checks; this package does not assert the original cleanup preserves MS143 phosphoThr308.

## Checks completed

Five orchestration tests passed: changed inputs reject; missing outputs reject; nonfinite scores reject; original entrypoint dispatch works; child-command failure propagates. Shell syntax and six-task dry run passed. Tests use synthetic dispatch stubs and do not substitute for an original-engine execution. Actual molecular execution is blocked by the installation findings above, not by MZ1 recovery from the replacement engine.
