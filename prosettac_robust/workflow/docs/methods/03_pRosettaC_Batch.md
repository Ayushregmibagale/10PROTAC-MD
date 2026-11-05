# Execution and completion contracts

Production submission is disabled pending verification of the exact Ernie engine,
scheduler/allocation requirements, and benchmark release. No direct nohup campaign
launch is available in the repaired script.

`02_batch_prosettac.py --plan` is read-only. `--build-manifest` performs a locked,
atomic union of discovered rows and existing rows, preserving unrelated filtered
entries. Changed input hashes invalidate prior state. Missing configs are included
as blocked pairs instead of silently omitted. Legacy completion flags are not
scientific validation receipts.

For an explicitly requested isolated local benchmark, `03_record_attempt.py` records
immutable inputs, command, UUID, start/end times and terminal status. Nonzero exit,
timeout, missing outputs, mismatched pose IDs, severe clashes, shifted scores or
possible -1 placeholder scores cause failure. Existing attempt directories cannot
be reused. A zero exit code alone is insufficient.

Rosetta scores are parsed from `Patchdock_Results/score.sc`, skipping SEQUENCE
and matching repeated headers. Every scored pose must have a matching combined
PDB. Completion is represented by `completion.json`, including all output hashes;
analysis rechecks it and the input hashes. A score file appearing during execution
never proves completion.

Geometry uses explicit VHL/POI/PROTAC selections from `selection.json`, heavy
atoms, and strict distances <5 A and <8 A. Fractions use 0-1; `bridge_pct` is
0-100. Missing or invalid data never becomes a zero-hit biological result.
The feature extractor accepts an explicit local attempt index; transfer must finish
before verification. No automatic partial rsync cache is treated as complete.

`n_dock` is a maximum requested global count, not a guaranteed number of valid
ternary complexes. Pilot sampling, independent seeds, runtime and storage must be
measured before the campaign budget is set.

The repaired controller now supports `--run --pilot`: foreground, sequential pairs,
explicit Python/seed and bounded global/local/conformer/worker/trial/time budgets.
Each attempt has a fresh UUID. Validated attempts with the exact same input and
protocol hashes are skipped; failed attempts are retained on retry. SIGTERM,
SIGHUP, timeout and keyboard interruption terminate the whole child process
group and leave a terminal receipt. Full mode requires source/input/evidence
hashes and a storage budget in a versioned release file, which has not been issued.

The original global/local docking stages remain. The former free conformer rigid
fit failed to keep both full ligand heads bound. Both-head-constrained ETKDG/UFF
now preserves those coordinates and rejects impossible linker geometries. Raw
scoring counts and the separate pose-QC acceptance counts are retained. A completed
search with no valid poses is COMPLETE_NO_VALID_POSES, not a biological negative.
