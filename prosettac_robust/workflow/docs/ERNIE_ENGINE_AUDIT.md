# Ernie read-only engine audit — 2026-11-04

SSH access succeeded after the network permission was approved. No remote files
were written and no molecular jobs were launched. Code snapshots and provenance
are under `engine_audit/ernie/`; they are audit artifacts, not a released engine.

## Observed host and environment

* Host: `ernie`, 96 logical CPUs, about 1,031,380 MiB total RAM.
* No sbatch/qsub/srun command was found on the SSH PATH; checked common absolute
  scheduler locations were also absent. This does not prove that no resource-use
  policy exists. The inspected runner uses a LOCAL thread/subprocess backend.
* Interpreter: `/home/s1638308/miniconda3/envs/protac/bin/python3`, Python 3.10.20.
* Metadata versions: NumPy 2.2.6, SciPy 1.15.2, RDKit 2025.9.5, pandas 2.3.3,
  PyRosetta 2026.3+releasequarterly.5e498f1409. openbabel-wheel distribution metadata
  was absent; this is not proof that the OpenBabel module is unavailable.
* `/home` was 81% full: about 5.7 T used and 1.4 T available on a 7.0 T filesystem.
  These are observation-time values, not a reservation for this project.

## Exact code findings

The runner hash matches the initially inspected local runner. The remote adapter,
main and resume script differ from their local counterparts; their actual contents
were therefore retrieved and hashed instead of assuming equivalence.

1. `pyrosetta_adapter.relax` catches loading/scoring RuntimeError, copies an
   unrelaxed structure and writes score -1.0. It can also catch repacking failure
   and continue with the original pose. Scientific failures are not propagated.
2. `LOCAL.runBatchCommands` prints failed subprocess return codes but returns
   `local_done`; exceptions from workers are printed and suppressed. A successful
   parent exit cannot establish that every child succeeded.
3. `resume_pipeline.score_docking_pdbs` assigns 0.0 after scoring exceptions.
   It also shifts score.sc for clustering and keeps a score.sc.orig/.shift file.
   Raw and shifted scores must not be mixed or compared without reconciliation.
4. The adapter disables RDKit chirality enforcement in two embedding paths.
   Complete output PROTAC identity and 3D stereochemistry must be checked, not
   inferred from the input SMILES or heavy-atom counts.
5. Global docking substitutes random rigid-body sampling for PatchDock. A local
   default_rng(42) does not seed the separate global NumPy RNG used for rotations.
   Fixed seed fragments do not establish complete campaign reproducibility.
6. Remote main fixes one absolute constraint-generation path problem, but still
   contains unguarded os.system calls and zero-exit error branches. That difference
   alone does not repair the other execution contracts.

The repaired local QC now rejects possible placeholder/shifted scores, incomplete
pose sets, failed/missing stage accounting, incorrect full-PROTAC connectivity or
3D stereochemistry, and modified completion artifacts. An engine_status.json stage
report is required; the current remote engine does not produce that report and
therefore cannot be certified by the new benchmark harness yet. This is deliberate
failure detection, not a claim that the remote engine has been repaired.

## Work required before engine release

Create a versioned isolated engine copy; propagate worker and scoring failures;
preserve raw scores; explicitly seed all relevant generators/workers; enforce and
verify chirality; report attempted/succeeded/rejected/failed counts for every stage;
test interruption, resume, no-hit and partial-failure cases. Benchmark known native
complex recovery and independent cross-docking. Validate behavior in this exact
Python/PyRosetta environment before any production launch.

Storage also needs a measured pilot and retention plan: millions of uncompressed
combined PDBs can exceed the observed free space. Requested pose counts are upper
sampling targets, not guaranteed valid output counts; budget from actual pilot
files and include MD trajectories, retries and backups.
