# Step 03 — pRosettaC Batch Submission

## Overview

All kinase × PROTAC pairs are submitted to Ernie (96-core HPC cluster) in two stages.
Stage 1 is a fast screen to identify any bridging hits; Stage 2 is a deep run on hits only.

---

## Pose counts

| Stage | Purpose | N_DOCK | Output poses | Flag |
|-------|---------|--------|--------------|------|
| Stage 1 | Screen | 50 | ~200 Rosetta-refined | (none) |
| Stage 2 | Deep | 1000 | ~20,000 Rosetta-refined | `Full: True` |

**Rationale:** matches the protac-md 19-POI training campaign (~20,000 poses per POI).
kinome_scan stage 1 used N_DOCK=50 → ~200 constrained poses; top-5 full runs used
N_DOCK=1000 + `Full: True` → ~20,000 poses (source: `25_run_top5_full_prostac.sh`).

---

## pRosettaC interface on Ernie

```
python3 RUNNER Protac_params.txt --n_dock N --n_workers W
```

- **RUNNER**: `/home/s1638308/protac/scripts/02_rosetta_utils/03_run_local_prosettac.py`
- **Conda env**: `s1638308/miniconda3/envs/protac` (has PyRosetta + RDKit)
- **Stage 2 params**: `Protac_params_full.txt` (adds `Full: True` and `Anchor atoms: 6 1`)

---

## Submission script

**Script:** `scripts/01_setup/02_batch_prosettac.py`

```bash
# Stage 1 — all pairs
python 02_batch_prosettac.py --stage1

# Check hits
python 02_batch_prosettac.py --filter-stage1

# Stage 2 — bridging hits only
python 02_batch_prosettac.py --stage2
```

### Parallelism

`MAX_PARALLEL = 10` concurrent jobs on Ernie, each with `N_WORKERS_PER = 2` Rosetta
workers. Total: 20 CPUs (~20% of Ernie's 96 cores — conservative to avoid queue penalties).

---

## Run counts

| PROTAC | N_DOCK source | Input count | Stage 1 ETA |
|--------|--------------|-------------|-------------|
| SJF8240 | kinome_scan (sym) | 461 | ~62 hrs |
| SJF_delta | kinome_scan (sym) | 461 | ~62 hrs |
| SJF1521 | KLIFS/AF2 | 644 | ~86 hrs |
| GSK215 | KLIFS/AF2 | 650 | ~87 hrs |
| JNJ1013 | KLIFS/AF2 | 644 | ~86 hrs |

Wall time per PROTAC at MAX_PARALLEL=10: ~8–9 hrs (8 min × N / 10).

---

## State tracking

Each run is tracked in `runs/manifest_stage{1,2}.tsv` with columns:
`protac, gene, uid, run_name, status, n_dock, bridge_poses, bridge_pct, submit_time, done_time, ernie_outdir`

Statuses: `PENDING → RUNNING → COMPLETE | FAILED`

---

## Completion detection

Script polls Ernie every 60 s for `{ernie_outdir}/score.sc`.
On completion, `score.sc` is rsynced locally and bridging poses counted.

---

## Notes

- Stage 2 `Protac_params_full.txt` is written locally and rsynced alongside other inputs.
  The original `Protac_params.txt` is preserved for auditing.
- `Anchor atoms: 6 1` is required with `Full: True` (matches top5_full kinome_scan format).
- SJF8240 and SJF_delta inputs are symlinks to `~/kinome_scan/data/kinase_pRosettaC/`.
  rsync follows symlinks (`-az`) so the actual files land on Ernie.
