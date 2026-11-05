# 09 — Phase 2: N-lobe Fix and MD Candidate Selection (2026-06-22 to 2026-06-25)

## Overview

After Phase 1 (all 13 POIs completed pRosettaC + ML classification — negative result, LOO 0.46),
Phase 2 focuses on selecting the correct MD starting structures using the paper's exact workflow,
and re-running POIs whose combined PDBs were structurally invalid due to the N-lobe artifact.

---

## Paper's MD Candidate Selection Workflow

From the pRosettaC paper (Zaidman et al.), the correct procedure for selecting MD starting structures:

1. **Filter bridging poses**: For each pose in `Patchdock_Results/combined_*_0001.pdb`, compute the minimum distance from PROTAC (chain X) to VHL (chain A, res ≤ 103) and to POI (chain A, res ≥ poi_min):
   - **PASS**: PROTAC within 5 Å of VHL AND within 5 Å of POI
   - **MARGINAL**: PROTAC within 10 Å of VHL AND within 10 Å of POI
   - **FAIL**: otherwise
2. **Re-cluster bridging poses only** (PASS + MARGINAL): Kabsch RMSD on Cα atoms, hierarchical clustering with 4 Å threshold (average linkage)
3. **Rank clusters by size** (largest first)
4. **Select top 3 clusters** as MD candidates

Script: `scripts/analysis/recluster_feasible.py`

### Output
`feasible_md_candidates.csv` written to each POI's run directory:
```
rank, cluster_id, n_members, representative_pdb, quality
1, cluster_1, 88, combined_149_17_0001.pdb, PASS
2, cluster_2, 34, combined_107_45_0001.pdb, PASS
3, cluster_3, 21, combined_101_22_0001.pdb, PASS
```

---

## N-lobe Artifact Discovery and Fix

### Problem
For kinases with residue numbering starting below 104, pRosettaC's combined PDB builder
overwrites VHL chain A residues with kinase N-lobe atoms. This makes the VHL binding
cleft (res ~60–100) physically absent from the combined PDB. The PROTAC has no VHL to
dock to, so all "bridging" scores are meaningless.

**Affected POIs:**
| POI | PDB | Start res | VHL kept | Result |
|-----|-----|-----------|----------|--------|
| RIPK2 | 4C8B | 8 | res 1–7 only | VHL cleft absent — INVALID |
| MAPK14 | 3GCP | 4 | res 1–3 only | VHL cleft absent — INVALID |
| SLK | 3CQU | 21 | res 1–20 only | VHL cleft absent — SUSPECT |
| MAP4K5 | 2C5K | 1 | none | VHL completely absent — INVALID |

**Confirmed by:** PROTAC distance check showing VHL at 44.9 Å (should be <10 Å) while
Elongin B/C were at 1.2 Å — chain A was not VHL.

### Fix
Renumber affected POI PDBs so residues start at 200 (above VHL's 1–103 range).
Only residue numbers change; 3D coordinates are unchanged.

Script: `scripts/structure_prep/prepare_full_kinase.py`

```python
def renumber_pdb(src_pdb, dst_pdb, offset):
    for line in f:
        if line.startswith(('ATOM', 'HETATM', 'TER')):
            orig_resnum = int(line[22:26])
            new_resnum = orig_resnum + offset
            line = line[:22] + f'{new_resnum:4d}' + line[26:]
```

**New structures:**
- `RIPK2_protein_clean_full.pdb`: res 200–508 (offset +192)
- `MAPK14_protein_clean_full.pdb`: res 200–548 (offset +196)

**New run dirs on ernie:**
- `~/protac/pRosettaC/runs/high_affinity_degraders/RIPK2_full/`
- `~/protac/pRosettaC/runs/high_affinity_degraders/MAPK14_full/`

---

## recluster_feasible.py Results (Phase 1 runs)

Ran on the Phase 1 combined PDBs (generated with old structures). Results for valid POIs only:

### DDR2 (valid — res 553+, no N-lobe overlap)
- 88 bridging poses (PASS + MARGINAL) from 10,000 total
- 27 clusters formed at 4 Å threshold
- **Top 3 MD candidates**: cluster_1 (PASS), cluster_2 (PASS), cluster_3 (PASS)
- Saved to `DDR2/feasible_md_candidates.csv`

### EPHB2 (valid — res 615+, no N-lobe overlap)
- 0 bridging poses in Patchdock_Results (fallback: scanned Results/ cluster dirs)
- 5 bridging poses found in cluster directories
- 4 clusters formed
- **Top 3 MD candidates**: cluster_1 (MARGINAL), cluster_2 (MARGINAL), cluster_3 (MARGINAL)
- Saved to `EPHB2/feasible_md_candidates.csv`

### RIPK2 (INVALID — N-lobe artifact, res 8+)
- Fallback to existing clusters (>500 feasible poses found but results unreliable)
- **Existing clusters not valid for MD** — will be superseded by RIPK2_full re-run

### MAPK14 (INVALID — N-lobe artifact, res 4+)
- Fallback to existing clusters
- **Existing clusters not valid for MD** — will be superseded by MAPK14_full re-run

---

## MET Re-run

MET was previously run (May 2026) with the corrected Kabsch fix but had stale Results/.
Attempted partial restart: encountered "0 poses" bug (Bug 23) — needed to also delete score.sc.

**Final fix:** Full clean restart:
```bash
cd ~/protac/pRosettaC/runs/high_affinity_degraders/MET
rm -rf Patchdock_Results/ Results/ local.fasc
tmux new -s met
python ~/protac/scripts/rosetta_utils/run_local_prosettac.py Protac_params.txt --n_dock 50 2>&1 | tee ~/met_run.log
```

**Note:** MET was run with `--n_dock 50` (2,500 poses) instead of the standard 200 (10,000 poses).
This is inconsistent with other POIs. If the 2,500-pose run yields too few bridging poses for
robust clustering, MET should be re-run again with `--n_dock 200`.

Status as of 2026-06-23: 1200/2500 docking poses complete, ~13 minutes remaining.

---

## Standard Run Command (IMPORTANT)

**The canonical command for running any POI is:**
```bash
cd ~/protac/pRosettaC/runs/high_affinity_degraders/<POI>
tmux new -s <poi_lowercase>
python ~/protac/scripts/rosetta_utils/run_local_prosettac.py Protac_params.txt --n_dock 200 2>&1 | tee ~/<poi>_run.log
```

This produces **10,000 poses** (200 global PatchDock × 50 local Rosetta = 10,000).

**Do NOT use `run_priority1.sh` directly** — its default `N_DOCK=50` gives only 2,500 poses.
All 12 non-MET Phase 1 POI runs used n_dock=200.

---

## Run Commands for RIPK2_full and MAPK14_full (ernie)

```bash
# RIPK2_full
tmux new -s ripk2_full
cd ~/protac/pRosettaC/runs/high_affinity_degraders/RIPK2_full
python ~/protac/scripts/rosetta_utils/run_local_prosettac.py Protac_params.txt --n_dock 200 2>&1 | tee ~/ripk2_full_run.log

# MAPK14_full (separate tmux window)
tmux new -s mapk14_full
cd ~/protac/pRosettaC/runs/high_affinity_degraders/MAPK14_full
python ~/protac/scripts/rosetta_utils/run_local_prosettac.py Protac_params.txt --n_dock 200 2>&1 | tee ~/mapk14_full_run.log
```

**Expected runtime:** ~4–6 hours each (same as original runs with n_dock=200).

---

## Post-run Steps (after each POI finishes on ernie)

### 1. Rsync results to home
```bash
# From HOME:
rsync -av s1638308@ernie.bio.ed.ac.uk:~/protac/pRosettaC/runs/high_affinity_degraders/RIPK2_full/ \
  ~/protac/pRosettaC/runs/high_affinity_degraders/RIPK2_full/
rsync -av s1638308@ernie.bio.ed.ac.uk:~/protac/pRosettaC/runs/high_affinity_degraders/MAPK14_full/ \
  ~/protac/pRosettaC/runs/high_affinity_degraders/MAPK14_full/
```

### 2. Run find_feasible_poses.py to identify bridging poses
```bash
cd ~/protac
python scripts/analysis/find_feasible_poses.py RIPK2_full
python scripts/analysis/find_feasible_poses.py MAPK14_full
# Outputs: ~/ripk2_full_feasible.txt, ~/mapk14_full_feasible.txt
```

### 3. Run recluster_feasible.py to select top-3 MD candidates
```bash
python scripts/analysis/recluster_feasible.py RIPK2_full MAPK14_full
# Outputs: feasible_md_candidates.csv in each run dir
```

### 4. Same for MET after it finishes
```bash
python scripts/analysis/find_feasible_poses.py MET
python scripts/analysis/recluster_feasible.py MET
```

---

## SLK and MAP4K5 — Pending Decision

Both are non-degraders. SLK (starts res 21) and MAP4K5 (starts res 1) have the same
N-lobe artifact as RIPK2/MAPK14. However, since they are non-degraders and the Phase 1
ML result was already negative, the decision on whether to re-run these has been deferred.

Options:
1. Re-run with renumbered structures (same approach as RIPK2_full/MAPK14_full) — more
   accurate MD candidates but adds ~2 more compute runs
2. Keep Phase 1 results for SLK/MAP4K5 — acceptable if the goal is only to validate
   the degrader MD simulations

---

## Key Files Created / Modified This Session

| File | Action | Purpose |
|------|--------|---------|
| `scripts/structure_prep/prepare_full_kinase.py` | Created | Renumber RIPK2/MAPK14 PDBs to res 200+ |
| `scripts/analysis/recluster_feasible.py` | Created | Paper's MD candidate selection (filter → cluster → top 3) |
| `scripts/analysis/find_feasible_poses.py` | Existing | Scans combined PDBs for PASS/MARGINAL poses |
| `data/structures/POI/high_affinity_degraders/RIPK2/RIPK2_protein_clean_full.pdb` | Created | Renumbered RIPK2 (res 200–508) |
| `data/structures/POI/high_affinity_degraders/MAPK14/MAPK14_protein_clean_full.pdb` | Created | Renumbered MAPK14 (res 200–548) |
| `pRosettaC/runs/high_affinity_degraders/RIPK2_full/Protac_params.txt` | Created | pRosettaC params for RIPK2_full |
| `pRosettaC/runs/high_affinity_degraders/MAPK14_full/Protac_params.txt` | Created | pRosettaC params for MAPK14_full |
| `docs/05_bugs_and_fixes.md` | Updated | Bugs 21–25 (N-lobe, n_dock, score.sc, atomtypes, OpenCL/WSL2) |
| `docs/06_results_summary.md` | Updated | Phase 2 status, timeline entries |

---

## 2026-06-25 — DDR2 GROMACS MD Setup

### SJF8240 parameterization: COMPLETED
`~/protac-md/ligand_params/SJF8240_3d.acpype/SJF8240_3d_GMX.itp` generated with GAFF2 + AM1-BCC.
Run via `drMD` conda env (has AmberTools tleap); AM1-BCC took 55 minutes via sqm.

### GROMACS installation: COMPLETED (2025.4 CUDA)
- Installed in `gromacs` conda env
- **MUST use 2025.4**, not 2024.5 — 2024.5 fails on RTX 5090 (sm_120 Blackwell, only compiled to sm_90)
- 2025.4 from conda-forge: `conda install -c conda-forge "gromacs=2025.4=nompi_cuda_h39c90b0_0"`
- RTX 5090 maps PP+PME+update to GPU → **251 ns/day** for DDR2 ternary complex (231K atoms)
- WSL2 OpenCL does not work (CUDA-only passthrough in WSL2) — see Bug 25

### DDR2 rank-1 topology build: COMPLETED
- 3 starting structures split: chains A=VHL, B=ElonginB, C=ElonginC, D=DDR2
- Topology: AMBER99SB-ILDN protein + SJF8240 GAFF2 + TIP3P water + 0.15M NaCl
- Bugs fixed in `build_topology.sh`: atomtypes ordering (Bug 24), genion group 15, tc-grps
- Final system: `complex_ions.gro` — 231,535 atoms, 221 NA + 216 CL

### DDR2 rank-1 equilibration + production (COMPLETE as of 2026-06-26)
| Stage | Status | Notes |
|-------|--------|-------|
| EM (50000 steps, tol 1000) | ✓ DONE | 2219 steps, Epot=-3.808×10⁶ kJ/mol, CPU run |
| NVT (200 ps, 300 K) | ✓ DONE | 251 ns/day on RTX 5090, ~60 sec |
| NPT (1 ns, 300 K, 1 bar) | ✓ DONE | 178 ns/day, ~8 min |
| Production (100 ns) | ✓ **DONE** | 206 ns/day, ~10 h; `prod.xtc` 830 MB |

### tc-grps correction
All 3 .mdp files (nvt, npt, prod) fixed from `Protein_PROTAC Water_and_ions` to:
```
tc-grps = Protein  Other  Water_and_ions
```
"Other" = PROTAC (SJF8240_3d); these are the correct GROMACS default group names.

### Files changed this session
| File | Action |
|------|--------|
| `~/protac-md/scripts/setup/build_topology.sh` | Fixed: atomtypes ordering, genion group 15 |
| `~/protac-md/scripts/setup/run_equilibration.sh` | Fixed: removed -gpu_id 0 -nb gpu from EM |
| `~/protac-md/scripts/mdp/nvt.mdp` | Fixed: tc-grps to Protein Other Water_and_ions |
| `~/protac-md/scripts/mdp/npt.mdp` | Fixed: tc-grps to Protein Other Water_and_ions |
| `~/protac-md/scripts/mdp/prod.mdp` | Fixed: tc-grps to Protein Other Water_and_ions |
| `docs/05_bugs_and_fixes.md` | Added Bug 24 (atomtypes) and Bug 25 (WSL2 OpenCL) |

---

## 2026-06-29 to 30 — Phase 2 Finalisation

### ernie runs: all complete
All three POIs finished 10,000/10,000 local docking poses. pRosettaC's own clustering step
failed for all three (MET: working-dir bug; RIPK2_full/MAPK14_full: bad SDF path) — but
combined PDBs and score.sc are intact and that is all our analysis needs.

### Analysis results
`find_feasible_poses.py` + `recluster_feasible.py` run on all 5 degrader POIs:

| POI | Bridging poses | Rank-1 PDB | Score (REU) | Quality |
|-----|---------------|-----------|------------|---------|
| DDR2 | 88 | `combined_131_9_0001.pdb` | -30,221 | PASS |
| MET | 92 | `combined_31_8_0001.pdb` | -32,947 | MARGINAL |
| RIPK2_full | 72 | `combined_135_11_0001.pdb` | -34,530 | PASS |
| MAPK14_full | 51 | `combined_37_41_0001.pdb` | -21,528 | PASS |
| EPHB2 | 5 | `combined_38_26_0001.pdb` | -32,727 | MARGINAL |

### Scripts cleaned up
Old RIPK2 and MAPK14 (pre-N-lobe-fix, spurious >90% feasibility) removed from:
- `scripts/analysis/find_feasible_poses.py`
- `scripts/analysis/recluster_feasible.py`

Only `RIPK2_full` and `MAPK14_full` remain as the canonical runs for these two POIs.

### Next: MD setup for MET, RIPK2_full, MAPK14_full
Use `feasible_md_candidates.csv` rank-1 PDB for each. Same pipeline as DDR2:
`split_chains.py` → `build_topology.sh` → `run_equilibration.sh` → `run_production.sh`
