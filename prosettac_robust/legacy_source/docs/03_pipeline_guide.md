# 03 — Pipeline Reproduction Guide

This document describes how to reproduce the entire analysis from scratch on a fresh machine.

---

## Prerequisites

### Hardware
- Linux (tested on Ubuntu 22.04 / WSL2 on Windows)
- ≥ 16 GB RAM (32+ GB recommended for PyRosetta with large proteins)
- ≥ 24 CPU cores recommended; minimum 4

### Software
- Miniconda3 or Anaconda
- PyRosetta (academic license — free; register at https://www.pyrosetta.org)
- PatchDock binary (included in pRosettaC distribution)

---

## Step 0: Environment setup

```bash
# Install Miniconda if not present
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Create conda environment
conda create -n protac python=3.10
conda activate protac

# Install packages
conda install -c conda-forge rdkit biopython numpy pandas scipy -y
pip install pyrosetta  # or use the wheel from PyRosetta website

# Verify
python -c "import pyrosetta; print('PyRosetta OK')"
python -c "from rdkit import Chem; print('RDKit OK')"
python -c "from Bio.PDB import Superimposer; print('BioPython OK')"
```

---

## Step 1: Set up directory structure

```bash
# Create project
cd ~
mkdir "PROTAC degradation Prediction"
cd "PROTAC degradation Prediction"

# Create the symlink to avoid space-in-path issues with PyRosetta
ln -s ~/PROTAC\ degradation\ Prediction ~/protac

# Clone pRosettaC (the Ciulli lab's original tool)
cd ~
git clone https://github.com/gilbertlabucsd/pRosettaC.git PRosettaC
# pRosettaC's main.py should be importable from Python
export PYTHONPATH="$HOME/PRosettaC:$PYTHONPATH"
```

---

## Step 2: Prepare input structures

### 2a. VHL structure

1. Download `4W9H.pdb` from RCSB PDB
2. Extract chains A, B, C (VHL, ElonginB, ElonginC); remove HETATM, water
3. Save as `data/structures/VHL/4W9H_VBC_clean.pdb`

The VHL warhead SDF (`SJF8240_vhll_placed.sdf`) was placed by mapping SJF8240's VHL
warhead SMILES onto the 3JF ligand pose in 4W9H using MCS. This file is already present
in the repository.

### 2b. POI kinase structures

For each POI, download the PDB listed in `docs/02_structure_preparation.md`, extract
chain A (ATOM records only), remove ligands/water, save as `<POI>_protein_clean.pdb`.

Alternatively, `scripts/structure_prep/download_structures.py` automates most downloads.

### 2c. Place warheads

```bash
cd ~/protac
python scripts/structure_prep/place_warhead.py
```

This generates `<POI>_warhead_placed.sdf` and `<POI>_warhead_placed_H.sdf` for every POI.
Check the RMSD and anchor-atom distance printed to stdout:

- RMSD > 3.0 Å → flag for manual inspection (kinase fold divergence from MET)
- anchor-atom to nearest Cα > 10 Å → warhead is outside protein — placement failed

---

## Step 3: Prepare pRosettaC run directories

For each POI, create a run directory and `Protac_params.txt`:

```
pRosettaC/runs/<category>/<POI>/Protac_params.txt
```

Template for a high-affinity degrader (MET example):

```
Structures: /home/<user>/protac/data/structures/VHL/4W9H_VBC_clean.pdb /home/<user>/protac/data/structures/POI/high_affinity_degraders/MET/MET_protein_clean.pdb
Chains: ABC A
Heads: /home/<user>/protac/data/structures/VHL/SJF8240_vhll_placed.sdf /home/<user>/protac/data/structures/POI/high_affinity_degraders/MET/MET_warhead_placed.sdf
Anchor atoms: 6 1
Protac: /home/<user>/protac/data/structures/PROTAC/SJF8240/SJF8240.smi
Full: True
ClusterName: LOCAL
```

**Key fields:**
- `Structures`: absolute paths to VHL-VBC PDB and POI PDB, space-separated
- `Chains`: chain string for the combined system (ABC = VHL+ElonginB+ElonginC) followed
  by space and the POI chain identifier within that system (A = kinase chain)
- `Heads`: absolute paths to VHL warhead SDF and kinase warhead SDF
- `Anchor atoms`: atom indices (1-based) in PT0 and PT1 where the linker attaches.
  For SJF8240: PT0 anchor = atom 6 (N6, the amide nitrogen connecting to linker),
  PT1 anchor = atom 1 (C1, the first carbon of the foretinib scaffold connecting to linker)
- `Protac`: SMILES file for the full PROTAC molecule
- `Full: True`: use full PatchDock scoring (not fast approximation)
- `ClusterName: LOCAL`: use local clustering (as opposed to remote server)

---

## Step 4: Run pRosettaC for each POI

### Option A: Single POI (manual)

```bash
cd ~/PROTAC\ degradation\ Prediction/pRosettaC/runs/high_affinity_degraders/MET
python ~/protac/scripts/rosetta_utils/run_local_prosettac.py Protac_params.txt \
  --n_dock 200 --n_workers 18
```

### Option B: Sequential batch (all POIs, this machine)

```bash
cd ~/PROTAC\ degradation\ Prediction
N_DOCK=200 N_WORKERS=18 N_WORKERS_RESUME=18 \
  nohup bash scripts/rosetta_utils/run_batch_rerun.sh \
  > pRosettaC/runs/batch_rerun.log 2>&1 &

# Monitor
tail -f pRosettaC/runs/batch_rerun.log
```

### Option C: Resume from constraint_generation (if docking already done)

If PatchDock docking is already complete but constraint_generation/clustering hasn't run:

```bash
cd ~/PROTAC\ degradation\ Prediction/pRosettaC/runs/high_affinity_degraders/DDR2
N_WORKERS_RESUME=18 python ~/protac/scripts/rosetta_utils/resume_pipeline.py
```

### Option D: HPC cluster (parallel, recommended for all 11 remaining POIs)

See `docs/08_hpc_guide.md` (to be written once cluster setup is confirmed).

---

## Step 5: What the pipeline does internally

pRosettaC runs in 4 sequential stages:

### Stage 1: PatchDock global docking
- Takes the VHL-warhead half of the system and docks it against the POI kinase domain
- Generates N_DOCK = 200 global docking poses (rigid body, shape complementarity)
- Output: `Patchdock_Results/pd.*_docking_????.pdb` (200 files × 50 local poses = 10,000)

### Stage 2: Local Rosetta docking (FastRelax)
- For each of the 200 global poses, runs N=50 local Rosetta relaxations
- Optimises side-chain rotamers and backbone near the interface
- Output: 10,000 `pd.*_docking_????.pdb` structures in `Patchdock_Results/`

### Stage 3: Constraint generation (PROTAC linker modelling)
- For each of the 10,000 docked poses, attempts to build a physically reasonable PROTAC
  linker that connects the VHL warhead (already placed in VHL) to the kinase warhead
  (already placed in the kinase) using the SMILES-defined linker topology
- Uses RDKit conformer generation with distance constraints on the anchor atoms
- Runs PyRosetta PackRotamersMover to relax the protein-PROTAC interface
- Output: `combined_<pd>_<local>_0001.pdb` for each successful pose
- Typically 50–70% of poses succeed (others rejected: no valid linker conformation found)

### Stage 4: Clustering
- Clusters the successful combined PDB structures by interface Cα RMSD
- Top clusters become the final output in `Results/`
- Each cluster representative = `Results/cluster<N>/combined_<pd>_<local>_0001.pdb`

---

## Step 6: Analyse results

### Per-POI analysis

```bash
python scripts/analysis/run_poi_full_analysis.py MET
```

This generates two files in the run directory:
- `cluster_features.csv` — structural features per cluster (interface contacts, lysine distances, etc.)
- `premd_descriptors.csv` — quality checks and MD-readiness assessment for top clusters

### Cross-POI aggregation

```bash
python scripts/analysis/aggregate_features.py
```

Reads all `cluster_features.csv` and `premd_descriptors.csv` files; writes:
- `analysis/degradation_prediction/ml_features_best.csv` — one row per POI, best cluster

### ML classification

```bash
python scripts/analysis/degradation_classifier.py
# or
python scripts/analysis/degradation_classifier.py --high-affinity-only
```

---

## Step 7: Expected outputs per POI

After a complete successful run, each POI run directory should contain:

```
<POI>/
├── Protac_params.txt          <- input parameters
├── Init0.pdb                  <- VHL half-complex
├── Init1.pdb                  <- POI kinase domain
├── PT0.params                 <- Rosetta params for VHL warhead
├── PT1.params                 <- Rosetta params for kinase warhead
├── Patchdock_Results/
│   ├── pd.*_docking_????.pdb  <- 10,000 locally docked poses
│   ├── combined_*_0001.pdb    <- ~5,000-7,000 successful linker poses
│   ├── local.fasc             <- Rosetta score file
│   └── score.sc               <- constraint_generation scores
├── Results/
│   ├── cluster1/              <- top cluster (most poses)
│   │   ├── combined_*.pdb     <- all members of cluster1
│   │   └── ...
│   ├── cluster2/ ... cluster6/
│   └── Init.pdb
├── cluster_features.csv       <- per-cluster structural features
├── cluster_features.json
├── premd_descriptors.csv      <- MD-readiness assessment
└── premd_descriptors.json
```

---

## Troubleshooting

**"No conformations were generated" in constraint_generation**
- Normal for some poses — linker cannot span the distance between anchor atoms
- If > 90% of poses show this, the anchor atom indices may be wrong or the warhead
  placement has put the anchors in an unreachable geometry → re-check warhead.pdb

**score.sc all positive values → clustering fails**
- `resume_pipeline.py` handles this via `shift_scores_if_needed()`, which shifts all
  scores by `-(max_score + 100)` so they become negative (pRosettaC clustering requires
  negative scores)
- Original scores backed up as `score.sc.orig`; shift amount in `score.sc.shift`

**Results/ created but empty**
- Clustering ran but found no structures passing the score < 0 filter
- First check if `score.sc.shift` exists — if not, scores may still be positive
- Run `resume_pipeline.py` again; it will detect and fix the score shift

**PyRosetta init error with "cannot find extra_res_fa file"**
- Caused by space in path `PROTAC degradation Prediction`
- Verify the `~/protac` symlink exists: `ls -la ~/protac`
- All paths in params files must use `~/protac/...` not `~/PROTAC degradation Prediction/...`

**Two processes both running for the same POI**
- Run `ps aux | grep run_local_prosettac` to identify PIDs
- Kill all but one, or add an in-progress check to the batch script
- `run_batch_rerun.sh` includes a guard: skips if `Patchdock_Results/` is non-empty
