# PROTAC Degradation Prediction — SJF8240 / VHL Campaign

**Project started:** 2026-05-25  
**Lead:** Ayush Regmi Bagale  
**Objective:** Use pRosettaC ternary complex modelling to identify structural features that
discriminate kinases degraded by SJF8240 from those that are not.

---

## Quick Navigation

| Document | Contents |
|----------|----------|
| [docs/01_scientific_rationale.md](docs/01_scientific_rationale.md) | What SJF8240 is, why this approach, biological background |
| [docs/02_structure_preparation.md](docs/02_structure_preparation.md) | Every PDB used, why it was chosen, quality, what was cleaned |
| [docs/03_pipeline_guide.md](docs/03_pipeline_guide.md) | Complete step-by-step reproduction guide from scratch |
| [docs/04_parameter_justification.md](docs/04_parameter_justification.md) | Why each pRosettaC parameter was chosen |
| [docs/05_bugs_and_fixes.md](docs/05_bugs_and_fixes.md) | All bugs found, root causes, and exact fixes applied |
| [docs/06_results_summary.md](docs/06_results_summary.md) | Current results: MET (complete), DDR2 (pipeline check) |
| [docs/07_analysis_guide.md](docs/07_analysis_guide.md) | How to interpret cluster_features.csv, premd_descriptors.csv, ML output |

---

## PROTAC: SJF8240

| Property | Value |
|----------|-------|
| E3 ligase recruited | VHL (Von Hippel-Lindau) |
| Kinase warhead | Foretinib / XL880 / GSK1363089 pharmacophore |
| Total heavy atoms | 79 (PT0: 30 VHL warhead + PT1: 37 kinase warhead + linker: 12) |
| SMILES | `data/structures/PROTAC/SJF8240/SJF8240.smi` |
| Source | Huang et al. — pan-kinase degrader study |

---

## POI Panel (13 kinases)

### High-affinity degraders (Priority 1)
| POI | PDB used | Resolution | Ligand | Run status |
|-----|----------|-----------|--------|-----------|
| MET | 3LQ8 | 2.02Å | XL880/foretinib (88Z) | ✅ Complete — 6 clusters |
| DDR2 | **7AYM** | 2.12Å | IBZ3/4VD (DDR2-selective) | ⚠️ Pipeline check done; needs re-run with 7AYM |
| RIPK2 | 4C8B | 2.75Å | Ponatinib | ⏸ Queued |
| EPHB2 | **3ZFM** | 2.27Å | Apo | ⏸ Queued |
| MAPK14 | 3GCP | 2.25Å | SB203580 | ⏸ Queued |

### High-affinity non-degraders (Priority 2)
| POI | PDB used | Resolution | Ligand | Run status |
|-----|----------|-----------|--------|-----------|
| AXL | 5U6B | 2.84Å | Macrocyclic inhibitor | ⏸ Queued |
| SLK | 2J51 | 2.10Å | DKI | ⏸ Queued |
| ABL1 | 2HYY | 2.40Å | Imatinib (type II) | ⏸ Queued |
| EPHA2 | 5I9Y | 1.23Å | Dasatinib | ⏸ Queued |
| MAP4K5 | AF2 Q9Y4K4 v6 | — | AF2 model | ⏸ Queued |

### Low-affinity degraders (Priority 3)
| POI | PDB used | Resolution | Ligand | Run status |
|-----|----------|-----------|--------|-----------|
| TNIK | 5D7A | — | — | ⏸ Queued |
| PIP4K2C | 2GK9 | — | — | ⏸ Queued |
| CDK17 | AF2 O94921 v6 | — | AF2 model | ⏸ Queued |

> **Bold PDB IDs** = corrected during this campaign (see `docs/05_bugs_and_fixes.md`)

---

## Directory Structure

```
PROTAC degradation Prediction/
├── README.md                          <- this file
├── docs/                              <- all documentation
│   ├── 01_scientific_rationale.md
│   ├── 02_structure_preparation.md
│   ├── 03_pipeline_guide.md
│   ├── 04_parameter_justification.md
│   ├── 05_bugs_and_fixes.md
│   ├── 06_results_summary.md
│   └── 07_analysis_guide.md
├── data/
│   ├── structures/
│   │   ├── VHL/                       <- VHL:ElonginB:ElonginC + warhead SDF
│   │   ├── POI/                       <- kinase PDBs and placed warhead SDFs
│   │   │   ├── high_affinity_degraders/
│   │   │   ├── high_affinity_no_degradation/
│   │   │   └── low_affinity_degraders/
│   │   └── PROTAC/SJF8240/            <- SMILES, SDF, params
│   └── ligand_params/                 <- Rosetta .params files
├── pRosettaC/
│   └── runs/                          <- one sub-directory per POI
│       ├── high_affinity_degraders/MET/   <- full run + Results/ + analysis CSVs
│       ├── high_affinity_degraders/DDR2/
│       ├── high_affinity_no_degradation/
│       └── low_affinity_degraders/
├── scripts/
│   ├── structure_prep/
│   │   ├── place_warhead.py           <- MCS + Cα superposition warhead placement
│   │   └── download_structures.py
│   ├── rosetta_utils/
│   │   ├── run_local_prosettac.py     <- main pRosettaC wrapper
│   │   ├── resume_pipeline.py         <- resume from constraint_generation step
│   │   ├── run_batch_rerun.sh         <- sequential batch runner (all 12 POIs)
│   │   ├── watch_and_analyse.sh       <- auto-analysis watcher daemon
│   │   ├── constraint_generation.py   <- PROTAC conformer + Rosetta relaxation
│   │   └── pyrosetta_adapter.py       <- PyRosetta path/patch fixes for this env
│   └── analysis/
│       ├── run_poi_full_analysis.py   <- per-POI cluster feature extraction
│       ├── aggregate_features.py      <- cross-POI ML feature table builder
│       └── degradation_classifier.py  <- LOO logistic regression + decision tree
└── analysis/
    └── degradation_prediction/
        └── ml_features_best.csv       <- final ML feature table (one row per POI)
```

---

## Symlink note

The project directory name contains a space (`PROTAC degradation Prediction`) which breaks
PyRosetta's `-extra_res_fa` flag parser (it splits on whitespace). A symlink fixes this:

```bash
ln -s ~/PROTAC\ degradation\ Prediction ~/protac
```

All scripts that invoke PyRosetta use `~/protac/...` paths internally.

---

## How to re-run everything from scratch

See `docs/03_pipeline_guide.md` for the complete reproduction guide.

**Quick start** (environment already set up):

```bash
# 1. Place warheads for all POIs
cd ~/protac
python scripts/structure_prep/place_warhead.py

# 2. Run full sequential batch (18 workers)
cd ~/PROTAC\ degradation\ Prediction
N_DOCK=200 N_WORKERS=18 N_WORKERS_RESUME=18 \
  nohup bash scripts/rosetta_utils/run_batch_rerun.sh \
  > pRosettaC/runs/batch_rerun.log 2>&1 &

# 3. Analyse results as POIs complete
python scripts/analysis/run_poi_full_analysis.py MET   # or any POI name

# 4. Aggregate all POIs into ML feature table
python scripts/analysis/aggregate_features.py

# 5. Run classifier
python scripts/analysis/degradation_classifier.py
```
