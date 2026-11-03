# 10PROTAC-MD

**Project:** Proteome-wide pRosettaC + MD for 10 VHL-based PROTAC degraders  
**Goal:** Predict degradation across the druggable proteome by running pRosettaC for 10 PROTACs against all available kinase/protein structures, then applying ML models trained in PROTAC-MD.

---

## Overview

- **10 PROTACs:** TBD — select diverse VHL-based warhead chemotypes
- **Proteome targets:** All human proteins with available crystal structures (PDB) or AF2 models
- **Pipeline:** pRosettaC → bridging feature extraction → ML scoring (XGBoost/GNN from protac-md)

---

## Status

Project initialized: 2026-11-03

---

## Directory Layout

```
10PROTAC-MD/
├── analysis/          # Per-POI bridging + feature outputs
├── simulations/       # MD runs (top-scoring ternary complexes)
├── reference_structures/
├── ligand_params/     # ACPYPE/GAFF2 params per PROTAC
├── scripts/
│   ├── 01_setup/      # Warhead placement, topology prep
│   ├── 02_analysis/   # Bridging extraction, ML scoring
│   ├── 03_free_energy/
│   └── 04_gnn/
├── gnn/               # Graph models
├── free_energy/
├── logs/
└── docs/
```
