# 07 — Analysis Guide: How to Interpret Results

This document explains every output file produced by the analysis pipeline and how to
read the numbers.

---

## Output files per POI

After running `run_poi_full_analysis.py <POI>`, two files appear in the run directory:

### cluster_features.csv

One row per cluster. Key columns:

| Column | Meaning | Good range |
|--------|---------|------------|
| `cluster` | Cluster ID (cluster1 = largest) | — |
| `n_poses` | Number of docking solutions in this cluster | > 10 = well-populated |
| `score_orig_REU` | Rosetta energy (shifted — see below) | More negative = better |
| `interface_contacts_5A` | Number of atom pairs < 5 Å across POI–VHL interface | > 80 for degraders |
| `interface_area_est_A2` | Estimated buried surface area (contacts × 15 Å²) | > 1200 Å² |
| `poi_iface_atoms` | POI atoms contributing to interface | — |
| `vhl_iface_atoms` | VHL atoms contributing to interface | — |
| `polar_hbond_pairs` | Atom pairs where both atoms are N or O (H-bond capable) | > 2 preferred |
| `protac_disconnected` | Is the PROTAC linker broken? | Must be 0 |
| `protac_overlaps` | Severe PROTAC–protein steric overlaps | 0 ideal; < 5 acceptable |
| `clash_interface` | Hard clashes at POI–VHL interface (< 1.5 Å) | Must be 0 |
| `clash_protac_protein` | Hard clashes between PROTAC and protein | Must be 0 |
| `n_lys_within_25A_vhl` | POI lysines within 25 Å of VHL surface | Broader context |
| `n_lys_e2_optimal` | POI lysines 35–55 Å from modelled E2 active site | > 5 for degraders |
| `best_e2_lys_top3` | Top 3 E2-optimal lysines as "K{resid}({dist}A)" | — |
| `md_candidate` | Is this cluster selected for MD? (1=primary, 2=secondary, blank=no) | — |

**Score note:** `score_orig_REU` values are shifted by `score.sc.shift` relative to the
raw PyRosetta output. They are always negative (required by pRosettaC clustering). Relative
differences between clusters are meaningful; absolute values are not directly comparable
across POIs. To get the true Rosetta score: `true_score = score_orig_REU + shift_value`
where `shift_value` is in `Patchdock_Results/score.sc.shift`.

---

### premd_descriptors.csv

One row per MD candidate cluster. Key columns:

| Column | Meaning |
|--------|---------|
| `md_verdict` | READY / WARNING / REJECT — is this cluster suitable for MD? |
| `vhl_pct` | VHL structure completeness (%) — missing residues, gaps |
| `poi_pct` | POI structure completeness (%) — missing residues |
| `protac_n_atoms` | Should be 79 for SJF8240 |
| `protac_e2e_A` | End-to-end distance of PROTAC in Å (should be ~30–45 Å) |
| `pt0_atoms` / `pt1_atoms` | VHL / kinase warhead atom counts (should be 30 / 37) |
| `clash_iface_severe_15A` | Distances < 1.5 Å at interface — hard clashes |
| `clash_iface_bad_20A` | Distances < 2.0 Å — likely clashes |
| `clash_iface_warn_25A` | Distances < 2.5 Å — close contacts |
| `interface_contacts_5A` | Same as cluster_features (sanity check) |
| `closest_lys_A` | Distance to the nearest POI lysine from VHL surface (Å) |
| `n_lys_e2_optimal` | Same as cluster_features |
| `best_e2_A` | Distance of the best (closest optimal) E2 lysine (Å) |
| `needs_minimisation` | Should run energy minimisation before MD? |
| `issues` | Hard stop issues (empty = none) |
| `warnings` | Soft warnings (e.g., "POI 94.9% complete") |

#### md_verdict meanings
- **READY**: No hard clashes, PROTAC intact, structure complete enough → go to MD
- **WARNING**: Minor issues (a few close contacts, slight incompleteness) → can do MD but
  inspect first
- **REJECT**: Hard clashes at interface, disconnected PROTAC, or major structural gaps →
  do not use for MD; try the next cluster

---

## Cross-POI table: ml_features_best.csv

One row per POI. This is the input to the ML classifier.

Column naming convention:
- `best_*` = value from the best cluster (md_candidate = 1, or largest cluster if none)
- `mean_*` = mean across all clusters
- `max_*` = maximum across all clusters
- `md1_*` = value from the primary MD candidate (cluster1 or first md_candidate)

Key discriminative columns (used in the classifier):

| Column | Interpretation |
|--------|---------------|
| `best_interface_contacts_5A` | Quality of the best ternary complex geometry |
| `best_n_lys_e2_optimal` | How many accessible POI lysines can be ubiquitinated |
| `best_polar_hbond_pairs` | H-bonding at the best interface |
| `md1_best_e2_A` | Distance of the closest optimal E2 lysine — lower is better |
| `n_clusters` | Total number of productive ternary complex geometries found |

---

## Running the ML classifier

```bash
python scripts/analysis/degradation_classifier.py
```

Output:
1. Feature summary table (mean ± std per feature for degraders vs non-degraders)
2. Logistic regression LOO accuracy + confusion matrix + feature coefficients
3. Decision tree structure (depth 2) + feature importances
4. `analysis/degradation_prediction/classification_summary.csv`

**With high-affinity POIs only (most rigorous comparison):**
```bash
python scripts/analysis/degradation_classifier.py --high-affinity-only
```
This uses only the 10 high-affinity POIs (5 degraders + 5 non-degraders), removing the
3 low-affinity degraders that have different affinity confounders.

---

## How to read a ternary complex structure in PyMOL

```python
# Load a representative structure
load pRosettaC/runs/high_affinity_degraders/MET/Results/cluster2/combined_21_38_0001.pdb

# Colour by chain
util.cbc

# Show PROTAC as sticks
show sticks, resn PTC

# Show interface contacts
select interface, (chain A and polymer) within 5 of (chain C and polymer)
# (adjust chain IDs based on your structure)

# Check E2-lysine distance
# The VHL E2 binding site is approximately at the ElonginB/C surface
# K1232 in MET cluster2 is the primary candidate (37.4 Å from E2)
```

---

## Sanity checks to run on every completed POI

1. **n_clusters ≥ 2**: If only 1 cluster, the conformational sampling may be insufficient
2. **protac_disconnected = 0**: If disconnected, the PROTAC linker was not successfully built
3. **clash_interface = 0**: No hard clashes at the protein-protein interface
4. **vhl_pct > 85% and poi_pct > 85%**: Structure coverage is sufficient for analysis
5. **n_lys_e2_optimal**: Expected ≥ 5 for degraders; 0–3 for non-degraders
6. **protac_e2e_A ~30–45 Å**: The PROTAC should have a reasonable extended conformation

---

## Expected results summary once all POIs complete

Based on the literature and MET results, predicted discriminating thresholds:

| Feature | Degraders (expected) | Non-degraders (expected) |
|---------|---------------------|------------------------|
| best_interface_contacts_5A | > 80 | < 60 |
| best_n_lys_e2_optimal | > 5 | < 3 |
| md1_best_e2_A | < 45 Å | > 50 Å OR no optimal Lys |
| best_polar_hbond_pairs | > 2 | ≤ 1 |

These are hypotheses to be tested — actual results may differ.

---

## Frequently asked questions

**Why are Rosetta scores positive in some raw files?**
pRosettaC's constraint_generation uses PackRotamersMover which scores in absolute REU,
not relative to a reference. The scores are shifted to negative values by
`resume_pipeline.py::shift_scores_if_needed()`. See `docs/05_bugs_and_fixes.md`.

**What does it mean if a POI has 0 clusters?**
Clustering failed. Usually because: (a) all scores were positive and score.sc wasn't
shifted (fixed by rerunning resume_pipeline.py), or (b) the PROTAC couldn't form a linker
for any docking pose (the two warheads were too far apart or geometrically incompatible).

**Why does n_lys_within_25A_vhl differ from n_lys_e2_optimal?**
`n_lys_within_25A_vhl` counts lysines close to VHL generally (broad filter);
`n_lys_e2_optimal` specifically counts those in the 35–55 Å range from the E2 active
cysteine. A lysine can be near VHL but out of E2 reach, or far from VHL but in E2 range.

**What is the E2 modelled position?**
The E2 enzyme (UBE2D3/UbcH5c) is not explicitly present in the pRosettaC structures.
Its active cysteine position is modelled from published crystal structures of the
VHL-E2 complex. The approximate E2 cysteine position is estimated from the ElonginB/C
surface geometry in the 4W9H structure.
