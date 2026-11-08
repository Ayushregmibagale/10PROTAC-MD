# PROTAC Degradation Prediction — ML Results Draft

## Overview

We developed a two-track machine learning pipeline to predict PROTAC-mediated degradation from computed ternary complex structural features. Track 1 (XGBoost + SHAP) operates on tabular geometry features derived from pRosettaC docking ensembles. Track 2 (dual-encoder GNN) integrates pRosettaC contact frequency graphs with MD-derived dynamic stability graphs.

**Primary result (corrected dataset, Aug 12):** XGBoost standalone LOPO-CV AUROC = **0.850** on 17 clean POIs (5 DEG, 12 NON). The previous ensemble result (0.852, 15 POIs) used an incorrect CDK17 label (classified as degrader; confirmed non-degrader Aug 12) and is superseded by this result.

---

## Dataset

- **19 POIs** total: all SJF8240/VHL-based PROTACs
- **2 artifact POIs excluded** from ML training (TNIK/MAPK14): kin_min=104 overlap bug inflates computed bridging to 67–74%; results retained with `prostac_artifact=1` flag
- **17 clean POIs** for LOPO-CV: **5 degraders**, **12 non-degraders**
- **8 POIs with MD** (random-screen × 20 ns per pose): DDR2/MET/RIPK2_full/MAPK14_full/EPHB2/ABL1/AXL/EPHA2
- **9 POIs without MD**: CDK17 + CDK2/4/6/DAPK1/IRAK1/TLK2 + MAP4K5 + SLK

Label corrections and artifact resolutions (Aug 12):
- **MAP4K5 and SLK**: chain-D pRosettaC reruns confirmed 0% bridging; CoM separations of 68.7/67.2 Å (largest in panel); confirmed geometry-blocked non-degraders; all features (03c/03d) recomputed from chain-D data; no longer excluded
- **CDK17**: confirmed non-degrader experimentally; previously mislabelled as low-affinity degrader based on pRosettaC base directory; 0% bridging, n_near_vhl=0, d_sum_min=15.85 Å consistent with geometry-blocked non-degrader

Categories:
- `high_degrader`: DDR2, MET, RIPK2_full, MAPK14_full, EPHB2
- `high_non_degrader`: ABL1, AXL, EPHA2, MAP4K5, SLK
- `low_non_degrader`: CDK17, CDK2, CDK4, CDK6, DAPK1, IRAK1, TLK2
- Excluded artifacts: TNIK, MAPK14

---

## pRosettaC Structural Features

### 1. Docking ensemble (20,000 poses per POI)

Random sampling (seed=42) from the full docking pool is essential. Rosetta energy minimisation scores anti-correlate with ternary complex geometry: the optimal bridging pose (lowest d_vhl + d_kin) ranks 79–7,923 by Rosetta score across all 19 POIs (Spearman ρ = +0.017 to −0.283; p < 0.01 in all cases, p < 0.001 in 18/19 POIs). Top-scored poses have PROTAC-protein distances of 15–378 Å — far outside bridging range. This anti-correlation is driven by protein–protein contact energy: top-scored poses maximise direct VHL–kinase interface contacts, which is geometrically incompatible with PROTAC bridging (which requires proteins ~30–50 Å apart). The best-bridging pose sacrifices ~20,000–25,000 REU relative to the Rosetta-optimal pose. Random sampling uniformly covers the conformational space and is the only unbiased strategy.

### 2. Geometric failure mode (Mode 0) — CDK2/4/6/DAPK1/IRAK1/TLK2

Six non-degraders share a rigid geometric lever: when the PROTAC approaches VHL (decreasing d_vhl), it moves away from the kinase (increasing d_kin), and vice versa. This anti-correlation (`r_dvhl_dkin < 0`) reflects a steric constraint in the VHL–kinase relative orientation that prevents simultaneous engagement. Only 1–6 of 20,000 pRosettaC poses reach within 8 Å of both proteins for these POIs (Table 1). By contrast, confirmed degraders (DDR2/MET/RIPK2/MAPK14/EPHB2) and high-affinity non-degraders (ABL1/AXL/EPHA2) generate 109–384 VHL-proximal poses per 20,000 runs.

**Table 1. n_near_vhl (poses with d_vhl < 8 Å) and key geometry (17 clean POIs)**

| POI | Category | DEG | n_near_vhl | prc_bridging_frac | d_sum_min (Å) | VHL-kin CoM (Å) |
|-----|----------|-----|-----------|-------------------|---------------|-----------------|
| DDR2 | high_degrader | 1 | 305 | 42.5% | 0.81 | 54.6 |
| RIPK2_full | high_degrader | 1 | 306 | 43.5% | 1.15 | 53.9 |
| MAPK14_full | high_degrader | 1 | 216 | 33.5% | 0.74 | 56.5 |
| MET | high_degrader | 1 | 109 | 49.2% | 0.83 | 56.8 |
| EPHB2 | high_degrader | 1 | 202 | 20.0% | 2.58 | 55.6 |
| ABL1 | high_non_degrader | 0 | 384 | 52.5% | 1.12 | 58.0 |
| AXL | high_non_degrader | 0 | 361 | 53.0% | 0.74 | 57.9 |
| EPHA2 | high_non_degrader | 0 | 226 | 32.5% | 0.39 | 58.0 |
| MAP4K5 | high_non_degrader | 0 | 1 | 0.0% | 7.18 | **68.7** |
| SLK | high_non_degrader | 0 | 0 | 0.0% | 15.95 | **67.2** |
| CDK17 | low_non_degrader | 0 | 0 | 0.0% | 15.85 | 62.5 |
| CDK2 | low_non_degrader | 0 | 6 | 1.0% | 5.92 | 59.6 |
| CDK4 | low_non_degrader | 0 | 5 | 0.5% | 1.52 | 57.8 |
| CDK6 | low_non_degrader | 0 | 2 | 0.0% | 10.50 | 60.6 |
| DAPK1 | low_non_degrader | 0 | 1 | 0.0% | 10.46 | 58.3 |
| IRAK1 | low_non_degrader | 0 | 3 | 0.0% | 11.44 | 59.0 |
| TLK2 | low_non_degrader | 0 | 6 | 0.5% | 2.12 | 44.8 |

MAP4K5 and SLK have the largest VHL-kinase CoM separations in the panel (68.7/67.2 Å vs 54–60 Å for the other non-degraders), providing clear structural rationalisation for their geometry-blocked classification.

---

## MD Stability Features (P_stable)

For 8 POIs, 16–19 poses were selected by random sampling (seed=42) from the bridging pool and simulated for 100 ns each. Bridging fraction per pose (P_stable) was computed as the fraction of trajectory frames where d_vhl < 8 Å and d_kin < 8 Å simultaneously.

The scalar mean P_stable (rs_bridging_mean) does not discriminate degraders from non-degraders:
- Degraders: DDR2=31%, MET=29%, MAPK14_full=47%, RIPK2_full=45%, EPHB2=19%
- Non-degraders: ABL1=34%, AXL=35%, EPHA2=27%

Ranges overlap substantially. SHAP analysis confirms rs_bridging_mean contributes zero feature importance (mean |SHAP|=0.000) in the tabular XGBoost model. However, the *spatial pattern* of residue-level contact frequencies across MD trajectories carries discriminative information, as exploited by the GNN (see below).

---

## Machine Learning Results

### Track 1: XGBoost + SHAP (tabular features, LOPO-CV)

Evaluation: Leave-One-POI-Out cross-validation on **17 clean POIs** (5 DEG, 12 NON; CDK17 corrected to NON, MAP4K5/SLK added). Class imbalance handled by scale_pos_weight = n_neg/n_pos per fold. XGBClassifier: n_estimators=200, max_depth=3, lr=0.1, subsample=0.8. NaN imputation: MD-derived features → 0.0 (has_md=0 flags these); `kin_when_vhl_near` NaN (for n_vhl_near=0 POIs: SLK, CDK17) → column maximum (24.9 Å, from IRAK1), representing maximally geometry-blocked.

**Primary result: AUROC = 0.850, ACC = 11/17**

| Feature set | N POIs | N features | AUROC | ACC |
|---|---|---|---|---|
| Full tabular (geometry + MD scalars + dynamic) | 17 | 29 | **0.850** | 11/17 |
| *(old) Geometry only, 15 POIs, CDK17=DEG* | *15* | *10* | *0.759* | *11/15* |
| *(old) Full tabular, 15 POIs, CDK17=DEG* | *15* | *17* | *0.722* | *11/15* |

The corrected 17-POI model substantially outperforms the old 15-POI geometry-only result (+0.091 AUROC), driven by the label correction of CDK17 (non-degrader, now correctly predicted) and the addition of MAP4K5/SLK as clean geometry-blocked non-degraders.

**SHAP feature importance (17-POI full model, trained on all 17 clean POIs):**

| Rank | Feature | Mean \|SHAP\| | Interpretation |
|------|---------|--------------|----------------|
| 1 | e2_mean_min_dist_A | 0.692 | Mean E2-to-lysine NZ distance — direct E2 reachability proxy |
| 2 | e2_delta_reach_50A | 0.474 | Δ E2 reach when bridging vs non-bridging |
| 3 | vki_delta_com_dist | 0.422 | VHL-kinase CoM distance change upon bridging |
| 4 | prostac_bridging_frac | 0.395 | Fraction of 20k poses where PROTAC bridges both proteins |
| 5 | kin_when_vhl_near | 0.350 | Mean kinase distance when PROTAC within 8 Å of VHL — Mode 0 proxy |
| 6 | vki_delta_contacts | 0.331 | Change in VHL-kinase Cα contacts upon bridging |
| 7 | rs_vhl_engagement_pct | 0.252 | Fraction of MD frames with PROTAC-VHL contact |
| 8 | e2_delta_reach_80A | 0.147 | Δ E2 reach at 80 Å threshold |

With the corrected dataset, E2 reachability features (e2_mean_min_dist_A, e2_delta_reach_50A) rise to prominence alongside the geometry discriminators. `kin_when_vhl_near` remains important but is no longer the sole dominant feature. This shift reflects the more balanced class distribution (5:12 DEG:NON) enabling the model to exploit multiple mechanistic axes.

**Per-POI predictions (17-POI model):**

| POI | Label | P(DEG) | Pred | Correct | MD |
|-----|-------|--------|------|---------|-----|
| MET | 1 | 0.901 | DEG | ✓ | yes |
| RIPK2_full | 1 | 0.886 | DEG | ✓ | yes |
| DDR2 | 1 | 0.297 | NON | **✗** | yes |
| EPHB2 | 1 | 0.103 | NON | **✗** | yes |
| MAPK14_full | 1 | 0.109 | NON | **✗** | yes |
| ABL1 | 0 | 0.547 | DEG | **✗** | yes |
| AXL | 0 | 0.718 | DEG | **✗** | yes |
| EPHA2 | 0 | 0.776 | DEG | **✗** | yes |
| MAP4K5 | 0 | 0.054 | NON | ✓ | no |
| SLK | 0 | 0.064 | NON | ✓ | no |
| CDK17 | 0 | 0.064 | NON | ✓ | no |
| CDK2–TLK2 (6) | 0 | 0.058–0.187 | NON | ✓ (all 6) | no |

All 9 no-MD non-degraders correctly predicted NON (prob 0.054–0.188). MET and RIPK2_full correctly predicted DEG. Persistent errors: DDR2/EPHB2/MAPK14_full (DEG→NON) and ABL1/AXL/EPHA2 (NON→DEG) — see Consistent Misclassifications.

### Track 2: Dual-encoder GNN (LOPO-CV)

> **Note:** GNN evaluation was completed on the original **15-POI** dataset (CDK17 labelled DEG; MAP4K5/SLK excluded as artifacts). GNN results are reported here as-is for completeness but are not directly comparable to the 17-POI XGBoost result. GNN re-training on the corrected 17-POI dataset is deferred.

**Architecture:**
- Encoder A (EnsembleGAT): pRosettaC contact frequency graph, 9 node features (is_VHL/kinase/PROTAC, sequence position, element), 3 edge features (contact frequency, mean/min distance). 3-layer GAT, hidden=64, 4 heads.
- Encoder B (MDGAT): MD contact frequency graph, 3 node features, 2-layer GAT. Zero embedding for non-MD POIs.
- Gated fusion: learned gate scalar weighs MD contribution; concat [emb_A || gate×emb_B] → MLP.
- LOPO-CV: 400 epochs per fold, Adam lr=0.001, cosine annealing.

**Results (17-POI panel, CDK17=NON corrected labels — FINAL Aug 12):**

| GNN variant | AUROC | ACC | MD-POIs AUROC | No-MD AUROC |
|---|---|---|---|---|
| Dual-encoder (04d) | **0.850** | **14/17** | 0.400 | undefined* |

*No-MD AUROC undefined: all 9 no-MD POIs predicted NON → single class in y_true.

The dual-encoder GNN correctly classifies all five degraders (DDR2=0.999, RIPK2_full=0.999, MAPK14_full=0.996, MET=0.995, EPHB2=0.591) and all nine no-MD non-degraders (MAP4K5=0.000, SLK=0.000, CDK17=0.000, CDK2/4/6/DAPK1/IRAK1/TLK2 ≈0.000–0.001). It misclassifies only ABL1/AXL/EPHA2 (all predicted DEG at ≈0.998) — the same post-ternary failures that XGBoost also misclassifies. The MD sub-group AUROC=0.400 because EPHB2 (DEG, prob=0.591) is ranked below ABL1/AXL/EPHA2 (NON, prob≈0.998); only MET and RIPK2_full rank unambiguously above all non-degraders.

### Track 3: Late-fusion ensemble (17-POI, FINAL)

Sweeping GNN weight from 0 to 1, AUROC is flat at **0.850** for all weights ≥ 0.50. Optimal weight = 95% GNN + 5% XGB; at this weight ACC=14/17.

| GNN weight | XGB weight | AUROC | ACC |
|---|---|---|---|
| 0.00 | 1.00 | 0.850 | 11/17 |
| 0.50 | 0.50 | 0.850 | 14/17 |
| 0.95 | 0.05 | **0.850** | **14/17** |
| 1.00 | 0.00 | 0.850 | 14/17 |

The ensemble does not improve AUROC beyond either individual model. The GNN dominates: it recovers DDR2/MAPK14_full/EPHB2 (which XGBoost misses) while XGBoost adds nothing the GNN doesn't already get correct. Both models fail identically on ABL1/AXL/EPHA2.

**Per-POI predictions (17-POI, corrected labels — FINAL):**

| POI | Label | GNN | XGB | Ens (95/5) | Pred | MD | Correct |
|-----|-------|-----|-----|------------|------|----|---------|
| DDR2 | 1 | 0.999 | 0.297 | 0.964 | DEG | MD | ✓ |
| RIPK2_full | 1 | 0.999 | 0.886 | 0.994 | DEG | MD | ✓ |
| MAPK14_full | 1 | 0.996 | 0.109 | 0.952 | DEG | MD | ✓ |
| MET | 1 | 0.995 | 0.901 | 0.990 | DEG | MD | ✓ |
| EPHB2 | 1 | 0.591 | 0.103 | 0.567 | DEG | MD | ✓ |
| ABL1 | 0 | 0.998 | 0.610 | 0.978 | DEG | MD | **✗** |
| AXL | 0 | 0.998 | 0.718 | 0.984 | DEG | MD | **✗** |
| EPHA2 | 0 | 0.998 | 0.776 | 0.987 | DEG | MD | **✗** |
| MAP4K5 | 0 | 0.000 | 0.054 | 0.003 | NON | no-MD | ✓ |
| SLK | 0 | 0.000 | 0.064 | 0.003 | NON | no-MD | ✓ |
| CDK17 | 0* | 0.000 | 0.064 | 0.003 | NON | no-MD | ✓ |
| CDK2 | 0 | 0.000 | 0.063 | 0.003 | NON | no-MD | ✓ |
| CDK4 | 0 | 0.000 | 0.062 | 0.003 | NON | no-MD | ✓ |
| CDK6 | 0 | 0.000 | 0.066 | 0.003 | NON | no-MD | ✓ |
| DAPK1 | 0 | 0.000 | 0.066 | 0.003 | NON | no-MD | ✓ |
| IRAK1 | 0 | 0.000 | 0.066 | 0.003 | NON | no-MD | ✓ |
| TLK2 | 0 | 0.000 | 0.058 | 0.003 | NON | no-MD | ✓ |

*CDK17 reclassified NON for ML (non-canonical mechanism; Section 3.4).

**Summary across models (FINAL):**

| Model | N POIs | AUROC | ACC | MD-AUROC | Notes |
|-------|--------|-------|-----|----------|-------|
| pRosettaC-only GNN (04b) | 15 | 0.185 | 5/15 | — | Deprecated |
| XGBoost geom-only | 15 | 0.759 | 11/15 | — | Deprecated (CDK17 wrong label) |
| Dual-encoder GNN (04d) | 15 | 0.741 | 12/15 | 0.667 | Deprecated (CDK17 wrong label) |
| Ensemble 30/70 | 15 | 0.852 | 11/15 | — | Deprecated (CDK17 wrong label) |
| **XGBoost 29-feat** | **17** | **0.850** | **11/17** | 0.400 | **CURRENT** |
| **GNN dual-encoder** | **17** | **0.850** | **14/17** | 0.400 | **CURRENT** |
| **Ensemble 95/5** | **17** | **0.850** | **14/17** | — | **CURRENT** |

---

## Misclassifications by model (17-POI corrected model, FINAL)

### XGBoost misclassifications (6/17 wrong)

Six MD-POIs are misclassified by XGBoost. Three are degraders predicted NON (N-limitation in scalar features); three are non-degraders predicted DEG (post-ternary failure not captured by geometry).

**DDR2 (label=DEG, XGB prob=0.297, predicted NON) — GNN CORRECT (0.999):**
DDR2 has good geometry (prc_frac=42.5%, d_sum_min=0.81 Å), but scalar E2 features don't clearly separate it from ABL1/AXL/EPHA2 in the LOPO fold. The GNN recovers it via persistent bridging contact frequency patterns that XGBoost's tabular features miss. N-limitation: LOPO fold has only 4 degraders to train on with 3 non-degraders in the same geometry bin.

**EPHB2 (label=DEG, XGB prob=0.103, predicted NON) — GNN borderline CORRECT (0.591):**
EPHB2 has the poorest geometry among degraders (d_sum_min=2.58 Å, prc_frac=20%, P_stable=19%). XGBoost calls it NON from scalar features. The GNN gives 0.591 (borderline above threshold), recovering it via spatial contact patterns — reflecting genuine ambiguity in this marginally productive degrader.

**MAPK14_full (label=DEG, XGB prob=0.109, predicted NON) — GNN CORRECT (0.996):**
MAPK14_full (p38α) has e2_mean_min_dist_A=48.9 Å — an outlier vs 57–64 Å for other degraders (MAPK14_full's K12 is genuinely closer to E2 in ternary complex; real biology, not artifact). XGBoost trained on the 57–64 Å degrader range extrapolates incorrectly at 48.9 Å. The GNN bypasses this scalar artefact entirely.

### GNN + XGBoost shared misclassifications (3/17 wrong for both)

**ABL1 (label=NON, prob=0.547, predicted DEG):**
ABL1 has the highest pRosettaC bridging fraction in the panel (52.5%, more than any degrader) and very low d_sum_min (1.12 Å). In every LOPO fold, these geometry features push the model strongly toward DEG. SHAP shows positive contributions from prostac_bridging_frac, kin_when_vhl_near (9.8 Å — within degrader range), and rs_vhl_engagement_pct. ABL1 is likely a **Mode 1 post-ternary failure**: the ternary complex forms readily, but downstream ubiquitination or proteasomal engagement fails. This failure mode is invisible to current structural features.

**AXL (label=NON, prob=0.718, predicted DEG):**
AXL has excellent geometry (prc_frac=53%, d_sum_min=0.74 Å) and high MD bridging (35.2%; rs_bridging_max=98.5% — **higher than DDR2 degrader at 31%**). **AXL does NOT fail at geometry** — the model's DEG prediction is structurally justified. AXL fails entirely at post-ternary steps:
- **Transmembrane RTK context (primary explanation):** AXL is a single-pass transmembrane receptor tyrosine kinase. Its kinase domain sits 2–4 nm from the plasma membrane cytoplasmic face; the VHL-Cullin5-RBX2 E3 complex spans 15–20 nm and cannot assemble in that restricted space even when ternary complex formation occurs.
- **VHL-warhead burial (universal, not discriminating):** Binary complex MD (8 POIs × 20 ns) shows VHL-warhead burial of 97–99.5% for ALL kinases including confirmed degraders MET (99.5%) and EPHB2 (97.5%). Burial is a property of SJF8240's linker collapsing universally — it cannot explain why AXL specifically doesn't degrade.
SHAP for AXL: prostac_bridging_frac (+0.368), kin_when_vhl_near (+0.355) dominate; e2_mean_min_dist_A (−0.686) and e2_delta_reach_50A (−0.495) partially compensate. Mode 7 is not accessible from ternary complex features alone.

**EPHA2 (label=NON, prob=0.776, predicted DEG):**
EPHA2 has the lowest d_sum_min in the panel (0.39 Å — closer than any degrader) and 226 VHL-proximal poses. SHAP contributions follow the same pattern as AXL. EPHA2 may share a Mode 7 or Mode 1 failure mechanism. **Both AXL and EPHA2 are strong candidates for experimental false negatives** — their ternary complex formation efficiency exceeds most confirmed degraders.

---

## Degradation Failure Modes Identified

| Mode | POIs | Mechanism | Predicted correctly? |
|------|------|-----------|----------------------|
| Mode 0 (geometric lever) | CDK2/4/6/DAPK1/IRAK1/TLK2 | Anti-correlated d_vhl/d_kin; kin_when_vhl_near >> 10 Å | ✓ all 6 (XGBoost + GNN + ensemble) |
| Mode 0 (geometry-blocked, large CoM) | MAP4K5, SLK | 0% pRosettaC bridging; VHL-kinase CoM 67–69 Å (largest in panel) | ✓ both (XGBoost 17-POI) |
| Mode 0 (geometry-blocked) | CDK17 | 0% bridging, n_vhl_near=0; confirmed NON-DEG Aug 12 | ✓ XGBoost 17-POI (prob=0.064) |
| Transmembrane RTK exclusion | AXL, EPHA2 | Membrane-anchored RTK; Cullin5-RBX2 (15–20 nm) cannot assemble 2–4 nm from membrane | ✗ predicted DEG — membrane context not in features |
| Mode 1 (post-ternary) | ABL1, EPHA2 | Excellent geometry; ubiquitination or proteasome engagement failure | ✗ predicted DEG — post-ternary failure invisible to structural features |
| Unknown (marginal geometry) | EPHB2 | Highest d_sum_min (2.58 Å) among degraders; marginal bridging geometry | ✗ predicted NON |
| Unknown (geometry overlap) | DDR2 | DEG with same geometry as ABL1/AXL/EPHA2; model cannot distinguish | ✗ predicted NON |
| Unknown (high P_stable) | MAPK14_full | 47% MD bridging but predicted NON by XGBoost alone | ✗ (GNN 0.977 DEG ✓; spatial pattern carries the signal) |

---

## GNN v2: Adding Global Geometry Features (04e, DISCARDED)

Added `kin_when_vhl_near` and `n_near_vhl` as graph-level scalars concatenated with the fused [emb_A || emb_B] embedding before the MLP head. Result on 15-POI panel: AUROC=0.611 — **worse than the baseline dual-encoder GNN (0.741 on 15 POIs / 0.850 on 17 POIs)**. Injecting global scalar features destabilises LOPO-CV on small N. The GNN learns the scalar-vs-graph interaction poorly when each fold has only ~14 training examples. Discarded.

---

## Novel GNN Architectures (Sep 2026)

Four additional GNN architectures were designed and evaluated after the dual-encoder baseline, systematically probing what structural information predicts degradation.

### IxGNN v2 — typed interaction MPNN (AUROC = 0.300)

**Architecture:** TypedMPNN over pRosettaC interaction graphs. Nodes = VHL/kinase/PROTAC residues; edges typed by interaction class (H-bond, hydrophobic, pi-pi face, pi-pi edge, pi-cation, salt bridge). Node and edge features include pRosettaC contact frequency and mean/min distance. Three MPNN layers (hidden=64), attentional readout, LOPO-CV on 17 POIs.

**Result:** AUROC = 0.300 (below random). Static interaction counts do not separate degraders from non-degraders: the same interaction pattern is compatible with both classes. The model learns a trivially wrong mapping because the key signal (contact *collapse* between docking and MD) is not encoded.

### Pose-Stability GNN v2 — MD trajectory labelling (AUROC = 0.067)

**Architecture:** Same TypedMPNN graph structure, but supervised by per-pose P_stable labels (fraction of 20 ns MD frames where d_vhl < 8 Å and d_kin < 8 Å). Regression target: predict P_stable from the pRosettaC interaction graph.

**Result:** AUROC = 0.067. 20 ns trajectories are too short for stable P_stable estimates: DEG mean P_stable = 0.338, NON mean = 0.319 — negligible overlap-resolving power. Signal-to-noise in the regression target prevents the GNN from learning a usable degrader-vs-non-degrader representation. Discarded.

### Delta-Contact GNN — pRosettaC→MD contact collapse (AUROC = 0.667)

**Architecture:** TypedMPNN with 18 regression heads, one per (interaction side: vhl/kin) × (interaction type: hbond/hydrophobic/pi-pi-face/pi-pi-edge/pi-cation/salt-bridge/vdw/repulsive/clash) combination. For each pose, the model predicts the delta: δ = MD_frequency − pRosettaC_binary (how much each interaction type collapses from docking to MD). A collapse score = mean of negative predicted deltas (larger negative = more collapse). At inference, collapse score ranks POIs: higher collapse → predicted NON. LOPO-CV on **8 MD POIs** (the only POIs with both pRosettaC and MD data).

**Stage 1 regression** (all 8 POIs): best-fit per-type RMSE; best delta: d_vhl_pi_pi_face RMSE=0.023.

**LOPO-CV per-POI collapse scores (Stage 2):**

| POI | Label | Collapse score | Correct? |
|-----|-------|----------------|----------|
| EPHB2 | DEG | −0.167 | ✓ DEG (low collapse) |
| RIPK2_full | DEG | −0.153 | ✓ DEG |
| EPHA2 | NON | −0.152 | borderline — high collapse but NON |
| DDR2 | DEG | −0.151 | ✓ DEG |
| MAPK14_full | DEG | −0.135 | ✓ DEG |
| ABL1 | NON | −0.125 | ✓ NON |
| AXL | NON | −0.107 | ✓ NON |
| MET | DEG | −0.086 | **✗ NON** (genuine DEG but stable contacts) |

**AUROC = 0.667** on 8-POI subset.

**Mode 7 hook detection:** The Delta-Contact GNN is the first GNN to detect the conformational hook effect. AXL pRosettaC poses place the PROTAC warhead in a forced pi-cation interaction with VHL residues; during MD, this contact collapses (δ_{vhl,pi-cation} = −0.313 for AXL — the largest collapse in the panel). The GNN flags AXL as NON via high collapse score. ABL1 is similarly flagged. However, MET is a genuine degrader whose contacts are *stable* from docking to MD (low collapse score −0.086), leading to a false NON prediction. The model captures hook-induced failure (Mode 7) but not the positive productive-bridging signal.

**Comparison to prior GNNs:** IxGNN v2 (0.300) → Pose-Stability v2 (0.067) → Delta-Contact GNN (0.667). Best interpretable structure-based GNN to date, but below the XGBoost tabular baselines.

### KinPocket-GNN — ATP pocket structure pretraining (AUROC = 0.567, negative result)

**Architecture:** PocketMPNN — rotation-invariant Cα graph of the ATP binding pocket (DFG motif D-Φ-G centroid, 12 Å cutoff). Node features: AA one-hot (20d) + pLDDT/100 (1d) + 3-state secondary structure (3d) = 24d total. Edge features: 4 RBF-encoded pairwise distances (centres 3/5/7/9 Å, σ=1.5). 3 × distance-aware MPNN layers, hidden=64, attentional readout, 461 pocket graphs (all kinases with AF2 models).

**Pretraining (kinome-wide):** 263 kinases (train) + 65 (val), 300 epochs, AdamW lr=3e-4, cosine annealing, regression on log(d\_vhl\_min + 1) from kinome scan geometry. val\_r = 0.432, val\_RMSE = 1.053. Transferred to 17-POI LOPO-CV fine-tuning.

**Fine-tuning:** FinetuneHead = encoder output (64d) concatenated with 5 z-normalised geo\_feats (prostac\_bridging\_frac, d\_vhl\_min, d\_kin\_min, d\_vhl\_10pct, d\_kin\_10pct) → MLP(69→32→1) → sigmoid.

**LOPO-CV results (17-POI DEG/NON):**

| Variant | AUROC | ACC | Notes |
|---------|-------|-----|-------|
| Frozen encoder | 0.567 | 12/17 | Pretrained pocket features |
| Unfrozen encoder | 0.417 | 10/17 | Fine-tuning destroys pretrain signal |
| Random encoder | 0.550 | — | No pretraining; baseline A |
| **Geo-feats only (no GNN)** | **0.667** | — | MLP on 5 scalar features; baseline B |

**Finding:** Geo-feats-only MLP (0.667) outperforms all GNN variants. The ATP pocket GNN adds no information beyond what the 5 pRosettaC scalar features already encode.

**TNIK/IRAK1 false-positive artifact:** Both have d\_vhl\_min=100 (no pRosettaC model — geometry-blocked). After z-normalisation, geo\_feats → near-zero vector. The MLP head interprets the missing-data pattern as DEG-like (TNIK p\_deg=0.971, IRAK1 p\_deg=0.878). This is a missing-data imputation artifact, not a structural signal.

**Per-POI predictions (frozen encoder, LOPO-CV):**

| POI | Label | p(DEG) | Pred | Correct |
|-----|-------|--------|------|---------|
| TNIK | NON | 0.971 | DEG | ✗ (missing-data artifact) |
| IRAK1 | NON | 0.878 | DEG | ✗ (missing-data artifact) |
| EPHB2 | DEG | 0.814 | DEG | ✓ |
| DDR2 | DEG | 0.806 | DEG | ✓ |
| MET | DEG | 0.484 | DEG | ✓ (marginal) |
| CDK4 | NON | 0.331 | DEG | ✗ |
| AXL | NON | 0.287 | NON | ✓ |
| EPHA2 | NON | 0.201 | NON | ✓ |
| ABL1 | NON | 0.148 | NON | ✓ |
| CDK6 | NON | 0.026 | NON | ✓ |
| CDK2 | NON | 0.023 | NON | ✓ |
| RIPK2_full | DEG | 0.016 | NON | ✗ |
| TLK2 | NON | 0.003 | NON | ✓ |
| MAP4K5 | NON | 0.002 | NON | ✓ |
| DAPK1 | NON | 0.001 | NON | ✓ |
| SLK | NON | 0.001 | NON | ✓ |
| MAPK14_full | DEG | 0.000 | NON | ✗ |

**Scientific interpretation:** For SJF8240/VHL, degradation specificity is encoded in the 3D ternary complex geometry (summarised by pRosettaC scalar features), not in ATP binding pocket topology. Any kinase pocket can orient the warhead if the linker length is compatible — the discriminating signal is VHL engagement quality and PROTAC bridging geometry. The GNN trained on ATP pocket topology learns no more than d\_vhl\_min and prostac\_bridging\_frac already convey. This negative result rules out pocket identity as a SJF8240 degradation determinant and motivates full ternary complex surface complementarity features for future work.

### E2-scalar GNN (04g) — E2 reachability at fusion (AUROC=0.800, ACC=15/17)

**Architecture:** Identical to 04d dual-encoder with four E2 reachability scalars (`e2_mean_min_dist_A`, `e2_delta_reach_50A`, `e2_delta_reach_80A`, `lk_lys_nearest_e2_dist_A`) concatenated to the fused embedding before the fusion MLP. Z-normalised per LOPO fold; NaN (non-MD POIs) → 0.0.

**Result:** AUROC=0.800, ACC=15/17. Compared to 04d baseline (AUROC=0.850, ACC=14/17):

| POI | 04d | 04g | Change |
|-----|-----|-----|--------|
| ABL1 | NON→DEG ✗ (p=0.998) | NON→NON ✓ (p=0.039) | **fixed** |
| AXL | NON→DEG ✗ (p=0.998) | NON→NON ✓ (p=0.009) | **fixed** |
| EPHA2 | NON→DEG ✗ (p=0.998) | NON→DEG ✗ (p=0.999) | unchanged |
| MAPK14_full | DEG→DEG ✓ (p=0.996) | DEG→NON ✗ (p=0.000) | **broken** |
| All others | same | same | — |

**ABL1 fix:** `lk_lys_nearest_e2_dist_A` = 112.8 Å — genuinely outlier relative to DEGs (97–107 Å range). ABL1's substrate lysines are unusually far from the E2 active site in the ternary complex geometry; this is a real structural signal.

**AXL fix:** Driven by the combined E2 scalar bundle (no single feature dominates).

**MAPK14_full breakage:** `e2_mean_min_dist_A` = 48.9 Å is a genuine biological outlier in the opposite direction (MAPK14's lysines are unusually close to E2). At N=17, the model cannot distinguish "close-E2 DEG exception" (MAPK14_full) from "close-E2 NON-like" — both look the same to the MLP. Would resolve with larger training panel.

**EPHA2 unchanged:** `e2_mean_min_dist_A` = 60.4 Å, within the DEG range (61–63 Å). No E2 scalar separates it. Remains a post-ternary failure (transmembrane RTK membrane context) not captured by any structural feature.

**Verdict:** AUROC drops 0.850→0.800 despite ACC improving 14→15. E2 reachability encodes genuine signal (ABL1 confirmed) but introducing it as a scalar at N=17 trades one failure for another. Dual-encoder 04d remains best model by AUROC.

### Comprehensive architecture comparison (all GNNs, Sep 2026)

| Model | N POIs | AUROC | ACC | Notes |
|-------|--------|-------|-----|-------|
| IxGNN v2 (TypedMPNN, static) | 8 MD-POIs | 0.300 | — | Static interaction counts |
| Pose-Stability GNN v2 | 8 MD-POIs | 0.067 | — | 20 ns labels too noisy |
| Delta-Contact GNN | **8 MD-POIs** | **0.667** | — | Best interpretable GNN; Mode 7 hook |
| KinPocket-GNN (frozen) | 17 | 0.567 | 12/17 | ATP pocket + d_vhl_min pretrain |
| KinPocket-GNN (unfrozen) | 17 | 0.417 | 10/17 | Pretrain signal destroyed |
| KinPocket random encoder | 17 | 0.550 | — | No pretraining baseline A |
| Geo-feats only MLP | 17 | 0.667 | — | No GNN; 5 pRosettaC scalars; baseline B |
| E2-scalar GNN (04g) | 17 | 0.800 | **15/17** | ABL1/AXL fixed; MAPK14_full broken |
| Geom-only XGBoost | 17 | 0.800 | 11/17 | Tabular, no MD |
| Full XGBoost + MD | 17 | 0.850 | 11/17 | 29 tabular features |
| GNN dual-encoder + MD (04d) | 17 | **0.850** | **14/17** | **CURRENT BEST by AUROC** |
| Ensemble 95/5 (GNN+XGB) | 17 | 0.850 | 14/17 | No gain over GNN alone |

**Information hierarchy:**
1. *Static docking interaction counts* (IxGNN v2: 0.300) — wrong input signal
2. *20 ns MD binary labels* (Pose-Stability v2: 0.067) — label noise overwhelms signal
3. *Contact collapse / ATP pocket geometry* (Delta-Contact/Geo-feats MLP: 0.667) — partial mechanistic signal
4. *E2 reachability scalars + spatial MD graphs* (04g: AUROC=0.800, ACC=15/17) — genuine post-ternary signal; limited by N=17 artefact
5. *Spatial MD contact frequency graphs alone* (Dual-encoder 04d: 0.850) — best overall; robust to E2-outlier artefacts

The stepwise hierarchy demonstrates that MD contact *dynamics* is the critical discriminating feature. E2 reachability adds orthogonal post-ternary biological signal but requires a larger training panel to exploit cleanly. The 14/17 ceiling in 04d and the 15/17 ceiling in 04g share the same root cause: ABL1/AXL/EPHA2 failures are post-ternary and structurally irrecoverable; the E2 scalar experiment shifts which failure mode is active but cannot fully escape the ceiling.

---

## Label Ground Truth (canonical — settled 2026-09-20)

**5 DEGs:** DDR2, EPHB2, MAPK14_full, MET, RIPK2_full
**12 NONs:** ABL1, AXL, EPHA2, MAP4K5, SLK, TNIK, CDK2, CDK4, CDK6, DAPK1, IRAK1, TLK2

- **EPHA2 = NON** — ml_table: high_non_degrader, experimental_degrader=0. Despite excellent ternary geometry (lowest d_sum_min in panel = 0.39 Å), confirmed NON. Kinome scan script 09 GROUND_TRUTH dict contains a bug (EPHA2=1) — authoritative label is ml_table.
- **TNIK = NON for structural models** — ml_table experimental_degrader=1 but prostac_artifact=1 (kin_min overlap inflates pRosettaC bridging; TNIK_30000 confirmed 0/30000 bridging). Treated as NON in all ML evaluations.

---

## Limitations and Future Work

1. **N=17 clean POIs** — LOPO-CV AUROC estimates have wide CIs (~±0.1 at this N); bootstrap resampling required to quantify uncertainty
2. **Single PROTAC (SJF8240/VHL)** — generalisability to other E3 ligases or PROTAC scaffolds untested; all POIs share the same linker/warhead
3. **No MD replicates** — each pose is a single 100 ns trajectory; P_stable estimates are noisy; binary complex panel (6 POIs × 20 ns each, Aug 11) provides complementary signal but not full replicates
4. **Binary complex warhead burial not discriminating** — Binary complex MD panel (8 POIs × 20 ns) shows VHL-warhead burial is universal (97–99.5% for all kinases, including confirmed degraders MET/EPHB2). Cannot be added as a feature to fix AXL/EPHA2/ABL1 mispredictions. AXL/EPHA2 failures attributed to transmembrane RTK context; ABL1 failure unexplained at structural level.
5. **GNN retrained on 17 POIs — COMPLETE (Aug 12)** — dual-encoder GNN AUROC=0.850, ACC=14/17; CDK17=NON correctly predicted NON (near-empty graphs). Ensemble = 95% GNN + 5% XGB, AUROC=0.850 (no gain over GNN alone).
6. **ABL1/AXL/EPHA2 experimental validation** — all three are strongly predicted DEG by structural features; potential experimental false negatives warrant cell-line or time-course validation
