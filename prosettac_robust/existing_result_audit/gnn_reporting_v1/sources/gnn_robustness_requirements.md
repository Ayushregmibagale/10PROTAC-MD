# GNN Robustness Requirements
*Written: 2026-08-14*

---

## Is N=17 Robust?

**No — N=17 is proof-of-concept only.**

- Bootstrap 95% CI for GNN: **[0.393, 1.000]** — lower bound nearly touches chance (0.5)
- 5 DEG / 12 NON: LOPO-CV trains on as few as 4 degraders per fold
- GNN has hundreds of parameters (GATConv ×3 layers, 4 heads, hidden=64) — severely overparameterised at N=17 without pretraining
- No true held-out test set — LOPO-CV on the same 17 POIs used for pipeline development

**What IS robust at N=17:**
- No-MD subgroup (geometry-blocked, 9 POIs): AUROC=1.000 [1.000, 1.000] — pRosettaC perfectly separates Mode 0/1; this holds at any N
- Mechanistic findings (AXL hook, TNIK triply improbable, failure modes) are structural results independent of ML sample size
- XGBoost and GNN agree (AUROC=0.850 both) — two independent model families converging is meaningful signal

---

## How Many POIs Are Needed for Robustness?

| N (labelled POIs) | Status |
|---|---|
| 17 | Proof-of-concept; CIs too wide for generalisation claims |
| ~50 | Minimum for meaningful GNN training; CIs tighten to ~±0.15 |
| ~200 | Publication-grade robustness with tight CIs |
| ~500+ | Field-standard for structure-based PROTAC ML |

**Practical next step:** 5-PROTAC dataset — train on 4 PROTACs × 15 POIs = 60 points, test on 1 held-out PROTAC × 15 POIs = 15 points. Leave-one-PROTAC-out tests generalisation across chemical scaffolds — more meaningful than LOPO-CV.

---

## Pretraining vs MD — Completely Different Signals

**MGDB pretraining (3,921 graphs, one-time, already done):**
- Teaches the encoder general ternary complex topology from *static* structures (PDB + Boltz-2 predicted)
- Transfer learning: encoder starts from sensible weights rather than random noise
- Does **NOT** reduce the MD requirement

**Why pretraining alone is not enough:**
> pRosettaC-only GNN (after MGDB pretraining) scored **AUROC = 0.185** — essentially useless.

Pretraining is on *static* graphs. The discriminating signal between degraders and non-degraders is *dynamic*: which residue contacts are stable when the PROTAC is bridging. Pretraining cannot teach this because static structures contain no temporal information.

**Summary:**
- Pretraining → "what does a ternary complex look like structurally" ✓
- MD → "which contacts hold when PROTAC bridges in dynamics" ✓ (required, cannot substitute)

---

## Which POIs Actually Need MD?

| POI type | MD needed? | Reason |
|---|---|---|
| Geometry-blocked (Mode 0/1) | **No** | pRosettaC 0% bridging → already perfectly classified by pRosettaC graph |
| Bridging-capable (degraders + Mode 2/3/4 non-degraders) | **Yes** | MD contact frequencies needed to separate stable from unstable ternary |

In a 5-PROTAC expansion (~75 POIs): ~25–30 would be bridging-capable and need MD P_stable screen (20 poses × 20 ns each ≈ 500 ns total). The majority (~70%) are geometry-blocked and need only pRosettaC.

**Rule: do NOT run MD for geometry-blocked POIs** — pRosettaC correctly predicts them and MD is wasted compute.

---

## 5-PROTAC Expansion — Full Requirements

| Component | Required? | Notes |
|---|---|---|
| MGDB pretraining | Reuse existing weights | One-time; already done on 3,921 graphs |
| Boltz-2 inference | No | Only needed if expanding MGDB with new glue structures |
| pRosettaC | Yes — 60 new runs | 4 new PROTACs × 15 POIs; each PROTAC needs its own runs (different linker geometry) |
| MD P_stable screen | Yes — ~25 POIs only | Bridging-capable subset; 20 poses × 20 ns each |
| GNN fine-tuning | Yes | Load pretrained weights → fine-tune on 60 train + 15 held-out test PROTAC |

**Bottleneck:** pRosettaC (weeks of ernie compute), not pretraining or Boltz-2.

---

## Defence Framing

> "N=17 is sufficient for hypothesis generation and mechanism discovery, but not for a deployable predictive model. The consistent AUROC=0.850 across two independent model families (XGBoost and GNN) under LOPO-CV suggests the signal is real. The geometry-blocked subgroup result (AUROC=1.000) is fully robust at N=9. Expanding to the 5-PROTAC dataset with a genuine held-out test PROTAC is the planned validation."
