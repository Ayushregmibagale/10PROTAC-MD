# 06 — Results Summary

**Last updated:** 2026-06-03  
**Status:** MET complete; DDR2 pipeline-check only; 11 POIs queued

---

## MET — Complete (high-affinity degrader, reference)

### Run details
- Structure: 3LQ8 (2.02 Å, foretinib/XL880 = actual SJF8240 kinase warhead)
- N_DOCK = 200, N_WORKERS = 4 (initial run)
- Total docking poses: 10,000 (200 global × 50 local)
- Successful combined PDBs: ~6,000 (60%)
- Clusters produced: **6**
- MD candidates (md_verdict = READY): **3** (cluster1, cluster2, cluster6)

### Cluster summary

| Cluster | Poses | Interface contacts (5Å) | Polar H-bonds | Score (REU) | MD candidate |
|---------|-------|------------------------|---------------|-------------|--------------|
| cluster1 | 50 | 88 | 1 | 138.7 | Primary (1) |
| cluster2 | 49 | **194** | 3 | 129.0 | Secondary (2) |
| cluster3 | 33 | 106 | 4 | 164.1 | — |
| cluster4 | 28 | 115 | 2 | 159.2 | — |
| cluster5 | 21 | **243** | 5 | 161.1 | — |
| cluster6 | 19 | 170 | 6 | 159.0 | Tertiary (3) |

### Key features (primary MD candidate = cluster1)
- VHL structure completeness: 98.1% (2 small gaps, none near interface)
- MET kinase completeness: 94.9% (15 residues missing, none near interface)
- PROTAC: 79 atoms, end-to-end distance = 37.1 Å, not disconnected ✓
- Interface clashes: 0 (severe), 0 (bad), 0 (warn) ✓
- Closest lysine: K1161 at 3.5 Å from VHL surface
- n_lys_e2_optimal: **16** (16 MET lysines in 35–55 Å range from E2 active site)
- Primary E2-optimal lysine: K1161 at 35.9 Å
- ElonginB present: True ✓, ElonginC present: True ✓

### Interpretation
MET shows excellent ternary complex geometry:
- Cluster2 (secondary MD candidate) has 194 interface contacts and the lowest Rosetta
  score (129.0 REU) — this is likely the most stable and physically relevant geometry
- 16 optimal lysines in ubiquitination range provides multiple options for ubiquitination
- No clashes, complete structure, PROTAC intact — all quality checks pass
- **This is exactly what we expect for a strong degrader** ✓

---

## DDR2 — Pipeline check only (NOT for analysis)

### Run details
- Structure used: DDR1/3ZOS (WRONG — replaced by 7AYM, needs clean re-run)
- N_DOCK = 200, N_WORKERS = 4
- Completed: 2026-06-03 14:31 BST
- Results/: 186 combined PDB structures across clusters

### Status
This run **cannot be used for scientific conclusions** because:
1. The protein structure was DDR1 (3ZOS), not DDR2 (7AYM)
2. The warhead SDF was changed mid-run (at ~86% completion) to the 7AYM-based placement
3. The Results/ contains structures from a mixed DDR1/warhead run

**This run was useful for:** confirming the full pipeline (docking → constraint_gen →
clustering → Results/) runs to completion without errors on this system.

**Action required:** Re-run DDR2 with corrected 7AYM structure.

---

## All other POIs — Queued

RIPK2, EPHB2, MAPK14, AXL, SLK, ABL1, EPHA2, MAP4K5, TNIK, PIP4K2C, CDK17 have not
been run yet. Input files (Protac_params.txt, warhead SDFs, protein_clean.pdb) are ready.

---

## What to expect when runs complete

### Expected features for degraders (based on MET and pRosettaC literature)
- Interface contacts: > 80–100 contacts at 5 Å
- n_lys_e2_optimal: ≥ 5 lysines in the 35–55 Å range
- Primary E2 distance: 35–45 Å (tighter = better)
- md_verdict: READY (no severe clashes)

### Expected features for non-degraders
- Interface contacts: < 50 (less complementary surfaces)
- n_lys_e2_optimal: 0–3 (few or no accessible lysines)
- OR: interface contacts moderate but no lysine in range
- May still have some clusters but with worse geometry

### The low-affinity degraders (TNIK, PIP4K2C, CDK17)
- May show moderate interface contacts (degradation is possible just slow/partial)
- Likely fewer lysines in range, or larger E2–lysine distances
- Will serve as an independent test set for the high-affinity model

---

## ML classification plan (to run after all POIs complete)

1. **Primary model**: logistic regression on 10 high-affinity POIs (5 degraders + 5
   non-degraders) with LOO cross-validation
2. **Features**: `CORE_FEATURES = ['best_interface_contacts_5A', 'best_n_lys_e2_optimal',
   'best_polar_hbond_pairs', 'md1_best_e2_A']`
3. **Validation**: LOO accuracy on the 10 high-affinity POIs; apply to 3 low-affinity
   degraders as held-out test
4. **Decision tree**: depth-2 tree for interpretable rules
5. **Expected output**: classification summary CSV + printed confusion matrix

---

## Timeline

| Date | Event |
|------|-------|
| 2026-05-25 | Project started |
| 2026-05-25 to 06-02 | MET run and analysis completed |
| 2026-06-03 | BioPython rotation bug found and fixed; all warheads regenerated |
| 2026-06-03 | DDR2 structure corrected (DDR1 → DDR2/7AYM) |
| 2026-06-03 | EPHB2 structure updated (ADP-bound 2HEN → apo 3ZFM) |
| 2026-06-03 | DDR2 pipeline-check run completed (using old DDR1 structure) |
| TBD | All remaining POIs run (HPC or local) |
| TBD | Cross-POI analysis and ML classification |
