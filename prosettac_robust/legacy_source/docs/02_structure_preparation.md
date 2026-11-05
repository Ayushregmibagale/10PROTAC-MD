# 02 — Structure Preparation

## Overview

For each POI we need two cleaned protein structures:
1. **VHL:ElonginB:ElonginC complex (VBC)** — the E3 ligase half, with SJF8240 VHL warhead placed
2. **POI kinase domain** — with SJF8240 kinase warhead placed

Both structures feed into pRosettaC as rigid bodies. The warhead positions define the
pharmacophoric constraints that the PROTAC linker must satisfy.

---

## VHL structure

**PDB: 4W9H** — 2.10 Å resolution  
Title: pVHL:ElonginB:ElonginC in complex with VHL ligand (2S,4R)-4-hydroxyl-proline scaffold  
Chains used: A (VHL), B (ElonginB), C (ElonginC)  
Ligand in PDB: 3JF (VHL-targeting small molecule — same chemotype as SJF8240 VHL warhead)

**Why 4W9H:**
- High resolution (2.10 Å) — excellent electron density for the binding pocket
- Contains the entire VBC E3 complex needed for pRosettaC (VHL + both elongins)
- The bound ligand (3JF) is structurally nearly identical to SJF8240's VHL warhead
  (same (2S,4R)-4-hydroxyproline core, same thiazole-benzyl group); the warhead's
  binding pose is directly transferable
- Standard structure in the published PROTAC modelling literature (used by Ciulli group)
- Alternative 4W9H derivatives with different VHL ligands give essentially identical
  protein coordinates (RMSD < 0.3 Å in the binding pocket)

**Cleaning steps applied:**
1. Extracted chains A, B, C (VBC complex)
2. Removed water molecules and crystallographic additives
3. Retained bound ligand position as reference for warhead placement (then removed ligand)
4. Saved as `data/structures/VHL/4W9H_VBC_clean.pdb`

**VHL warhead placement:**
- SJF8240 VHL warhead SMILES mapped onto the 3JF binding pose via maximum common
  substructure (MCS) overlap
- Both molecules share the (2S,4R)-4-hydroxyproline core + thiazole-benzyl group
- MCS coverage: the core pharmacophore atoms are fixed; flexible substituents were
  energy-minimised (MMFF94) around the fixed MCS positions
- Output: `data/structures/VHL/SJF8240_vhll_placed.sdf` (no H) and
  `SJF8240_vhll_placed_H.sdf` (with H, for constraint_generation)

---

## Kinase structures

The foretinib pharmacophore (SJF8240 kinase warhead) was extracted from MET/3LQ8 and
mapped into each target kinase frame via Cα superposition.

### MET — reference structure

**PDB: 3LQ8** — 2.02 Å  
Title: Structure of the kinase domain of c-MET bound to XL880 (GSK1363089 / foretinib)  
Ligand: 88Z (XL880 = foretinib analog)  
Chain used: A (kinase domain, residues ~956–1382; used residues 1053–1345 for analysis)

**Why 3LQ8:**
- Contains the exact same pharmacophore as SJF8240's kinase warhead (XL880 / 88Z)
- MET is the primary target of foretinib; this structure has the highest-fidelity
  binding pose for this warhead class
- Good resolution (2.02 Å); electron density clearly defines warhead orientation
- Used as the REFERENCE for all other warhead placements via Cα superposition

**Warhead extraction:**
- 88Z heavy-atom coordinates extracted from 3LQ8 HETATM records
- MCS computed between 88Z SMILES and SJF8240 kinase warhead SMILES
- 73% warhead coverage by MCS (27/37 heavy atoms share the same scaffold core)
- Non-MCS atoms placed by ETKDG embedding followed by MMFF94 minimisation

---

### DDR2

**PDB: 7AYM** — 2.12 Å (UPDATED from 3ZOS during this campaign)  
Title: Structure of DDR2 kinase domain in complex with IBZ3 (compound 4VD)  
Ligand: 4VD (IC50 = 3.3 nM DDR2 inhibitor)  
Chain used: A (kinase domain, residues 553–858)

**Why 7AYM (and why 3ZOS was wrong):**
- 3ZOS (previously used) is the **DDR1** kinase domain — a distinct paralogue with ~52%
  kinase domain sequence identity to DDR2. Using it introduces systematic error in the
  binding pocket geometry.
- 7AYM is a genuine **DDR2** crystal structure at excellent resolution (2.12 Å)
- The bound inhibitor 4VD is a DDR2-selective compound — the binding pocket is in a
  relevant, closed inhibitor-bound conformation
- Superposition RMSD onto MET/3LQ8 with 7AYM: **0.83 Å** (162 Cα pairs) — very good

**Cleaning:** ATOM records only (chain A), HETATM (4VD inhibitor) removed.

---

### RIPK2

**PDB: 4C8B** — 2.75 Å  
Title: Structure of the kinase domain of human RIPK2 in complex with ponatinib  
Ligand: 0LI (ponatinib)  
Chain used: A

**Notes:**
- Correct RIPK2 crystal structure — no issue with paralogue confusion here
- The bound ligand is ponatinib (a type II kinase inhibitor), not foretinib
- However, the protein structure is what matters for warhead placement; the bound
  ligand is removed and the warhead is placed by Cα superposition from MET
- 2.75 Å resolution is adequate; B-factors are acceptable
- Superposition RMSD: moderate (RIPK2 is a serine/threonine kinase, more divergent
  from MET than receptor tyrosine kinases)

---

### EPHB2

**PDB: 3ZFM** — 2.27 Å (UPDATED from 2HEN during this campaign)  
Title: Crystal structure of EphB2 (apo kinase domain)  
Ligand: none (apo structure)  
Chain used: A (residues 615–894)

**Why 3ZFM (and why 2HEN was problematic):**
- 2HEN (previously used) has ADP bound in the ATP site. ADP-bound kinases are in the
  active (DFG-in) conformation. While not catastrophic, ADP:
  (a) occupies the binding site where the foretinib warhead should be placed
  (b) biases the activation loop towards the active conformation, which may not reflect
  how foretinib actually binds EphB2
- 3ZFM is the first published wild-type human EphB2 kinase domain structure (apo)
- Without a bound nucleotide or inhibitor, the binding pocket is "open" and unbiased
- 2.27 Å resolution; no gaps in the kinase domain
- There is NO EphB2 crystal structure with a type I or type II ATP-competitive inhibitor
  in the PDB as of June 2026 — 3ZFM is the best available option
- Superposition RMSD: 2.47 Å (209 Cα pairs) — reasonable for EphB2 vs MET

---

### MAPK14 (p38α)

**PDB: 3GCP** — 2.25 Å  
Title: Human p38 MAP kinase in complex with SB203580  
Ligand: SB2 (SB203580, a type I p38 inhibitor)  
Chain used: A (residues 4–352; first 103 residues overlap with VHL numbering — only
residues 104–352 contribute to POI-side interface in analysis)

**Notes:**
- Correct protein (MAPK14 = p38α); well-validated structure
- SB203580 is a different scaffold than foretinib but shares the same ATP-binding pocket
- The protein conformation in 3GCP is DFG-in (active-like, type I inhibitor bound)
- Foretinib's binding mode to p38 is less well-characterised than to MET/EphB2;
  the DFG-in conformation is a reasonable assumption for initial modelling
- Warhead placement note: MAPK14 has the largest structural divergence from MET
  (RMSD ~12 Å over the full kinase domain), so warhead placement relies on the
  sequence alignment picking the correct ATP-pocket residues. Manual verification
  in PyMOL is recommended after the run completes.

---

### AXL

**PDB: 5U6B** — 2.84 Å  
Title: Structure of the AXL kinase domain in complex with a macrocyclic inhibitor  
Ligand: 7YS  
Chain used: A (residues 521–813)

**Notes:**
- Correct AXL crystal structure with an ATP-competitive inhibitor bound
- 2.84 Å resolution is on the lower end; acceptable for warhead placement purposes
- AXL is a receptor tyrosine kinase closely related to MET; superposition quality good
- Macrocyclic inhibitor keeps the ATP pocket in a defined conformation

---

### SLK

**PDB: 2J51** — 2.10 Å  
Title: Crystal structure of human STE20-like kinase bound to a kinase inhibitor  
Ligand: DKI  
Chain used: A (residues 21–308; first 20 residues overlap with VHL — only 104–308 used)

**Notes:**
- Good resolution (2.10 Å) with an ATP-competitive inhibitor bound
- SLK has limited structural data; 2J51 is the most commonly used structure
- Residues 21–103 overlap numerically with VHL chain A (1–103) in the combined pRosettaC
  structure — the analysis script uses poi_min=104 to exclude these

---

### ABL1

**PDB: 2HYY** — 2.40 Å  
Title: Human ABL kinase domain in complex with imatinib (STI571/Gleevec)  
Ligand: STI (imatinib — a classic type II kinase inhibitor)  
Chain used: A (residues 235–498)

**Notes:**
- Imatinib is a type II inhibitor (DFG-out conformation) — same binding mode class as
  foretinib. This makes 2HYY particularly appropriate: the binding pocket geometry
  matches what foretinib would induce.
- Very well-validated structure (imatinib binding to ABL is textbook structural biology)
- Good resolution (2.40 Å)

---

### EPHA2

**PDB: 5I9Y** — 1.23 Å (excellent)  
Title: Crystal structure of Ephrin A2 (EphA2) receptor protein kinase with dasatinib  
Ligand: 1N1 (dasatinib)  
Chain used: A (residues 600–896)

**Notes:**
- Outstanding resolution (1.23 Å) — the best in the panel
- Dasatinib is a type I inhibitor (ATP-binding site, DFG-in)
- EphA2 kinase domain is closely related to EphB2 (which is a non-degrader in our panel)
- EphA2 vs EphB2 both bind foretinib; using an EphA2 structure for EphA2 is correct
- Initial concern: 3MBW was present in the folder but is the EphA2 **extracellular**
  domain (LBD/CRD) — completely wrong. 5I9Y (kinase domain) is correct and was
  confirmed to be what EPHA2_protein_clean.pdb is derived from.

---

### MAP4K5

**PDB: AlphaFold2 Q9Y4K4 v6** — no experimental structure available  
Mean pLDDT: 82.2 (kinase domain); 258/320 residues > 70 pLDDT  
Chain used: A (residues 1–320; first 103 overlap with VHL — only 104–320 used)

**Notes:**
- No crystal structure for MAP4K5 exists in the PDB as of June 2026
- The AF2 v6 prediction has high confidence in the kinase domain (pLDDT > 70 throughout)
- Kinase domains are among the best-predicted fold classes by AF2
- This limitation must be explicitly disclosed in any publication
- The warhead placement will be less reliable than for experimentally determined structures;
  the ATP-binding pocket geometry may deviate from the true structure

---

### TNIK

**PDB: 5D7A** — used for TNIK kinase domain  
Chain used: A (residues 11–310; first 10 overlap with VHL — only 104–310 used)

---

### PIP4K2C

**PDB: 2GK9**  
Chain used: A (residues 45–417; first 44 overlap with VHL — only 104–417 used)

**Notes:**
- PIP4K2C is a phosphoinositide kinase, not a protein kinase — different fold than the
  rest of the panel. The ATP-binding site is analogous but not identical.
- The foretinib warhead placement via Cα superposition will be less reliable here
- Manual inspection of the placed warhead is particularly important for this target

---

### CDK17

**PDB: AlphaFold2 O94921 v6**  
Mean pLDDT: 65.6 (kinase domain) — notably lower confidence than MAP4K5  
188/350 residues > 70 pLDDT

**Notes:**
- No crystal structure available; AF2 model has lower confidence than MAP4K5
- CDK17 (PCTK3) is a cyclin-dependent kinase — should bind foretinib
- Results from this target should be interpreted with the most caution in the panel

---

## Summary table

| POI | PDB | Resolution | Ligand class | Category | Issue |
|-----|-----|-----------|-------------|----------|-------|
| MET | 3LQ8 | 2.02Å | Foretinib (reference) | Degrader | None — ideal |
| DDR2 | 7AYM | 2.12Å | DDR2-selective (type I) | Degrader | None after fix |
| RIPK2 | 4C8B | 2.75Å | Ponatinib (type II) | Degrader | Resolution |
| EPHB2 | 3ZFM | 2.27Å | Apo | Degrader | No inhibitor in PDB |
| MAPK14 | 3GCP | 2.25Å | SB203580 (type I) | Degrader | High RMSD vs MET |
| AXL | 5U6B | 2.84Å | Macrocyclic | Non-degrader | Resolution |
| SLK | 2J51 | 2.10Å | DKI | Non-degrader | None |
| ABL1 | 2HYY | 2.40Å | Imatinib (type II) | Non-degrader | None |
| EPHA2 | 5I9Y | 1.23Å | Dasatinib (type I) | Non-degrader | None — best resolution |
| MAP4K5 | AF2 | — | AF2 model | Non-degrader | No exp structure |
| TNIK | 5D7A | — | — | Low-aff. degrader | — |
| PIP4K2C | 2GK9 | — | — | Low-aff. degrader | Non-kinase fold |
| CDK17 | AF2 | — | AF2 low-conf | Low-aff. degrader | Low pLDDT |

---

## Warhead placement script

All placements are done by `scripts/structure_prep/place_warhead.py`.

**Algorithm:**
1. Extract 88Z coordinates from `MET_3LQ8_chainA.pdb`
2. Find MCS between 88Z SMILES and SJF8240 kinase warhead SMILES; map 88Z 3D coordinates
   onto warhead via MCS atom correspondence; embed remaining atoms (ETKDG → MMFF94)
3. For each non-MET target: extract all Cα atoms; run pairwise sequence alignment
   (BLOSUM62, global); collect aligned Cα pairs; iteratively reject outliers (5 cycles,
   2σ cutoff); compute Kabsch superposition with remaining pairs
4. Apply inverse transform: `target_coord = rot @ (met_coord - tran)` where `rot` and
   `tran` are from BioPython Superimposer called as `set_atoms(fixed=MET, moving=target)`
5. Write warhead as SDF (no H) and SDF (with H) for use in pRosettaC

**BioPython convention (critical):**
`set_atoms(fixed=MET, moving=target)` returns `(rot, tran)` where the convention is
RIGHT-MULTIPLYING ROW VECTORS: `rot.T @ target + tran ≈ MET`. The inverse transform
(MET coordinates → target frame) is therefore `rot @ (met_coord - tran)`, NOT
`rot.T @ (met_coord - tran)`. See `docs/05_bugs_and_fixes.md` for the full story.

**To re-run placements:**
```bash
cd ~/protac
python scripts/structure_prep/place_warhead.py              # all POIs
python scripts/structure_prep/place_warhead.py --force DDR2 EPHB2  # specific POIs only
```
