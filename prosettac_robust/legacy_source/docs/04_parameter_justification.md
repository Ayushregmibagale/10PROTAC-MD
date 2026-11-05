# 04 — Parameter Justification

Every significant parameter choice is documented here with its rationale.

---

## pRosettaC core parameters

### N_DOCK = 200

**What it controls:** Number of global PatchDock docking poses retained and then subjected
to local Rosetta relaxation. Each global pose generates 50 local relaxations, so
N_DOCK=200 → 10,000 total locally docked structures.

**Why 200:**
- Directly matches the original pRosettaC publication (Zaidman et al., *Nature Chemical
  Biology*, 2020), which validated the method with N_DOCK=200
- Published retrospective studies using pRosettaC uniformly use 200 global poses
- Preliminary testing with N_DOCK=50 shows reduced cluster diversity (fewer distinct
  ternary complex geometries sampled)
- N_DOCK=500 is feasible but increases compute time ×2.5 with diminishing returns in
  cluster diversity for the size of kinase·VBC systems (~500 residues combined)
- **For a publication, N_DOCK=200 is the standard and defensible choice.**

### Full: True

**What it controls:** Toggles between full PatchDock scoring (Shape Complementarity +
Atomic Contact Energy + Geometric Desolvation) and fast approximation.

**Why True:**
- Full scoring produces higher-quality, more geometrically reasonable docking poses
- The fast mode skips the Atomic Contact Energy term, which is important for
  protein-protein interfaces
- Compute cost difference is modest (< 20% for the PatchDock stage)

### Anchor atoms: 6 1

**What it controls:** The atom indices (1-based) in the VHL warhead (PT0) and kinase
warhead (PT1) where the PROTAC linker attaches. These define the pharmacophoric
constraints: the linker must start at atom 6 of PT0 and end at atom 1 of PT1.

**Why 6 and 1 (for SJF8240 specifically):**
- PT0 (VHL warhead, 30 heavy atoms): atom 6 is N6, the amide nitrogen that connects
  to the PEG linker. This is the chemically correct attachment point based on the
  SJF8240 SMILES and crystal structure.
- PT1 (kinase warhead, 37 heavy atoms): atom 1 is C1, the carbon at the foretinib
  scaffold's point of linker attachment. This was confirmed by matching the SDF atom
  ordering against the SJF8240 SMILES decomposition.
- These values are **constant across all POIs** because the PROTAC molecule is the same;
  only the protein structures change between POIs.

**Verification:** The anchor atoms in the placed SDF files should be within a chemically
reasonable distance (~3–6 Å) from the nearest protein Cα. `place_warhead.py` reports this
distance for every POI. All POIs showed anchor-to-Cα distances of 3–8 Å (acceptable).

### Chains: ABC A

**What it controls:**
- First token (`ABC`): the combined chain string for the ternary complex input — A = VHL,
  B = ElonginB, C = ElonginC
- Second token (`A`): the chain identifier of the E2-bearing component within the VBC
  complex. This is used by clustering to compute E2–lysine distances. The E2 (UBE2D/UbcH5)
  is modelled as associated with ElonginB/ElonginC in the VHL pathway, and chain A (VHL)
  is used as the reference for distance measurements.

**Why ABC A:**
- Directly from the 4W9H VBC complex structure: chains A/B/C are VHL/ElonginB/ElonginC
- The kinase chain in the combined structure is labelled A by pRosettaC convention
  (takes the chain ID of the input POI structure)

---

## Warhead placement parameters

### Superposition algorithm: BLOSUM62, global alignment, 5-cycle outlier rejection, 2σ cutoff

**Why global alignment (not local):**
- Kinase domains are well-conserved in their overall fold; global alignment captures
  the full structural correspondence
- Local alignment risks aligning only the most similar region (e.g., the hinge) while
  missing the overall orientation, which would misplace the warhead

**Why BLOSUM62:**
- Standard substitution matrix for protein sequence alignment
- Appropriate for kinase domain sequences (30–80% identity)
- Gap penalties: open=-10, extend=-0.5 (conservative, prevents excessive gaps in
  well-aligned kinase domains)

**Why 5-cycle outlier rejection at 2σ:**
- Rigid superposition is sensitive to loop regions that diverge in position even when
  the overall fold is similar (activation loop, P-loop, αC-helix can move by 5–15 Å)
- Iterative outlier rejection removes these flexible regions from the superposition,
  giving a much better alignment of the core ATP-binding site residues
- 2σ cutoff balances removing genuine outliers without over-pruning conserved residues
- 5 cycles is sufficient for convergence in all tested kinase pairs

### MCS-based warhead mapping (ETKDG + MMFF94)

**Why MCS (Maximum Common Substructure) rather than direct coordinate transfer:**
- The reference ligand (88Z / foretinib analog from 3LQ8) has a slightly different
  scaffold than SJF8240's kinase warhead (morpholinopropoxy group differs from the
  warhead's simpler substituent)
- MCS identifies the shared pharmacophore core (27/37 atoms = 73% coverage)
- The shared atoms get coordinates directly from the crystal structure (highest fidelity)
- The non-shared atoms (4 atoms in the warhead not present in 88Z) are placed by
  ETKDG conformer generation followed by MMFF94 force-field minimisation

---

## Analysis parameters

### poi_min / poi_max (per-POI residue ranges)

**What they control:** Which residues in the combined ternary complex PDB are attributed to
the POI (kinase) for interface analysis, as opposed to VHL or elongins.

**Why non-trivial:** pRosettaC outputs a combined PDB where:
- Chain A residues 1–~160 = VHL
- Chain B = ElonginB
- Chain C = ElonginC
- Chain A residues starting from the POI's original numbering = kinase

For kinases with low start residue numbers (e.g., SLK starts at residue 21, MAP4K5 at
residue 1), the kinase numbering overlaps with VHL numbering (residues 1–103). The analysis
script uses `poi_min` to exclude these shared-numbering residues.

**Standard rule applied:** If the first kinase residue in the structure is ≤ 103, set
`poi_min = 104` to avoid including VHL residues. This affects: SLK, MAP4K5, TNIK,
PIP4K2C, CDK17, MAPK14.

### E2 lysine distance cutoffs

- `n_lys_within_25A_vhl`: counts POI lysines within 25 Å of VHL surface — a very generous
  filter that captures all potentially relevant lysines
- `n_lys_e2_optimal`: counts POI lysines within **35–55 Å** of the modelled E2 active
  cysteine position. This is the mechanistically relevant range:
  - < 35 Å: too close (steric clash with the ubiquitin chain or E2 enzyme)
  - 35–55 Å: optimal range for E2 → lysine ubiquitin transfer
  - > 55 Å: too far for productive ubiquitination

These cutoffs are based on:
- Crystal structures of E2·ubiquitin·substrate complexes
- The structural biology of ubiquitin chain extension (Ub is ~4 nm extended)
- Published pRosettaC analysis criteria (Zaidman et al., 2020)

### Interface contact distance: 5 Å

**Why 5 Å:**
- Standard cut-off for "residue in contact" at protein-protein interfaces
- Captures direct van der Waals contacts without including second-shell residues
- Directly comparable to other published PROTAC interface analyses

---

## Score shifting

**Why score.sc values need shifting:**
The original pRosettaC `clustering.py` hardcodes `if score < 0` as the filter for valid
structures. PyRosetta's `PackRotamersMover` (used in constraint_generation) produces
positive Rosetta energy scores because it starts from an unrelaxed structure and the
absolute scale differs from a full-system `ref2015` score.

**Fix:** `resume_pipeline.py::shift_scores_if_needed()` shifts all scores by
`-(max_score + 100)` so the maximum score becomes -100 REU. This makes all structures
pass the clustering filter while preserving relative score differences. The original
scores are backed up in `score.sc.orig` and the shift value is recorded in `score.sc.shift`
so downstream analysis can recover original values.

**This affects the score values reported in cluster_features.csv** — they are shifted
scores (always negative). The `score_shift` field in `run_poi_full_analysis.py`
`POI_CONFIGS` records the shift applied, allowing recovery of the true Rosetta score.

---

## Machine-level parameters

### N_WORKERS = 4 → 18

- Initial runs used N_WORKERS=4 (17% of 24 CPUs) to avoid overloading the system
- After the first POI (MET) completed successfully, switched to N_WORKERS=18 (75% of 24
  CPUs) for all subsequent runs
- 75% CPU leaves enough headroom for the OS and monitoring processes
- N_WORKERS controls the **local Rosetta docking** parallelism
- N_WORKERS_RESUME controls the **constraint_generation** parallelism

### N_WORKERS_RESUME = 18 (defaulting to cpu_count - 2 = 22)

- constraint_generation.py spawns one Python subprocess per pose
- PyRosetta has high per-process startup cost (~3–5 seconds) so more workers ≠ linearly
  faster, but 18–22 workers provides good throughput
- Setting N_WORKERS_RESUME=18 matches the docking worker count for consistent resource use
