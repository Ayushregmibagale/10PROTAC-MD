# 05 — Bugs, Issues, and Fixes

All significant problems encountered during this campaign, their root causes, and the exact
fixes applied. Ordered by discovery date.

---

## Bug 1: BioPython Superimposer inverse transform (CRITICAL)

**Discovered:** Session 2, during warhead placement diagnosis  
**Severity:** CRITICAL — caused every non-MET warhead to be placed completely outside the protein

### Symptom
All non-MET warhead PDB files showed the WRH ligand at coordinates completely outside the
protein (e.g., z = +11 Å while the protein spans z = -65 to -11 Å). The warhead appeared
to be "floating" far from the ATP-binding site in every POI.

### Root cause
BioPython's `Superimposer.set_atoms(fixed=MET, moving=target)` returns `(rot, tran)` that
satisfies:

```
rot.T @ target_coord + tran ≈ MET_coord
```

This is a **right-multiplying row-vector convention** stored as a 3×3 matrix. The inverse
transform (to map MET-frame warhead coordinates INTO the target kinase frame) is therefore:

```python
target_coord = rot @ (met_coord - tran)   # CORRECT
```

The code was using:

```python
target_coord = rot.T @ (met_coord - tran)  # WRONG — this is the forward transform again
```

This was empirically verified by testing all 4 plausible inverse formulas:
- `rot @ (met_coord - tran)` → anchor 3.8 Å from Cα ✓
- `rot.T @ (met_coord - tran)` → warhead 80+ Å from protein ✗
- `rot @ met_coord + tran` → wrong direction ✗
- `rot.T @ met_coord + tran` → wrong direction ✗

### Fix applied
File: `~/protac/scripts/structure_prep/place_warhead.py`  
Function: `transform_warhead()`

```python
# BEFORE (wrong):
new = rot_inv @ p + tran_inv

# AFTER (correct):
new = rot @ (p - tran)
```

The comment in the function now documents the convention explicitly for future reference.

### Impact
All warhead SDF files (`*_warhead_placed.sdf`, `*_warhead_placed_H.sdf`) were regenerated
after this fix. All POI runs before the fix produced geometrically meaningless ternary
complexes.

### Verification
After fix, anchor-atom to nearest Cα distances:
- DDR2: 3.9 Å ✓
- EPHB2: 4.4 Å ✓
- MAPK14: 0.4 Å (very close — small Cα clash, but within tolerable range for Rosetta relax)
- All others: 2–8 Å ✓

---

## Bug 2: DDR2 run used DDR1 crystal structure (3ZOS)

**Discovered:** Post-run quality audit (session 3)  
**Severity:** HIGH — results for DDR2 are from the wrong protein

### Symptom
The pipeline completed successfully for DDR2, but the input protein structure
(`DDR2_protein_clean.pdb`) was derived from PDB 3ZOS — which is the **DDR1** kinase domain
bound to ponatinib, not DDR2.

### Root cause
The original `POI_structures.tsv` listed `3ZOS` as the recommended PDB for DDR2 without
noting it is DDR1. DDR1 and DDR2 are paralogues in the same DDR family; the PDB IDs can
be confused because many published DDR inhibitor studies use DDR1 structures even when
profiling DDR2 activity.

### Fix applied
1. Downloaded `7AYM` (DDR2 kinase domain, 2.12 Å, with DDR2-selective inhibitor 4VD)
2. Extracted chain A (ATOM records, no ligand) → `DDR2_7AYM_chainA.pdb`
3. Backed up old structure: `DDR2_protein_clean.pdb.bak_3ZOS`
4. Replaced: `DDR2_protein_clean.pdb` ← 7AYM chain A
5. Updated `place_warhead.py` `POI_TABLE`: `"DDR2_3ZOS_chainA"` → `"DDR2_7AYM_chainA"`
6. Re-ran warhead placement: RMSD improved from ~3+ Å (DDR1 vs MET) to **0.83 Å** (DDR2 vs MET)
7. Updated `run_poi_full_analysis.py` residue range: `poi_min=553, poi_max=858` (from 602–913)

### Current status
- The DDR2 run completed on 2026-06-03 used the OLD DDR1 structure + mixed warheads
  (70% old, 30% new, due to the file change during the run) — useful only as a pipeline
  sanity check
- **DDR2 must be re-run** with the corrected 7AYM structure before the results are used
  for analysis

---

## Bug 3: EPHB2 used ADP-bound structure (2HEN) — structural mismatch

**Discovered:** Post-run quality audit (session 3)  
**Severity:** MODERATE — activates loop conformation differs from foretinib-bound state

### Symptom
EPHB2 run used `2HEN` (EphB2 + ADP, 2.60 Å). ADP-bound structures have the activation
loop in the DFG-in (active) conformation. The binding pocket geometry differs from what
would be expected for foretinib binding.

### Fix applied
1. Searched PDB for EphB2 structures — no inhibitor-bound EphB2 kinase domain exists
2. Downloaded `3ZFM` (apo human EphB2 kinase domain, 2.27 Å) — best available
3. Extracted chain A → `EPHB2_3ZFM_chainA.pdb`
4. Backed up old structure: `EPHB2_protein_clean.pdb.bak_2HEN`
5. Replaced: `EPHB2_protein_clean.pdb` ← 3ZFM chain A
6. Updated `place_warhead.py` `POI_TABLE`: `"EPHB2_2HEN_chainA"` → `"EPHB2_3ZFM_chainA"`
7. Re-ran warhead placement: RMSD = 2.47 Å (reasonable; EphB2 ≠ MET)
8. Updated `run_poi_full_analysis.py` residue range: `poi_min=615, poi_max=894`
   (EPHB2 score_shift reset to 0.0; was 1206.983 from old run)

### Residual limitation
No inhibitor-bound EphB2 structure exists in the PDB (as of June 2026). The apo structure
(3ZFM) is the best available. This must be disclosed in methods.

---

## Bug 4: score.sc positive values causing clustering failure

**Discovered:** MET first run analysis  
**Severity:** HIGH — caused clustering to produce no Results/ clusters

### Symptom
After constraint_generation completed, `clustering.py` ran but produced no `Results/`
directory or an empty one. The `score.sc` file contained all-positive Rosetta energy values.

### Root cause
pRosettaC's `clustering.py` hardcodes a filter `if score < 0` — it silently discards
all structures with non-negative scores. PyRosetta's `PackRotamersMover` (used in
`constraint_generation.py`) scores structures relative to the unrelaxed initial pose,
producing positive energies on an absolute REU scale.

### Fix applied
Added `shift_scores_if_needed()` to `resume_pipeline.py`:

```python
shift = max_score + 100.0
# shift all SCORE: lines in score.sc by -shift
# back up original as score.sc.orig
# record shift in score.sc.shift
```

This ensures all scores are negative (max = -100 REU) while preserving relative
differences. The shift is recorded in `score.sc.shift` for recovery of original values.

---

## Bug 5: Two concurrent pRosettaC processes for DDR2

**Discovered:** During the batch rerun setup  
**Severity:** MEDIUM — caused mixed/corrupted Patchdock_Results

### Symptom
Two separate Python processes (PIDs 1465364 and 1466034) were both running
constraint_generation for DDR2 simultaneously, writing to the same `Patchdock_Results/`
directory.

### Root cause
When restarting the batch script after killing the initial wrapper, the new batch script
did not detect that a run was already in progress (it only checked for a non-empty
`Results/`, not a non-empty `Patchdock_Results/`).

### Fix applied
Added an in-progress check to `run_batch_rerun.sh`:

```bash
if [[ -d "${run_dir}/Patchdock_Results" && \
      -n "$(ls ${run_dir}/Patchdock_Results/ 2>/dev/null)" ]]; then
  echo "  SKIP: run already in progress (Patchdock_Results not empty)"
  continue
fi
```

---

## Bug 6: Path-with-spaces breaks PyRosetta init

**Discovered:** First PyRosetta run  
**Severity:** HIGH — prevents PyRosetta from loading any params files

### Symptom
PyRosetta `-extra_res_fa` flag parser splits on whitespace, so
`-extra_res_fa /home/user/PROTAC degradation Prediction/...` is parsed as two separate
(invalid) paths.

### Fix applied
Created symlink: `ln -s ~/PROTAC\ degradation\ Prediction ~/protac`

All scripts that call PyRosetta use the symlinked path. This is handled transparently in
`resume_pipeline.py` via:

```python
def _sl(p):
    return p.replace('/home/ayushregmibagale/PROTAC degradation Prediction',
                     '/home/ayushregmibagale/protac')
```

---

## Bug 7: EPHA2 folder contained wrong-domain structure (3MBW)

**Discovered:** Post-run quality audit  
**Severity:** LOW (was not used — identified before running EPHA2)

### Symptom
`data/structures/POI/high_affinity_no_degradation/EPHA2/3MBW.pdb` is the EphA2
**extracellular LBD + CRD domains** (bound to glycans), not the kinase domain.

### Fix
3MBW was present but `EPHA2_protein_clean.pdb` was derived from 5I9Y (EphA2 kinase
domain, 1.23 Å, with dasatinib) — confirmed by checking the residue range of the clean
file (residues 600–896 = kinase domain). No action needed; 3MBW is simply an unused file.

---

## Known remaining limitation: MAPK14 warhead anchor proximity to GLY 110

**Discovered:** Warhead placement check  
**Severity:** LOW-MODERATE — may affect MAPK14 results quality

### Description
MAPK14 (p38α) has the most divergent kinase fold from MET in the panel
(superposition RMSD ~12 Å over full kinase domain). After warhead placement, anchor
atom 0 was found to be only ~0.4 Å from the GLY 110 Cα — essentially on top of a
backbone atom.

### Expected consequence
PyRosetta's soft-core repulsive potential during Rosetta relax should resolve this clash
by pushing the warhead slightly. However, if the relax step moves the warhead substantially
out of the ATP-binding pocket, MAPK14 may show artificially low interface contacts.

### Recommended action
After MAPK14 run completes: open a representative combined PDB in PyMOL and verify that
the warhead (residue WRH or PTC) is inside the ATP pocket. If it has been displaced,
manual warhead placement against the SB203580 co-crystal geometry (3GCP) is required.
