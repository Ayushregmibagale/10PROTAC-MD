# Ernie DDR2 Run Troubleshooting — 2026-06-22

## Summary of session

Attempting to re-run pRosettaC DDR2 ternary complex modeling on ernie (s1638308@ernie.bio.ed.ac.uk) with the Kabsch alignment fix for PROTAC placement.

**Final outcome:** All three bugs resolved. DDR2 run active (PID 3521046) with PROTAC correctly placed. Confirmed PASS and MARGINAL results across 459 checked poses (~1.5% feasibility rate, expected for pRosettaC).

---

## Bug 1: PROTAC always at origin (CRITICAL — pre-existing)

**Symptom:** `check_protac.py` on all combined_*_0001.pdb files shows PROTAC centroid at [0.50, 0.29, 0.13] — 69.9 Å from VHL, 73.4 Å from DDR2. FAIL.

**Root cause:** `gen_const_conf_free` in `pyrosetta_adapter.py` generates ETKDGv3 conformers at the origin and never aligns them to the docked binding site. The `docked_X_Y.pdb` file contains the correct warhead target positions but was only used for virtual atom center-of-mass, not rigid-body alignment.

**Fix applied (pyrosetta_adapter.py):** Added Kabsch rigid-body alignment block in `gen_const_conf_free` (approx. lines 1001-1035). After ETKDGv3 embedding, substructure-match HeadA_ha and HeadB_ha in the full PROTAC mol, get target coordinates from `docked_heads`, compute Kabsch R,t, apply to all PROTAC atoms before writing SDF.

**Verified working in isolation on ernie:**
- `_can_align: True`
- Before alignment: [-0.82, -0.09, 0.02]
- After alignment: [5.76, 31.7, 106.07]
- Target: [5.97, 31.25, 106.2]  ✓

---

## Bug 2: Stale confs_*.sdf files masking the fix

**Symptom:** After applying the Kabsch fix and restarting the pipeline, `confs_100_1.sdf` still had centroid at [0.5, 0.29, 0.13]. The fix was confirmed present (`grep -c "_can_align"` = 2) but not applied.

**Root cause:** When restarting, we only deleted `combined_*_0001.pdb` but left the old `confs_*.sdf` files (generated before the fix). The pipeline found these stale confs files and reused them, skipping `gen_const_conf_free` entirely.

**Fix:** Always delete BOTH `confs_*.sdf` AND `combined_*_0001.pdb` before restarting. Use `find` instead of `rm` glob to ensure all files are caught:
```bash
find ~/protac/pRosettaC/runs/high_affinity_degraders/DDR2/Patchdock_Results/ -name "confs_*.sdf" -delete
find ~/protac/pRosettaC/runs/high_affinity_degraders/DDR2/Patchdock_Results/ -name "combined_*_0001.pdb" -delete
```

---

## Bug 3: 0-byte confs_*.sdf + Kabsch in wrong function (RESOLVED 2026-06-22)

### Part A: Kabsch fix was in the wrong function

**Root cause discovered:** `_pl.GenConstConf` in `apply_patches()` (line 1126 of `pyrosetta_adapter.py`) was mapped to `gen_const_conf_anchored`, NOT `gen_const_conf_free`. The Kabsch alignment fix from Bug 1 was added to `gen_const_conf_free` but the pipeline never called that function.

The pycache on ernie had been compiled from an older version of the code where `_pl.GenConstConf = gen_const_conf_free` (or similar), so old runs appeared to use ETKDGv3. Once pycache was cleared, the current code ran `gen_const_conf_anchored` for the first time.

### Part B: gen_const_conf_anchored DG embedding fails for DDR2 poses

**Symptom after pycache clear:** All `confs_*.sdf` files were 0 bytes. `gen_const_conf_anchored` uses distance-geometry (DG) embedding with a crystal-structure `coordMap`. For DDR2 docking poses, the DG embedding fails for all 500 seeds (protein-protein geometry is geometrically infeasible for this anchor approach). No conformers are accepted → empty SDF.

### Fix applied (2026-06-22)

Changed the function mapping in `apply_patches()` from `gen_const_conf_anchored` to `gen_const_conf_free`:

```python
# OLD (line 1126):
_pl.GenConstConf = gen_const_conf_anchored   # DG embedding — fails for DDR2 poses

# NEW:
_pl.GenConstConf = gen_const_conf_free        # ETKDGv3 + Kabsch alignment to docked site
```

Also added Kabsch alignment to `gen_const_conf_anchored` as well (lines 935-967) as a safety measure, though this function is no longer the active code path.

### Verification

After fix + clean restart (PID 3521046, 16 workers):
- `confs_100_10.sdf` centroid: **[11.76, 44.08, 91.38]** (NOT at origin) ✓
- `combined_101_22_0001.pdb` check_protac result:
  - Closest PROTAC atom to VHL: **6.8 Å** (MARGINAL)
  - Closest PROTAC atom to DDR2: **7.0 Å** (MARGINAL)
- Batch check across 459 combined PDBs: **PASS: 2, MARGINAL: 5, FAIL: 452** (~1.5% feasibility)

**~1.5% feasibility is expected and correct** — most PatchDock poses have VHL and DDR2 too far apart for the PROTAC to bridge. With 10,000 total poses, ~150 PASS/MARGINAL structures are expected, which is sufficient for clustering.

---

## Correct full restart procedure (use this going forward)

```bash
# Kill any running run
kill $(pgrep -f resume_pipeline) 2>/dev/null

# Delete ALL stale outputs (use find, not rm glob, to be safe)
find ~/protac/pRosettaC/runs/high_affinity_degraders/<POI>/Patchdock_Results/ -name "confs_*.sdf" -delete
find ~/protac/pRosettaC/runs/high_affinity_degraders/<POI>/Patchdock_Results/ -name "combined_*_0001.pdb" -delete

# Clear pycache to ensure updated pyrosetta_adapter.py is used
find ~/protac/scripts/rosetta_utils/__pycache__ -name "pyrosetta_adapter*.pyc" -delete 2>/dev/null

# Restart with 16 workers
cd ~/protac/pRosettaC/runs/high_affinity_degraders/<POI>
N_WORKERS_RESUME=16 python ~/protac/scripts/rosetta_utils/resume_pipeline.py 2>&1 | tee ~/poi_run.log &
echo "PID: $!"
```

**Verify after ~3-5 min:**
```bash
# Check confs centroid (should NOT be near [0, 0, 0])
python3 -c "
from rdkit import Chem; import numpy as np, glob
sdfs = [s for s in sorted(glob.glob('Patchdock_Results/confs_*.sdf')) if __import__('os').path.getsize(s) > 0]
m = Chem.SDMolSupplier(sdfs[0])[0]
print(sdfs[0], 'centroid:', m.GetConformer().GetPositions().mean(axis=0).round(2))
"

# Check PROTAC placement in a combined PDB
python3 ~/check_protac.py $(ls Patchdock_Results/combined_*_0001.pdb | head -1)

# Batch feasibility check
pass=0; marginal=0; fail=0
for pdb in Patchdock_Results/combined_*_0001.pdb; do
  result=$(python3 ~/check_protac.py "$pdb" 2>/dev/null | grep -E "PASS|MARGINAL|FAIL")
  case "$result" in
    *PASS*) ((pass++)) ;;
    *MARGINAL*) ((marginal++)) ;;
    *FAIL*) ((fail++)) ;;
  esac
done
echo "PASS: $pass  MARGINAL: $marginal  FAIL: $fail"
# Expected: ~1-2% PASS+MARGINAL for most POIs
```

---

## Docking file naming convention (ernie)

PatchDock output files: `pd.{global}_docking_{local:04d}.pdb` (e.g., `pd.100_docking_0001.pdb`)

`resume_pipeline.py` glob: `*_docking_????.pdb`

Suffix extraction: `pd.100_docking_0001.pdb` → suffix `100_1`

`constraint_generation.py` calling convention (from `resume_pipeline.py`):
```
python constraint_generation.py <VHL_SDF> <POI_SDF> <PROTAC_SMI> <suffix> <pdb_filename> <chains>
# e.g.: constraint_generation.py SJF8240_vhll_placed_H.sdf DDR2_warhead_placed_H.sdf SJF8240.smi 101_1 pd.101_docking_0001.pdb ABCA
# Must be run from Patchdock_Results/ directory
```

---

## check_protac.py

Utility script at `~/check_protac.py` (both home and ernie). Takes a combined PDB path and reports:
- PROTAC centroid vs VHL (chain A res ≤ 103) and POI (chain A res ≥ 553)
- Closest PROTAC atom distance to each protein
- PASS (<5 Å to both), MARGINAL (5-10 Å to both), or FAIL

---

## Key file locations (ernie)

| File | Path |
|------|------|
| pyrosetta_adapter.py | `~/protac/scripts/rosetta_utils/pyrosetta_adapter.py` |
| resume_pipeline.py | `~/protac/scripts/rosetta_utils/resume_pipeline.py` |
| constraint_generation.py | `~/protac/scripts/rosetta_utils/constraint_generation.py` |
| DDR2 run dir | `~/protac/pRosettaC/runs/high_affinity_degraders/DDR2/` |
| Patchdock_Results | `~/protac/pRosettaC/runs/high_affinity_degraders/DDR2/Patchdock_Results/` |
| VHL head SDF | `~/protac/data/structures/VHL/SJF8240_vhll_placed_H.sdf` |
| DDR2 head SDF | `~/protac/data/structures/POI/high_affinity_degraders/DDR2/DDR2_warhead_placed_H.sdf` |
| PROTAC SMILES | `~/protac/data/structures/PROTAC/SJF8240/SJF8240.smi` |
| Run log | `~/poi_run.log` |
