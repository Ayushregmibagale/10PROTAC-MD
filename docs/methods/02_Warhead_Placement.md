# Step 02 — Warhead Placement

## Overview

Each PROTAC's warhead is placed into every kinase structure frame by Cα
superposition onto a reference co-crystal, then transferring the warhead
coordinates with the inverse rotation.

---

## Foretinib (SJF8240 / SJFδ) — Reused from kinome_scan

**Source:** `~/kinome_scan/data/kinase_pRosettaC/{UNIPROT_ID}/`

The kinome_scan project (`07_prepare_kinome_pRosettaC.py`) already placed the
foretinib warhead into 461 kinase frames using:
- AlphaFold2 structures (all DFG-in predicted)
- DFG motif found by sequence (`rfind('DFG')` / `rfind('DLG')`)
- 240-residue kinase domain window extracted around DFG motif
- Cα superposition onto DDR2 reference (`DDR2_protein_clean.pdb`)
- Warhead transferred using inverse Kabsch rotation

These 461 `protein.pdb` + `warhead.sdf` pairs were symlinked into
`prosettac_inputs/SJF8240/` and `prosettac_inputs/SJF_delta/`.

**Note:** These use AF2 structures (DFG-in predicted), not crystal structures.
The 144 kinases with KLIFS DFG-out/out-like structures in the merged_manifest
are a higher-confidence subset. pRosettaC samples pose space regardless, so
AF2 DFG-in inputs remain valid — the DFG motif region is correct geometrically
even if not in the open conformation.

**Original warhead SDF:** `~/protac/data/structures/PROTAC/SJF8240/SJF8240_warhead.sdf`

---

## Lapatinib (SJF1521), VS-4718 (GSK215), JH-I-25 (JNJ-1013)

**Script:** `scripts/01_setup/01_place_warhead_batch.py --run`

### Warhead SDFs (3D conformers from PubChem)

| Warhead | PROTAC | PubChem CID | Source |
|---------|--------|-------------|--------|
| Lapatinib | SJF1521 | 208908 | PubChem 3D conformer |
| Defactinib (VS-4718) | GSK215 | 25151352 | PubChem 3D conformer |
| JH-I-25 | JNJ-1013 | 44449334 (CAS 1042673-20-0) | PubChem 3D conformer |

Saved to: `warheads/{lapatinib,VS4718,JHI25}_warhead.sdf`

### Reference co-crystal PDBs

| PROTAC | Reference PDB | Kinase | Chain |
|--------|--------------|--------|-------|
| SJF1521 | 1XKK | EGFR | A (289 Cα) |
| GSK215 | 5JC7 | FAK | A (641 Cα) |
| JNJ-1013 | 6EGA | IRAK1 | A (267 Cα) |

Saved to: `reference_structures/{1XKK,5JC7,6EGA}_kinase.pdb`

### Algorithm

For each (PROTAC, kinase) pair:
1. Load Cα atoms from kinase structure and reference PDB
2. Truncate to `min(N_kin, N_ref)` residues (N-terminal alignment)
3. Kabsch align kinase Cα onto reference Cα → rotation R, translation t
4. Apply inverse transform to warhead atoms: `new_pos = (warhead_pos - t) @ R`
5. Write `warhead.sdf`, `protein.pdb` (chain D, residues 200+), `Protac_params.txt`

Minimum 30 Cα atoms required; pairs with fewer are skipped.

### Results

| PROTAC | Total | Placed | Failed |
|--------|-------|--------|--------|
| SJF1521 | 645 | 644 | 1 |
| GSK215 | 651 | 650 | 1 |
| JNJ-1013 | 645 | 644 | 1 |

**Failed:** PRKCH — too few Cα atoms in AF2 model (stub prediction).

### Limitations

- N-terminal Cα truncation is a simple approximation; MCS-based alignment
  would be more accurate for divergent kinase folds but requires RDKit MCS
  which is slow at this scale.
- PubChem 3D conformers are single low-energy conformers, not necessarily
  the bioactive binding conformation. For lapatinib and VS-4718 this is
  acceptable since the warhead is flexible and pRosettaC samples linker space.
- FAK reference (5JC7) is a large multi-domain kinase (641 Cα); the FERM
  domain Cα will dominate alignment. This may reduce accuracy for non-FAK kinases.

---

## Foretinib warhead origin note

SJF8240_warhead.sdf was originally prepared as the foretinib moiety of SJF8240
(the c-MET/p38α degrader), extracted and 3D-minimised from the full PROTAC
structure. It is identical to the kinase-binding head used in SJFδ since both
PROTACs share the foretinib warhead with different linker lengths.
