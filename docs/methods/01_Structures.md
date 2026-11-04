# Step 01 — Structure Acquisition

## Overview

We use a tiered strategy combining experimental crystal structures (KLIFS) and
AlphaFold2 predictions (kinome_scan), selected per-PROTAC based on required
DFG conformation.

---

## Source A — KLIFS Crystal Structures

**Database:** KLIFS (Kinase–Ligand Interaction Fingerprints and Structures)  
**URL:** https://klifs.net  
**API:** OpenCADD-KLIFS Python package (`conda-forge opencadd=1.0.1`, Python 3.8)  
**Conda env:** `klifs_screen`

### Query

- Species filter: Human
- DFG state filter: `in`, `out`, `out-like`, `inter` (excludes `na`, `ab-loop`)
- Result: **13,969 structures** across **328 kinases** (KLIFS naming)

### Deduplication

Best structure per (KLIFS gene name, DFG state), sorted by:
1. `curated = True` first
2. `missing_residues` ascending
3. `resolution` ascending

Result: **497 structures** across **328 kinases** (KLIFS naming, ~400+ HGNC).

### Download

Script: `scripts/01_setup/00_download_klifs_structures.py --download`  
Method: `session.coordinates.to_pdb(structure_klifs_id, out_dir, entity='complex')`  
Output: `structures/{KLIFS_GENE}/HUMAN_{GENE}_{PDB}_{chain}_complex.pdb`

**Result:** 496/497 downloaded (RSK3 KLIFS ID 2777 / PDB 2z7r: KLIFS API failure;
recovered directly from RCSB PDB as fallback).

### Why KLIFS

- Experimental coordinates (crystal/cryo-EM), not predicted
- DFG conformation state annotated per structure
- Critical for foretinib (Type II, DFG-out): placing a DFG-out warhead on a
  DFG-in structure causes activation-loop clashes and incorrect ternary geometry
- Curated pocket alignment across 85 conserved residues

---

## Source B — AlphaFold2 Structures

**Source:** EBI AlphaFold DB v4, via UniProt REST API  
**Script:** `~/kinome_scan/scripts/02_download_af2_kinases.py`  
**Kinase list:** 641 reviewed human kinases from UniProt
  (`keyword:Kinase AND organism_id:9606 AND reviewed:true`)

**Result:** 631 / 641 AF2 structures downloaded to
`~/kinome_scan/data/structures/{UNIPROT_ID}.pdb`

### Why AF2 as supplement

- Full kinome coverage — ~200 kinases have no PDB crystal structure
- AF2 almost always predicts DFG-in (ground-state) conformation
- Appropriate for Type I inhibitor warheads (lapatinib, VS-4718, JH-I-25)
- Not used for foretinib-based PROTACs (DFG-out required, AF2 DFG-in incorrect)

---

## Tiered Selection per PROTAC

| PROTAC | Warhead | Type | DFG needed | KLIFS first | AF2 fallback |
|--------|---------|------|------------|-------------|--------------|
| SJF8240 | Foretinib | II | out / out-like | ✓ | ✗ |
| SJFδ | Foretinib | II | out / out-like | ✓ | ✗ |
| SJF1521 | Lapatinib | I | in | ✓ | ✓ |
| GSK215 | VS-4718 | I.5 | any | ✓ | ✓ |
| JNJ-1013 | JH-I-25 | I | in | ✓ | ✓ |
| MS143 | MK-2206 | Allosteric | — | — | — (excluded) |

---

## Merged Manifest

Script: `scripts/01_setup/00b_build_merged_manifest.py`  
Output: `structures/merged_manifest.tsv`

**Mapping key:** KLIFS UniProt IDs (`kinase.uniprot`) → AF2 filenames (`{UNIPROT_ID}.pdb`)

Columns per row (one row = one HGNC kinase):
- `hgnc_gene`, `uniprot_id`, `klifs_name`
- `klifs_in_path`, `klifs_out_path`, `klifs_out-like_path`, `klifs_inter_path` (+ resolution)
- `af2_path`
- `{PROTAC}_path`, `{PROTAC}_source` — recommended structure and source for each PROTAC

---

## Notes

- KLIFS uses non-standard gene names (e.g. `AurA` → AURKA, `AMPKa1` → PRKAA1).
  Mapping to HGNC done via `kinase.gene_name` field in KLIFS kinase metadata API.
- RSK3 (KLIFS ID 2777, PDB 2z7r): permanent KLIFS API failure; substituted with
  direct RCSB download (`curl https://files.rcsb.org/download/2z7r.pdb`).
- MS143 targets AKT via its PH domain (allosteric site, not ATP pocket).
  This binding mode does not generalise across the kinome; excluded from screen.
- JNJ-1013 reference co-crystal not yet confirmed; flagged as skip in
  `01_place_warhead_batch.py` until IRAK1 reference PDB is identified.
