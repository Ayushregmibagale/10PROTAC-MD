# Kinome geometry repairs — 14 November 2026

Repairs and verification are complete for the bounded ten-pose diagnostic and VHL references. This does not complete a validated SJF8240 kinome reranking.

## Completed repairs

- Preserved every original docking PDB, trajectory, label and analysis output.
- Performed two recorded local repair stages for the original ten leading poses: MMFF94s plus explicit positional restraints and overlap penalties, fixed proteins, deterministic seed42 and at most two starting perturbations per stage. No docking or MD production was run. The second stage increased steric penalties without changing the evaluation gates. Its finite-difference gradient checks passed. Optimization energies are repair objectives, not binding energies or cross-target ranking scores.
- Regenerated 198 VHL reference heads by chemical atom maps from deposited4W9H/3JF bound-parent coordinates and fitted them through sequence-verified ELOC. All198 pass bond-geometry checks. They are separate derivatives, not overwrites of the old PT0 files.
- Attempted sequence and CE structural transfer of a chemically correct foretinib head from3LQ8 into all198 shortlisted archived kinase models. All198 failed the existing domain/pocket fit requirements. No threshold was relaxed to invent a successful placement. The same alignment code passes an identical-structure positive control with near-zero pocket RMSD.
- Saved chemically correct native MET/foretinib and VHL parent references as a separate control bundle. These are independently bound parents, not an observed ternary complex.

## Before/after result

The repair assessment requires maximum intercomponent and nonbonded ligand overlap≤0.6Å; the internal-ligand check is stricter than the earlier0.8Å diagnostic limit. It also requires sensible bonded distances1.0–2.1Å, preserved archived graph/3D stereo, optimizer convergence and at least three atoms from each chemical head within4Å. Fixed protein–protein clashes are included separately, so a ligand-only repair does not incorrectly rescue an overlapping protein arrangement.

3 historical poses pass these operational geometry/contact checks: FER, TGFBR1, PRKG1. They are useful only for studying the archived molecular model. They are not validated binders or supplier-matched SJF8240 candidates, and their pocket engagement remains unverified.

| Kinase | Ligand overlap after repair, Å | Internal overlap after repair, Å | Fixed protein overlap, Å | Smaller head contact count at4Å | Geometry/contact result |
|---|---:|---:|---:|---:|---|
| TXK | 0.527 | 0.362 | -0.192 | 1 | Fail |
| MAP3K14 | 0.401 | 0.357 | 0.626 | 5 | Fail |
| FER | 0.432 | 0.348 | 0.517 | 5 | Pass for archived-model diagnostics only |
| EPHA1 | 0.419 | 0.350 | 0.513 | 3 | Fail |
| RPS6KB1 | 0.405 | 0.369 | 1.175 | 3 | Fail |
| PDIK1L | 0.405 | 0.347 | 1.049 | 12 | Fail |
| TGFBR1 | 0.381 | 0.354 | 0.170 | 4 | Pass for archived-model diagnostics only |
| MAP3K13 | 0.404 | 0.350 | 0.417 | 0 | Fail |
| PRKG1 | 0.406 | 0.352 | 0.252 | 10 | Pass for archived-model diagnostics only |
| PIK3CG | 0.432 | 0.375 | 1.341 | 0 | Fail |

EPHA1 did not meet the optimizer convergence criterion in the stronger repair stage, despite its reduced overlaps. TXK remains strained: one covalent bond is approximately2.107Å, above the unchanged2.1Å bound, and it lacks enough4Å contacts. MAP3K14 retains a protein–protein overlap above0.6Å. RPS6KB1, PDIK1L and PIK3CG retain substantial fixed protein–protein overlaps. MAP3K13 loses correct-head proximity. These failures are structural outcomes, not non-binding or non-degradation labels. Alternative cutoffs4.5Å and5Å are recorded in REPAIR_COMPARISON.tsv; they do not fix molecular provenance or pocket-validation failures.

## Molecular identity — do not relabel existing results

The saved historical ligand graph has InChIKey `MCJFROJTDORRBB-UHFFFAOYSA-N` without stereo; the historical intended SMILES is `MCJFROJTDORRBB-ZWUDBOSPSA-N`. The actual archived3D coordinates in these ten poses imply `MCJFROJTDORRBB-OVEQSRFISA-N`. Repairs preserve that actual3D graph/stereo. The supplier-verified SJF8240 is `SCUZPTNDZHEHLB-ZWUDBOSPSA-N`. These have different chemical connectivity, not merely different conformations. A geometry minimization cannot convert one into the other.

The historical reference-preparation code assumes that raw PDB coordinate ordering matches an independently parsed reference SMILES ordering, then overwrites part of an independently embedded molecule. That is an unsafe atom-correspondence procedure and provides a concrete mechanism for the head geometry failures. The new reference route uses CCD atom IDs and graph checks. The exact full historical execution chain has not been reconstructed, so the code finding is not claimed as the sole proven cause.

## What remains before an actual SJF8240 top ten

The VHL reference geometry and bounded local repair scripts are fixed and checked. A defensible supplier-matched top ten still requires suitable target-pocket structures/conformations and chemically correct ternary construction. Existing archived kinase models cannot be certified through the tested native-head transfer. That work is a new candidate-generation branch, not a repair of completed MD labels or a simple rerun of the old contact ranking. Do not launch new MD from these historical repaired poses as though they were supplier SJF8240.

## Files

- REPAIR_COMPARISON.tsv: all before/after results and three cutoff checks.
- stronger_relaxation/legacy_repairs/<gene>/: repaired ligand SDF, observed-fragment PDB, complete reconstructed-E3 PDB, atom maps and first-party receipts.
- legacy_repairs/ and probe_archive/: preserved weaker-stage repairs and initial probe.
- supplier_reference_candidates/:198 repaired VHL heads and explicit kinase-placement failure receipts.
- supplier_native_MET_control/: correct chemical/native-parent control inputs; production_ready remains false.
- chemical_identity_audit.json, reference_rebuild_summary.tsv, independent_repair_verification.json: provenance, failures and checks.
- scripts/: complete source for identity auditing, reference rebuilding, local repair, stronger refinement and verification/reporting.

Numerical replay, file integrity, actual coordinate stereochemistry, fixed protein coordinates, all198 rebuilt VHL head geometries and a reference-alignment positive control pass. This is software/geometry validation, not prospective biological validation.
