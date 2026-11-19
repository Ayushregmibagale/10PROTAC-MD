# Historical DFG preparation and SJF8240/MAPK14 control

The historical source is `~/kinome_scan/scripts/07_prepare_kinome_pRosettaC.py`. Its helpers find the last DFG (or, if absent, DLG) in a structure, take up to 120 residues on either side, pair C-alpha atoms by offset from the motif, and transform the supplied head into the target frame. It logs the fit RMSD and writes input files without a fit or steric rejection. The standard runner parameters in the legacy source hardcode anchors 6 and 1, chain sets ABC and A, and the historical SJF8240 SMILES.

## What the historical ledger establishes

The archive contains 461 successful preparations from 641 rows, with 164 missing a DFG/DLG, 10 missing a PDB, and six lacking enough residues. Of the 461 successful historical kinase IDs, 397 also occur in the current 435-kinase panel. Recorded successful fit RMSDs span 6.06–51.75 Å, with median 11.00 Å. The original SJF8240/MAPK14 placement had 10.85 Å C-alpha fit RMSD and 3.00 Å maximum atom-type-aware vdW overlap.

Independent checks measured contacts for 460 of 461 archived placements. Of those, 454 had a minimum heavy-atom distance below the current 1.5 Å gross-contact criterion, and 457 exceeded 0.6 Å maximum overlap. One had no contacts within 4.5 Å. These counts describe the historical placements and the stated geometric check. The historical SJF8240 SMILES InChIKey is `MCJFROJTDORRBB-ZWUDBOSPSA-N`; the current supplier-verified representation is `SCUZPTNDZHEHLB-ZWUDBOSPSA-N`. Reference substitution therefore needs identity and anchor verification.

## Direct test of the proposed six-reference substitution

The benchmark reused the original DFG functions and the six current pilot binary references. It generated 397 placements per compound, for 2,382 of 2,610 requested pairs. Thirty-eight target IDs were missing from this method's output for each compound, including 11 without the legacy AF2 file and 27 with no eligible motif window. Of the 2,382 placements, 381 passed the existing 1.5–4.5 Å contact screen. Forty-five also had worst vdW overlap at or below 0.6 Å, including 18 pairs absent from the latest prepared-candidate ledger. The 0.6 Å value is a steric review band. These 18 are geometry candidates. They have not been staged or designated release-ready.

## SJF8240/MAPK14

MAPK14/p38α is the campaign's primary SJF8240 control. Smith et al. reported enhanced p38α selectivity and activity against p38δ under their assay conditions. The MAPK14 reference protein comes from the deposited 1W82 structure. Its inhibitor was replaced in the model using the foretinib scaffold from 5IA4 and explicit KLIFS positions. This creates a modeled input with the supplier-verified SJF8240 connectivity. [Smith et al.](https://www.nature.com/articles/s41467-018-08027-7) · [PDB 1W82](https://www.rcsb.org/structure/1W82)

The whole-pocket fit gives 1.81 Å RMSD and 1.55 Å over the head neighborhood. The transferred binary initially had 1.64 Å worst overlap. A rigid-head refit and fixed-head protein side-chain repack reduced that to 0.71 Å, with 2.63 Å minimum heavy-atom distance. The established SJF8240 VHL head remains in the input, and the modeled foretinib head maps to the exact supplier ligand. A separate bounded 12-global × two-local test is running on Ernie.
