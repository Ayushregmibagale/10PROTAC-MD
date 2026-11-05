# Bound-head preparation and coordinate transfer

Use the exact Tocris full molecular graph recorded in `config/molecules.json`.
Do not substitute a parent inhibitor or reuse SJF8240 SMILES across degraders.
The inherited file named VS4718 contains pexidartinib; it is not a valid GSK215 head.

GSK215 is extracted from 7PI4 using the 7QB CCD bonds and deposited coordinates.
Its reconstructed stereochemical InChIKey matches Tocris catalog 7818 exactly.
The explicit head split cuts CCD bonds N4-C20 and C21-N5, excludes CH2-C(=O)
from both heads and preserves all rings. The atom mapping and one-based anchors
are recorded in `reference_structures/validated/GSK215/heads.json`.
The split remains a candidate until verified against the exact pRosettaC engine.

Protein alignment uses BLOSUM62 sequence correspondence or explicit CE structural
correspondence across divergent kinases, followed by iterative
core fitting, with separate pocket coverage and atom RMSD diagnostics. Ligand
coordinates start in the reference protein frame. The inverse row-vector transform
is `(reference_head - t) @ R` for `target @ R.T + t = reference`.
No N-terminal truncation or unfitted PubChem conformer transfer remains in this path.

Current conservative candidate checks: 100-400 domain C-alpha atoms, at least
100 matched residues, sequence identity >=25%, fitted core >=60% of aligned
residues, core RMSD <=2.5 A, mapped pocket coverage >=80% with >=10 residues,
pocket RMSD <=2 A, minimum head/protein heavy distance 1.5-4.5 A.
These are provisional preparation filters, not benchmark-certified accuracy cutoffs.

A native GSK215 self-placement control verifies coordinate transfer and identity.
It does not demonstrate independent docking recovery or cross-target generalization.
SJF1521, JNJ1013 and SJF_delta have mapped modeled conjugated heads; their
scientific review remains explicit before release. New candidates use `prepared_inputs/`; legacy files are preserved. SJF8240 now has
exact bound foretinib and unmethylated VHL heads extracted from 3LQ8/88Z and
4W9H/3JF, respectively, with a separate MET self-placement preparation control.

The supplementary SJFδ/MAPK13 control follows the template pair described by
Smith et al. (doi:10.1038/s41467-018-08027-7): foretinib from EPHA2 5IA4 is
overlaid on p38δ 4EYJ. Explicit KLIFS pocket positions 7189/2949 give 83 common
observed Cα atoms, RMSD 1.608 Å and minimum heavy contact 1.794 Å. Neither structure
is an experimental p38δ/foretinib complex. The model does not reproduce the
paper's OPLS3/Jaguar torsion treatment and is not an affinity prediction. Its
original-template evidence is retained in the preparation receipt; identity
self-placement RMSD is not substituted for the cross-template modeling error.
The original 112-candidate CE screen remains unchanged, so its limited coverage
must be reported and addressed before broad kinome claims.

A separate explicit-KLIFS-pocket comparison uses 346 cached official maps.
Observed positions are mapped to canonical UniProt residues using verified source
coordinates; AF2 inherits maps only from consistent evidence for the exact same
domain. 556 domains have such maps; 214 remain explicitly unavailable. Bijective
pocket fitting retains the numerical core, coverage, RMSD and clash limits above.
It produces 161 candidates, including all five known-target controls. Forty-one
overlap the prior 112 CE candidates. This comparison is not silently adopted,
pooled with the old screen or interpreted as validation of biological binding.
