# Metric names for the existing SJF8240 experiment

| Quantity | Defensible description | Do not infer without additional evidence |
|---|---|---|
| Historical “VHL”–ligand distance | Elongin-B–ligand minimum distance for the 141 audited MD systems | VHL-pocket retention or VHL unbinding |
| Historical “bridging <5 Å” | Same-frame ligand proximity to the legacy-named E3 component and POI groups | Productive ternary-complex probability or degradation |
| Corrected true-VHL/POI contact occupancy | Sampled simultaneous proximity to sequence-identified VHL and target, with stated atom selection and cutoff | Correct native head/pocket engagement or equilibrium binding probability |
| Head/pocket contact retention | Contacts between explicitly mapped head atoms and defined pocket residues, compared with a bound reference | Affinity, cooperativity, ubiquitination or cellular degradation by itself |
| GNN marginal-frequency product | Product of separate sampled contact frequencies | Observed same-frame co-contact frequency |
| E2 reachability descriptor | Geometry/proximity to modeled E2 positions under documented frame transformations | Catalytic ubiquitination probability |
| Multiplicative Pdeg-style score | Composite of explicitly defined sampling/proximity/reachability descriptors | Calibrated degradation probability without a validated probabilistic model |
| Saved GNN score | Historical trained-model output on its recorded feature/label cohort | A calibrated cellular probability or validation of the features' mechanistic labels |

The existing corrected E2 script (`15_e2_reachability_corrected.py`) explicitly anchors alignment to sequence-verified Elongin C (chain B, 85 Cα atoms). That is consistent with the compiled-topology role audit. Do not blanket-label that coordinate transformation wrong because another script mislabeled chain A. However, its bridging-conditioned subsets still require checking against the exact contact definition it reads.

Preserve original structures, trajectories, metrics and model outputs. Correct descriptions and audit derivatives separately; re-analysis does not imply rerunning PRosettaC or MD.
