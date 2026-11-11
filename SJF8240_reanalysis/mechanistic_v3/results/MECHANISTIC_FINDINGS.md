# Per-lysine geometric findings

141 existing poses; 41 sampled frames each = 5,781 frames. Heavy-atom solvent occlusion, 1.4 A water probe, 960 sphere points. All original inputs unchanged. Time windows are 0–5, 5–10, 10–15 and 15–20 ns, last including 20. Each pose has equal weight. E2 quantities reproduce previously audited ligand-connected images (maximum absolute distance replay error recorded per pose).

These measurements do not support the earlier claim that bridging buries all 18 ABL1 lysines. That claim incorrectly interpreted an E2 reach-fraction change as burial. NZ solvent-accessible area now measures exposure directly. The geometries remain highly overlapping across the two labelled classes. Exposed means NZ SASA >5 square angstroms; reachable means minimum distance to a modeled E2 point <80 angstroms. Thresholds were fixed before inspecting per-target exposure results or fitting the enriched model and are broad geometric proxies, not a ubiquitination reaction criterion.

- ABL1: mean NZ SASA 35.81 A² (pose SD 1.33); exposed fraction in bridge frames 0.984; joint exposed/reachable/bridge occupancy 0.319.
- AXL: mean NZ SASA 33.39 A² (pose SD 1.07); exposed fraction in bridge frames 0.970; joint exposed/reachable/bridge occupancy 0.344.
- DDR2: mean NZ SASA 40.62 A² (pose SD 1.32); exposed fraction in bridge frames 0.986; joint exposed/reachable/bridge occupancy 0.283.
- EPHA2: mean NZ SASA 39.41 A² (pose SD 1.61); exposed fraction in bridge frames 0.978; joint exposed/reachable/bridge occupancy 0.282.
- EPHB2: mean NZ SASA 40.78 A² (pose SD 1.34); exposed fraction in bridge frames 0.993; joint exposed/reachable/bridge occupancy 0.278.
- MAPK14_full: mean NZ SASA 29.39 A² (pose SD 1.43); exposed fraction in bridge frames 0.906; joint exposed/reachable/bridge occupancy 0.422.
- MET: mean NZ SASA 38.51 A² (pose SD 1.56); exposed fraction in bridge frames 0.984; joint exposed/reachable/bridge occupancy 0.360.
- RIPK2_full: mean NZ SASA 31.38 A² (pose SD 1.09); exposed fraction in bridge frames 0.915; joint exposed/reachable/bridge occupancy 0.317.

The accompanying tables retain every lysine, pose variability, bridge-state coverage and four time windows. Bridge and nonbridge descriptive averages can have different contributing poses; their difference is not a paired causal estimate. No feature threshold or direction was selected to maximize class separation; no new conventional classifier was fitted.

SASA represents solvent access and steric occlusion by retained complex heavy atoms. It does not model approach of the full E2–ubiquitin complex, catalytic orientation, CRL2 motion or ubiquitin transfer. Hydrogen atoms are omitted from occluders consistently. C-alpha residue IDs are simulation numbering, not verified ubiquitination-site annotations. Only two independently labelled negatives remain, regardless of pose/frame count.
