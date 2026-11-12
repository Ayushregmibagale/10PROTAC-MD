# Dynamic MD method comparison

Both the legacy per-lysine table and the mechanistic_v3 lysine analysis use MD trajectories, not merely static docking poses. The per-lysine claim check acknowledges that dynamic origin.

The source20_per_lysine_e2_distance.py first calculates each frame/lysine minimum distance over E2 ensemble points, then reach_50A is the fraction of bridge or nonbridge FRAMES for which that minimum is below50 A. It is not the fraction of E2 ensemble positions within50 A. Its ensemble JSON concatenates13 identical arrays of51 points, giving663 rows but51 distinct coordinates. This duplication does not change nearest-point distances; it does not explain away the observed class pattern.

Legacy masking uses the historically named metrics.csv VHL distance that the independent topology audit identifies as ELOB. Its per-protein unwrapping does not apply the later ligand-connected coherent whole-component placement. The corrected method fits E2 points through Elongin C as well, but uses actual VHL heavy-contact bridging and explicit coherent molecular images. Both methods measure modeled geometry, not a ubiquitin-transfer reaction.

mechanistic_v3 retains41 frames per pose at0.5 ns spacing for its per-lysine SASA/E2 descriptors, compared with the old per-lysine table's full sampled trajectory. MD contact baseline graphs use201 frames. The supplemental delta50 check uses those existing41-frame packets and equally averages paired within-pose contrasts. It is not a full201-frame matched recomputation and should not be presented as a definitive rejection of the legacy50-A signal.

The latest0.750 matched-pair VHL GNN did not directly consume per-lysine delta_reach_50A. Its engineered inputs were residue identity, NZ SASA/burial, per-lysine E2 distance,80-A point reach fraction and joint exposed/reachable/bridge occupancy across4 windows. It therefore tests a combined representation, not the exact specific50-A delta hypothesis.

The independent labels are on8 proteins (6:2), not132 lysines or thousands of frames. The legacy lysine-level p~3.03e-9 can be reproduced but is not an independent-protein significance test. The mean per-protein legacy delta descriptively separates the two known negatives; exact two-sided protein-level rank p=0.07143. The fraction-positive statistic overlaps (AXL2/16 exceeds MET2/19 and RIPK2 1/13), with protein-level rank p=0.28571. These retrospective comparisons do not establish classifier generalization, and dynamic origin does not by itself make the lysines independent.
