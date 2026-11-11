# Weighted MD GNN plus E2 geometry

Option2 adds two graph-level descriptors to the weighted MD GNN: audited ligand-imaged mean minimum E2/lysine distance and paired true-VHL heavy-contact delta reach80. All eight corrected-label targets have finite values. GAT layers and all other training settings match the weighted MD-only version; the first classifier layer expands from128 to130inputs. Normalize using only the seven training targets in each fold. No labels or held-out feature values enter the fitted scaler.

All141 core E2 cases and28,341 frames are reproduced against original contact/summary hashes and the reviewed target tables. Both graph roles receive the same true-VHL-based scalar values. The old ABL1 -0.0466 delta is not reused; corrected ligand-imaged delta is -0.036848. Delta80 measures a modeled E2-point reach fraction, not burial of18 lysines. Source features remain candidate geometric proxies, adopted explicitly for exploratory GNN fitting at the user's request.

SIK3/RPS6KB2 lack corrected audited E2 descriptors in this source packet and are not predicted here; missing descriptors are not fabricated or replaced with zeros. No new docking/MD or conventional model is run. Epochs300/seed42 and fold-local class balance remain fixed before fitting.
