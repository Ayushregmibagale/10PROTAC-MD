# Bai-style productive-fraction test

Read FINDINGS.md or FINDINGS.docx. Completed VHL adaptation of the16Å exposed-lysine/clash productive-fraction rule. Extra orientation and bridging requirements are removed. This is not a literal CRBN/CRL4A or full Rosetta-ensemble reproduction.

MD ensemble AUROC0.708;clustered retained starting-pose fraction AUROC0.583. The earlier adapted CRL2 score was0.750. These are retrospective eight-target scores; no new GNN training or simulations. Native catalytic-reference,missing endpoint,retained-ensemble and historical molecular-identity limitations remain unresolved.

Run numbered scripts in order using /home/ayushregmibagale/miniconda3/bin/python3. Dependencies and large source packets are supplied by sibling ../crl2_ubiquitination_v1 plus its original read-only MD sources.

1. scripts/01_score_productive_fraction.py — score the16Å rule and declared sensitivities from frozen141-pose packets.
2. scripts/02_cluster_retained_poses.py — reconstruct frame0 for141 retained MD poses,align trueVHL,cluster targetCA RMSD at3Å,and score cluster medoids. Retained member weights and equal-cluster scores are both saved. This uses complete-linkage/medoids,not the paper's MMTSB/interface-energy representatives.
3. scripts/03_verify_and_report.py — independently replay all Boolean classifications and metric arithmetic;write findings Markdown/Word.
4. scripts/04_check_package.py — check source integrity,cohort completeness and output receipts.

No target is relabelled. SIK3/RPS6KB2 remain unknown and excluded. All configurations stay in the denominator when a genuine target–scaffold clash makes them unproductive. Numerical PASS is not experimental validation. The full original Rosetta ensembles and catalytic endpoint calibration would be needed for a closer structural reproduction.
