# SJF8240 kinome-screen analysis

This folder contains the current analysis scripts for the broad pRosettaC kinome screen and its existing MD follow-up. It is separate from `SJF8240_renanalysis2`, which contains the smaller-cohort MD/GNN/E2 workflow.

Read [SCRIPT_GUIDE.md](SCRIPT_GUIDE.md) before running. This is an analysis workflow, not a docking or MD-production launcher. No degradation labels are inferred from bridging alone.

The raw screen contains 461 accession directories, rather than exactly 400. Availability and pose counts are recorded in `results/inventory.tsv`; missing runs are not negative biological results.

Current outputs, logs and validation evidence are linked to the active review workspace so the running screen can finish without being interrupted. The scripts and configuration are physical copies in this folder. The source workspace is `/mnt/c/Users/Ayush Regmi Bagale/SJF8240_kinome_pose_review_20261113`; retain it while these links are used.

Latest completed comparison: [engagement_v2/FINDINGS.md](engagement_v2/FINDINGS.md). The earlier top ten fails stronger checks and should not be advanced unchanged. See stages 6–8 in the guide.

Latest: [geometry repair findings](repair_v3/FINDINGS.md). Local repairs preserve historical chemistry; three pass the diagnostic. Supplier-matched SJF8240 selection remains unresolved because archived kinase-head transfers fail the reviewed criteria.
