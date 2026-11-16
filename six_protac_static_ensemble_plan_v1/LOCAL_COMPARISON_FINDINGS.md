# Matched nominal docking-budget comparison — 16 November 2026

Both completed local runs passed source/input/output integrity verification. Full raw outputs are archived losslessly outside Git with a verified gzip stream and SHA256 receipt. No MD or new Ernie job was run.

| Reference | Earlier50global x10local | New500global x1local |
|---|---|---|
| GSK215 | 3/8retained pass; best1.5628A; top1passes | 0/5pass; best5.6865A; top1fails |
| MZ1 external control | 0/187pass; best8.3451A | 0/201pass; best13.2798A |

Each scheme nominally attempted500local poses with up to20conformers each. Early feasibility rejection means these are not equal realized conformer totals. Different machines and derived task seeds prevent a paired runtime interpretation. One seed is insufficient to establish a generally optimal budget. Do not choose a production schedule from the better GSK result alone.

## Where recovery is lost

All500saved global/local orientation pairs per reference were evaluated retrospectively, without native-based selection. GSK215 best initial orientation3.4782A, best after local refinement5.6865A; the one initial orientation with RMSD<=5A was lost. MZ1 best initial5.5877A, best local10.9900A. Median local-minus-global RMSD was+1.5231A forGSK and+0.9861A forMZ1. These RMSD-only measurements do not certify the conjunction contact gate, linker feasibility or steric quality.

The trace implicates refinement as a stage requiring investigation, while MZ1 also lacks an initial orientation within5A. It does not prove a single cause. More global orientations alone did not repair recovery. Next inspect the actual Rosetta representation and search/refinement protocol in a separately recorded diagnostic, preserving the immutable baseline and all failures.

Independent chemistry/steric review reconstructs the expected molecular graph/stereo for all206retained poses.57/201MZ1 poses exceed their native ligand–protein maximum overlap by more than0.5A;0/5GSK do. This is an exploratory continuous comparison, not a newly validated acceptance threshold. MechanicalQC is insufficient structural qualification.

Full kinome production remains on hold pending recovery/pocket-fit qualification. A separately qualified bounded technical pilot can precede production, with explicit scope and final submission approval.
