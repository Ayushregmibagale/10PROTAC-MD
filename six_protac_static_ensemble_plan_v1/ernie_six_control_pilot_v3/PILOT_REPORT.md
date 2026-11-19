# Ernie six-control pRosettaC pilot v3

Approved scope is one static compound–kinase control per compound. Release SHA256 is `db2a56066b44a953ebf2013f39eb73426be78e3c279cece59ec59a6c3ae2064f`. The serial queue uses up to four workers, 50 global poses, two refinements per global pose, 20 conformers per feasible refinement and a 30-minute limit per case.

## SJF1521 EGFR input caveat

The frozen `preparation.json` records `production_ready: false`. The actual runner input supplies the one-based anchor indices `1 1`. The pilot preparation metadata does not document the atom identities, full-molecule correspondence or rationale for those selections. The pilot uses the approved frozen inputs and any resulting structures remain technical outputs.

## Raw results

All six approved controls reached a terminal outcome. Elapsed times are from the individual submission receipts. Retained counts are the frozen pilot's recorded mechanical QC counts.

| Control | Terminal result | Local PDBs | Linker logs | Scored | Retained | Elapsed |
|---|---|---:|---:|---:|---:|---:|
| GSK215 / PTK2 | Complete | 100 | 100 | 41 | 37 | 681.0 s |
| SJF8240 / MET | Complete, no feasible linker | 100 | 100 | 0 | 0 | 1,470.0 s |
| SJF1521 / EGFR | Complete | 100 | 100 | 58 | 55 | 925.6 s |
| JNJ1013 / IRAK1 | Failed during conformer processing | 100 | 91 | — | — | 549.0 s |
| SJF_delta / MET | Complete | 100 | 100 | 95 | 50 | 1,399.8 s |
| MS143 / AKT1 | Failed during conformer processing | 100 | 19 | — | — | 616.4 s |

JNJ1013 and MS143 logged `RuntimeError: Invalid conformer input` and exited 1 before final scoring. Their per-attempt linker logs and local PDBs are partial raw outputs. SJF8240 exited 0 with all 100 linker trials recorded as rejected. The original serial queue exited 1 on JNJ1013. The final two approved controls ran under the original approval and frozen settings in a separate continuation, with no repeated attempts.

The raw archive is `raw_results/`. Each attempt retains its runner log, submission and engine reports, raw stage outputs and scored structures where available.
