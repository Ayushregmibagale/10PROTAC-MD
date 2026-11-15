# Prepared coverage and geometry review — 15 November 2026

The broad manifest contains 690 identifiers. Current domain definitions cover 435 proteins in 770 structure/domain rows; all 770 files exist. This does not establish that all compounds bind those proteins.

| Preparation route | Candidate structures | Unique compound–protein pairs | Hash/identity checks |
|---|---:|---:|---:|
| Original CE transfer | 112 | 102 | 112/112 pass |
| Explicit-pocket comparison | 161 | 130 | 161/161 pass |

These are alternative preparation cohorts with overlapping protein/structure variants. Do not add them as independent biological samples. MS143 is not present in either cross-kinase cohort.

Most transferred kinase heads add steric overlap relative to their corresponding bound-head source. Overlap increases exceed 0.5 Å in 92/112 original candidates and 131/161 explicit-pocket candidates. This is an exploratory geometric diagnostic, not an affinity measurement or a degradation classifier. Protein repacking or a different receptor/ligand conformation could resolve some cases; suitability must be checked rather than assumed.

A universal 0.6 Å overlap rejection would also flag the native GSK215 VHL head (baseline 0.647 Å), so that absolute flag cannot independently qualify inputs. The report therefore retains continuous source baselines and overlap changes. None of the VHL transfers adds more than 0.2 Å to its corresponding baseline; the main transfer problem is on the kinase side. The modeled SJFδ VHL reference itself has 0.878 Å maximum overlap and needs review.

No candidate has been promoted to production-ready. Original inputs and rejection statuses are unchanged. Inspect PREPARED_HEAD_REVIEW.tsv and COVERAGE_BY_COMPOUND.tsv; source hashes and limitations are in COVERAGE_AUDIT.json. Initial uncalibrated diagnostics remain under ../coverage_audit_initial/ for provenance.

Next: review representative added overlaps, distinguish movable side-chain clashes from backbone/pocket incompatibility, assess receptor-state alternatives and pocket-local ligand refinement, and validate any repair against known controls before broad release. Do not loosen a threshold merely to reach 600 kinases. No MD is required or authorized.
