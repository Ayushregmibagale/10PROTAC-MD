# Remediation status

Publication repairs are installed in reviewed batches. **No Ernie production
campaign has been launched. Scientific release evidence is still incomplete.**
The original protac-md workflow is the baseline and its source remains unchanged.
The exploratory replacement engine is retained separately and excluded from the
production controller. All failed diagnostics and legacy inputs are preserved.

| Item | Current evidence |
|---|---|
| Supplier identity | All five Tocris stereochemical graphs verified; distinct heads and anchors |
| Legacy inputs | All 2,860 preserved and excluded from campaign release |
| KLIFS provenance | 496 exact coordinate fingerprints verified; RSK3 fallback excluded |
| Domain inventory | 770 domains across 435 UniProt proteins, explicit exclusions |
| Kinome transfer | Original CE: 112 candidates. Separate explicit-pocket comparison: 161 candidates, all five known-target controls represented; neither is a released complete-kinome study |
| Known-target preparations | PTK2/GSK215, MET/SJF8240, EGFR/SJF1521, IRAK1/JNJ1013 and a separate modeled MAPK13/SJFδ control |
| SJFδ MAPK13 control | Published 5IA4-to-4EYJ template pair; 83/85 mapped pocket residues, 1.608 Å pocket RMSD, 1.794 Å minimum heavy contact; no experimental p38δ/foretinib claim |
| Baseline integration | Original stage order with repaired both-head-constrained linker generation; valid native-frame GSK215 output and receipt-compatible no-valid-pose outcome |
| Chemistry components | All five full-molecule Rosetta parameter/load/stereochemistry controls passed |
| MD parameters | Five neutral GAFF2/AM1-BCC sets and GSK215 +1 set verified, including QM-output stereochemistry and atom-map hashes |
| Protonation | Provisional pH 7.4 enumeration plus 6.5–8 sensitivity; GSK215 neutral/+1 at 7.4, others neutral in this tool; no populations or final MD state selection claimed |
| MD completion checks | Actual final step, EM force convergence, exact TPR generation binding, parent output and prepared topology/force-field hashes required |
| MD systems | Neutral and +1 GSK215 strict topology builds passed; modeled protein loops and unchanged observed heavy atoms recorded |
| MD dynamics | 1,000-step EM and 0.04 ps NVT mechanical diagnostic only; EM did not converge and no scientific equilibration is claimed |
| MD analysis | Atom-mapped periodic contacts and trajectory/topology provenance verified on diagnostic trajectory |
| Execution | Bounded pilot controller, immutable attempts, raw scores, QC rejections, process-group interruption and verified resume/skip contracts |
| Storage | Full raw retention accepted; measured 112-pair model ~7.24 TB before MD/retries/backups; actual writable allocation still to be identified |
| ML | Receipt-derived geometry schema records rejections; old models are not assumed compatible |

Latest software validation: 77 regression tests passed in the staging checkout;
installation-specific output is in audits/installed_regression_<timestamp>.txt.
A test pass is a software result, not publication certification.

Before committing months of compute, remaining evidence includes representative
independent global recovery/cross-docking and sampling pilots; treatment of the
restricted transfer coverage; selection of relevant chemical states; measured
runtime/memory/storage on the intended Ernie allocation; and converged,
replicated MD with a prespecified analysis. No versioned campaign release has
been fabricated to bypass these requirements. See REMEDIATION_LOG.md for details.

See NATURE_READINESS.md for the requested Nature reporting standard, primary
guidance links and the evidence still needed to support scientific claims.
The 161-candidate comparison projects ~10.41 TB at the same measured budgets; it
has not replaced the original 112-pair plan or expanded the accepted 7.24-TB scope.
