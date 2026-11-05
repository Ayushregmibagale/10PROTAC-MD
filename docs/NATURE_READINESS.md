# Nature publication-readiness evidence

Target requested by the user: Nature-level scientific rigor and reporting.
This document is an evidence inventory, not an acceptance prediction or completed
reporting summary. An executable workflow alone does not establish novelty,
mechanism, broad interest or biological validity.

Nature's [initial-submission guidance](https://www.nature.com/nature/for-authors/initial-submission)
requires a reporting summary for relevant life-science research and an explicit
code-availability statement when custom code is central to the conclusions.
The [Nature Portfolio reporting policy](https://www.nature.com/nchem/editorial-policies/reporting-standards)
covers availability of supporting data, materials, code and protocols. GitHub
checkpoints provide inspectable source; the finalized supporting dataset still
needs an appropriate repository/accession and a frozen release.

The Communications Biology perspective
[Reliability and reproducibility checklist for molecular dynamics simulations](https://doi.org/10.1038/s42003-023-04653-0)
provides relevant scientific guidance: independent replicate simulations,
statistical treatment, convergence and enough methodological detail to assess
reliability. Its checklist asks about at least three simulations per condition.
This is guidance from that article, not a claim that every Nature journal imposes
a universal trajectory duration or that three replicas guarantee convergence.

| Claim/evidence area | Current result | Needed before making the claim |
|---|---|---|
| Correct purchased chemistry | Five Tocris identities, complete stereo graphs and separate heads verified | Preserve lot/assay conditions with experimental records |
| Correct modeled inputs | Exact bound fragments or explicitly labeled analog/template models; explicit kinase domains | Review model uncertainty, chemical states, missing loops and activation states |
| Broad kinase coverage | CE comparison: 112 candidates; explicit-pocket comparison: 161; neither is a complete kinome screen | Freeze one defensible cohort/method before inspecting production rankings; report exclusions and uncertainty |
| Real docking output | Valid native-frame control; raw-score and head-geometry QC | Independent global recovery/cross-docking and representative multi-seed sampling controls |
| Full PROTAC feasibility | Fixed both-head linker generation, exact stereo and strict clash checks | Sampling convergence and head/linker strain review; no assumed thermodynamic populations from generated-pose fractions |
| Force-field reproducibility | Five neutral and one GSK215 +1 GAFF2/AM1-BCC sets verified; mapped protein/ligand topology | Select assay-relevant states; check sensitivity of key conclusions to state/parameter/model choices |
| MD physical validity | Strict portable topology and bounded mechanical diagnostics | Converged minimization/equilibration, independent replicas, adequate time to assess relevant motions and uncertainties |
| Statistics | Receipt-derived distances and rejection denominators | Prespecified independent units, replica-level uncertainty, autocorrelation/block treatment and exclusions; frames are not independent biological replicates |
| Selectivity/degradation mechanism | Structures and published identities are available | Relevant experimental positives/negatives, orthogonal validation and evidence separating binding, ternary formation and degradation |
| Reproducibility | Versioned source, immutable commands/seeds/hashes, retained failures | Frozen code/data release, availability statements and environment/runtime details for the actual production run |
| Deployment | Full raw retention preferred; 112-pair estimate ~7.24 TB | Allocated writable storage, quota/reserve, measured Ernie pilot resources, backup/deposition plan |

Project-specific analysis rules: retain failed/no-valid-pose outcomes; never call
an excluded structure or incomplete search an inactive kinase; never report raw
Rosetta scores as binding free energies; do not count correlated MD frames or
related generated poses as independent replicates; separate native-frame input
checks from independent predictive recovery. These are scientific safeguards for
this project, not assertions of literal journal-specific mandatory cutoffs.

The primary scientific hypothesis and comparison set must be explicit. Related
compounds with different VHL attachment geometry can motivate a testable
orientation/selectivity hypothesis, but simulations alone must not be described
as proof of cellular degradation. Existing experimental evidence and any new
measurements must be linked to the exact compound/target/conditions used.

No experimental outcomes, convergence results, release DOI, data accession or
publication-readiness certificate have been invented. Open entries remain open
until their underlying evidence exists.
