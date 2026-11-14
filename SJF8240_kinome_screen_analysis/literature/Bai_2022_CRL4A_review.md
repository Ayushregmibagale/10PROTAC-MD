# Applying the CRL4A ubiquitination-ensemble idea to SJF8240

Review date: 14 November 2026. Literature review and proposed analysis only; no new structural scoring, MD, GNN fitting or biological validation was performed.

Bai et al., JBC (2022), DOI10.1016/j.jbc.2022.101653, combined ternary and CRL4A conformational ensembles and assessed exposed lysines near ubiquitin. They used16Å as a modeling cutoff and experimentally investigated ubiquitination. This is a geometry-based priming model, not a universal predictor of degradation. Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC9019245/

A directly relevant VHL structural reference is Crowe et al. (2024), Mechanism of degrader-targeted protein ubiquitinability, DOI10.1126/sciadv.ado6492, PDB8RWZ/8RX0. It provides experimentally determined CRL2–VHL machinery and a target-facing ubiquitination zone. Specific E2 identity, captured reaction stage and trapping strategy must be checked before using these structures as templates. Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC11468923/

## Proposed adaptation — our recommendation, not a result from this project

Keep pRosettaC for generating multiple chemically credible ternary arrangements. Add a distinct CRL2-VHL scaffold-compatibility stage. Fit through verified shared VHL/Elongin components, retain full scaffold/E2/ubiquitin geometry, and evaluate whether an exposed target lysine can approach the ubiquitin transfer site without obstruction. Define the actual atoms, for example lysine NZ and donor ubiquitin Gly76 carbonyl C, rather than transplanting a distance cutoff between different atoms. Model approach orientation as well as distance. Use alternative scaffold conformations without treating duplicate structures as independent evidence.

Report ternary contact persistence separately from ubiquitination compatibility. For each target, record the fraction of all evaluable samples that combine correct-head binding geometry, clash-free scaffold accommodation, lysine exposure and plausible approach. Also report the conditional fraction among ternary-contact frames. Invalid input/missing data remain separate QC failures; genuine target–scaffold clashes count against compatibility, not silently disappear from the denominator. These fractions are sampling-dependent structural scores, not thermodynamic or cellular probabilities. Report equal-pose and alternative cluster-weight sensitivity; do not tune weights on the eight labels.

Test the six confirmed degraders against ABL1/AXL using the existing MD where possible. SIK3/RPS6KB2 remain unknown prediction cases. Use target-level evaluation and group all poses/frames of a target in the same fold. Eight independent labelled targets remain eight examples, however many lysines or frames are scored. If improved graph inputs are later tested, represent scaffold/E2/ubiquitin spatial context and per-lysine exposure/accessibility explicitly, and use matched ablations without changing the held-out evaluation.

## Current local limitations

The existing full-frame reach50 calculation measures the fraction of frames whose nearest modeled E2 point lies within50Å, followed by a bridge-minus-nonbridge contrast. Its points represent E2C85 Cα positions, not ubiquitin Gly76. The corrected results overlap between degrader and non-degrader controls and do not validate the proposed full-machinery test. Source files: SJF8240_reanalysis/full_frame_e2_delta50_v1/scripts/extract.py and RESULTS.md.

Historical chemistry, actual-coordinate stereochemistry and static head-reference problems remain separate provenance limitations. Existing MD can support a retrospective historical-model test, but the new mechanistic interpretation cannot silently establish supplier-SJF8240 equivalence. Preserve the archived and supplier-matched branches separately. Current geometry repairs do not constitute an evaluated ubiquitination classifier.

## Concrete first experiment

Build a small, independently checked VHL scaffold ensemble and validate assembly against a known VHL ternary/catalytic reference. Project selected existing frames from all eight labelled proteins into it without fitting to their degradation outcomes. Compare productive-compatibility fractions with plain true-VHL bridging and the existing reach50 proxy, with cutoff/pose/scaffold sensitivity fixed before reading outcomes. Inspect ABL1/AXL failures and at least one degrader success visually. Proceed to broader screening or GNN changes only if the geometry and target-level comparison justify them. No universal threshold or guaranteed separation is assumed.
