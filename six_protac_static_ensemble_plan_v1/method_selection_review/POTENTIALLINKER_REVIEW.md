# PotentialLinker: preliminary source-checked review

Checked 16 November 2026. Paper: *Ubiquitination-Informed Reinforcement Learning Linker Generator Enables Discovery of Novel FGFR4 PROTAC Degraders*, DOI https://doi.org/10.1021/acs.jmedchem.6c02194. Crossref records online publication on 16 November 2026.

## Evidence actually accessed

Publisher-deposited abstract and bibliographic metadata retrieved through https://api.crossref.org/works/10.1021/acs.jmedchem.6c02194. Full publisher text was inaccessible during this check. This is an abstract-level review, not a full methods audit.

The authors describe PotentialLinker as reinforcement learning for linker optimization. A normalized weighted harmonic mean (NWHM) score predicts degradation activity. They report ternary-structure benchmarking without precomputed protein–protein complexes, plus STAT3/RIPK1 case studies recovering structure–activity relationships and suggesting E3 ligases. Guided design yielded FGFR4 degrader D11: DC50 56.23 nM and Dmax 89.58%. The abstract attributes selectivity to covalent interaction with FGFR4 Cys552 and reports oral bioavailability 21.3% and tumor growth inhibition 87.07%.

## Relevance and limits for our campaign

The experimental discovery makes this relevant evidence for mechanism-informed linker design. Our immediate question differs: characterize six fixed purchased molecules across kinase targets, with static ensembles and no new MD. We do not need to generate replacement linkers to answer that question.

A potentially reusable component is the ubiquitination-informed scoring function, as an external baseline or a source of structural features for later GNN comparison. That reuse remains conditional on inspecting the actual equation, structural construction, code/license and validation design. The abstract alone does not establish full CRL modeling, reference structures, distance atoms/cutoffs, clash handling, NWHM components, GNN use, or released implementation. A search found a repository listing named PotentialLinker on Gitee, but authorship, contents and paper linkage were not verified; do not describe it as verified author code.

No claim is made that their score distinguishes our confirmed binder/non-degraders from degraders, or that our pipeline reproduces theirs. No result from this paper validates the historical SJF8240 connectivity branch or supplies experimental labels for SIK3/RPS6KB2.

## Decision

Keep Ernie preparation and the running docking-repair qualification unchanged. Before adopting any scoring component, obtain the full article/supporting information, map the NWHM terms and ligase geometry, inspect author-linked implementation and assess applicability to VHL and fixed-molecule target profiling. Evaluate any later score using independent target-level labels, without selecting thresholds on held-out targets. Linker generation is outside the current six-drug campaign scope.
