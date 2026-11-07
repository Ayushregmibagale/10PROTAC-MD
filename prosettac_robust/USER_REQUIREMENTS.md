# Persistent user requirements — PROTAC-MD / SJF8240

Updated 2026-11-06. Read this before planning or changing this project. Later user messages take precedence. This records project instructions; it does not claim an account-wide memory feature.

## Current priority and scope

The user has already completed the SJF8240 experiment through PRosettaC, MD and GNN. Preserve and reuse the same existing structures and completed results. The primary task is to make the PRosettaC, MD, analysis and GNN scripts robust, reproducible and scientifically defensible, and audit the rationales (explicitly including bridging/contact criteria such as <5 Å). Do not treat robustness work as an instruction to repeat the expensive experiment.

The user subsequently clarified that helpful work already underway may continue, and explicitly authorized continuing useful bounded diagnostics, including the new minimization check. Retain helpful fixes, audits and small integration tests. This does not authorize a wholesale PRosettaC/MD/GNN rerun or Ernie production submission. Keep experimental method changes separate and do not retrospectively describe old results as if a new method produced them.

## Preservation and locations

- Original `/home/ayushregmibagale/protac-md` is a read-only reference; build on its successful workflow without modifying it.
- Original `D:\PROTAC degradation Prediction` (`/mnt/d/PROTAC degradation Prediction`) is a read-only source of the completed SJF8240 pipeline/results.
- Preserve prior top-level work in `/home/ayushregmibagale/10PROTAC-MD`.
- Put derivative work in `/home/ayushregmibagale/10PROTAC-MD/prosettac_robust/`, with scripts under `workflow/`. The user explicitly requested a separate folder and no changes to previous work.
- Stage changes separately, review/check them and preserve backups and source/input hashes. Do not overwrite original structures, trajectories, labels, predictions or model files.

## Scientific requirements

- Aim for rigorous, publication-quality work suitable for a high-impact submission; do not promise Nature acceptance or infer scientific validity from passing software tests.
- Focus on SJF8240 first because its end-to-end results already exist; extend justified fixes to other compounds afterward.
- The compounds were purchased from Tocris/Bio-Techne. Use supplier structures as chemical reference evidence, while honoring the user's instruction to retain the structures used in the existing experiment. Investigate discrepancies explicitly; do not silently substitute molecules, relabel trajectories, discard results or demand a full rerun as a default response.
- Describe precisely what was actually simulated. Existing chemistry/provenance concerns remain audit findings to resolve, not permission to claim unverified identities or to rerun the campaign automatically.
- Audit bridging/contact thresholds against actual code: atom selections, < versus <=, units, periodic boundaries, simultaneous versus marginal contacts, aggregation, sampling interval and missing-data behavior. A geometric contact proxy does not by itself establish degradation, thermodynamic stability, binding free energy or causal mechanism.
- Use existing trajectories for justified analyses or sensitivity checks wherever possible. Preserve original analysis outputs; write new audits separately.
- Keep labels/assay decisions independent of model predictions. Missing data and failed finite searches are not biological negatives. Historical labels must not be changed to improve accuracy.
- Preserve the existing GNN architecture/results and assess its role honestly. The user asked how essential GNN is; this is not permission to discard it or claim new predictive performance.
- Separate software fixes, bounded diagnostics, scientific-method changes and completed study results. Record exactly which source/structure/trajectory generated each artifact.

## Compute and deployment

- The user accepted a 7.24 TB storage allowance. This is a storage preference, not an allocation/quota guarantee or blanket authorization for new production calculations.
- Audit readiness before Ernie production. Do not guess the writable project/scratch allocation, submit a study campaign or deploy replacements into original directories.
- Helpful bounded local controls are authorized. State their purpose, budgets and limits. They must not become an unnoticed repeat of the original experiment.
- The new minimization tests local script/parameter handoff on a separate diagnostic pose. It does not replace the existing SJF8240 MD or prove degradation/stability. Preserve the interrupted attempt and any retry separately.

## Documentation, Git and communication

- Keep detailed notes of every material action: what changed, how it was done, why it was done, validation, failures and limitations.
- Keep committing and pushing appropriate source/audit checkpoints to the user's GitHub repository. Keep licensed binaries, private runtimes and bulky raw outputs out of Git.
- Continue authorized work without repeatedly asking whether to proceed. Do not stop at a plan or arbitrary checkpoint. Explicit pauses must be honored; later continue messages resume within the latest scope.
- Give concise, factual progress updates and tell the user when fixes are actually finished. Do not imply unattended work continues after a session ends or call unfinished scientific validation complete.
- Earlier 'continue until finished, even months' messages express persistence; they do not override the later instruction to prioritize reuse of the completed experiment.

## Outstanding decisions — do not invent answers

Assay-specific labels/endpoints and the correspondence of the historical 17-target set to primary experimental evidence still need review. The user has not selected a replacement assay or authorized relabeling. An Ernie storage/allocation location has not been confirmed. Keep independent read-only/code work moving while these remain unresolved.

## Findings to carry into future work

The retrospective compiled-topology audit now identifies all 141 core existing MD systems: the historical group labelled VHL is Elongin B/Q15370 (chain A); true VHL/P40337 is chain C, with 100% sequence identity/coverage. This is an analysis-role finding, not evidence that the proteins are absent or that MD must be repeated. Preserve the old metric/model outputs; evaluate and label any corrected analysis separately using the same trajectories. First-case distance reproduction and the ongoing 141-case contact re-analysis are documented under `existing_result_audit/`. Do not confuse successful numerical reproduction of the old metric with validation of its biological label.

The 141-case contact re-analysis is now complete: all stored-distance reproduction checks pass within 0.000022 Å. Corrected true-VHL and heavy-only contact definitions are separate derivatives under `existing_result_audit/contact_core_v1`; original metrics/models remain preserved. This was analysis of existing trajectories, not new MD. One static representative per target (17/17) also selects Elongin B under the old VHL label; do not overstate that as an exhaustive static-pose audit. The corrected E2 script uses the independently consistent Elongin-C frame anchor, so avoid blanket claims that all E2 geometry is wrong.

The user also requested continuing overnight. Continue useful work while the session is active, save frequent durable/Git checkpoints, and do not promise unattended work after the session ends. The requested overnight work does not alter the no-wholesale-rerun scope.

Latest clarification: focus on the GNN; the user does not want separate conventional ML models. Continue PRosettaC/MD integration hardening and preserve all completed results. Scalar-feature checks are part of GNN input validation, not a request to add another predictive model.
