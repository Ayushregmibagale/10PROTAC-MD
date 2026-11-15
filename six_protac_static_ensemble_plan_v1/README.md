# Six-compound static campaign execution checkpoint

Read PLAN.md, ACTIVE_WORK.md, STATUS.json and RESEARCH_CONTRACT.md. No MD will be run for these six compounds. Remote status files here are timestamped snapshots, not a live monitor.

## Current results

- Old seed42small-budget reference assessment: MZ1 passes0/62retained poses; GSK215passes1/3but top1fails. Metrics/gate/measurement controls are in results/REFERENCE_RECOVERY.json. This is not robust structural qualification or biological validation.
- MS143 verified Tocris8081graph. strict_review_v2contains a modified35-atom head with29atoms fixed to deposited4GV1/0XZscaffold coordinates. One of12trials passes the operational0.6ÅvdWoverlap screen. Thr308phosphorylation is retained. Binding-mode/state review and docking integration remain pending. The initial root-level modeled head is rejected; see MS143/CURRENT_STATUS.json.
- Bounded four-core pilot started on Ernie. Native-frame control completed without feasible linker poses; global GSK215is running and MZ1is queued. No full-kinome production has started. Read REMOTE_STATUS_SNAPSHOT.json and query remote STATUS.json for fresh state.

## Reproduction

From this directory, with the project analysis dependencies installed:

```bash
python3 test_pilot_controller.py
python3 test_ms143_geometry.py
python3 scripts/01_assess_reference_recovery.py --runs /path/to/downloaded/runs --references reference_evaluation --output /tmp/reference_recovery
python3 scripts/02_prepare_ms143_heads.py --project .. --sources references --output /tmp/ms143_fresh_review
python3 prepare_bundle.py --project .. --inputs reference_inputs --output /tmp/new_ernie_reference_bundle
```

The private remote engine depends on licensed PyRosetta/Rosetta helper components already present in the project/runtime; these are not included in GitHub. BUNDLE_MANIFEST.json records the exact private transfer. Newly rebuilt bundles must be separately frozen/reviewed, not mixed into an active attempt. Source/evaluation inputs are separate; native ternary coordinates are not required in generation inputs.

Poll: `ssh ernie 'cat /home/s1638308/ernie_static_reference_pilot_v1/STATUS.json'`

Stop: `ssh ernie 'touch /home/s1638308/ernie_static_reference_pilot_v1/STOP'`

Primary structure sources: https://files.rcsb.org/download/4GV1.cif and https://files.rcsb.org/ligands/download/0XZ.cif . Native reference identities are5T35/MZ1and7PI4/GSK215. Input and source hashes are retained in the manifests and result receipts. No source, pose count or geometry score is treated as an experimental degradation label.

## Verification at this checkpoint

Four controller checks and three MS143 geometry checks passed. Reference recovery measurement controls pass for native self-comparison, common rigid transformation and deliberately incompatible placement. These checks validate execution and measurement contracts; they do not certify kinome sampling sufficiency or biological prediction quality. All final pose sets still require the independent atom-type-aware steric review described in RESEARCH_CONTRACT.md.

## Preparation coverage audit

See [coverage_audit/FINDINGS.md](coverage_audit/FINDINGS.md) for current coverage, hash integrity and source-calibrated steric review. Ten relevant checks passed in the isolated checkout. The full campaign remains unreleased; ask for submission approval only when ready. No further Ernie submissions beyond the previously authorized pilot queue. No MD.
