# SIK3 and RPS6KB2 prediction-only assessment

| model                       | role   | unknown_target   |   mean_score |   min_score |   max_score |   degrader_vote_fraction |
|:----------------------------|:-------|:-----------------|-------------:|------------:|------------:|-------------------------:|
| legacy_corrected_md_only_v1 | elob   | RPS6KB2_full     |     0.562268 |    0.52303  |    0.659893 |                        1 |
| legacy_corrected_md_only_v1 | elob   | SIK3_full        |     0.561903 |    0.523786 |    0.659404 |                        1 |
| legacy_corrected_md_only_v1 | vhl    | RPS6KB2_full     |     0.579751 |    0.541105 |    0.634672 |                        1 |
| legacy_corrected_md_only_v1 | vhl    | SIK3_full        |     0.579771 |    0.540586 |    0.634567 |                        1 |

These are predictions from the older corrected MD-only v1 model, independently replayed on CPU. Neither unknown was labelled or used in fitting. This is not experimental degradation evidence. ELOB is an alternative historical arm, not true VHL.

The AUROC 0.8333 temporal catalytic model cannot yet score these targets: their matching four-window 29-feature graphs and five catalytic-geometry annotations have not been built. Existing aggregated 3-feature graphs are not compatible. No missing catalytic features were silently zero-filled. Newest-model unknown classification remains unavailable until matching inputs are prepared and verified from existing trajectories.
