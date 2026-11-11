# SJF8240: matched holdouts and lysine-aware MD GNN

Three user-authorized stages: original GNN leave-pair-out, per-lysine mechanistic analysis, and a lysine-aware four-window GNN under identical matched holdouts. Six confirmed/user-labelled positives and two negatives; SIK3/RPS6KB2 remain experimentally unknown. No new MD or docking.

Read docs/STEP_BY_STEP.docx (Word) or docs/STEP_BY_STEP.md. RESULTS.md will record the final comparison. Compact tables are under results/. Full local weights, receipts and CPU replays are under analysis/{baseline,temporal}/{elob,vhl}/; lysine packets are under analysis/lysines/poses/.

config/study_plan.json is the new study specification. config/protocol.json and user_labels.json are unchanged v1 parameter/label contracts required by the reused loader. Active source freezes are pair_source_hashes.json and temporal_source_hashes.json. The original v1 execution manifest is archived under provenance/ and is not the current executable manifest.

Training: scripts/pair_train.py and scripts/temporal_train.py, each --role elob|vhl. Independent CPU replay: add --device cpu --verify, with the Python 3.10 GPU PYTHONPATH overlay removed. Extraction: scripts/extract_lysines.py --workers 2; graph aggregation: scripts/build_temporal.py. Physical scales are fixed; no scaler is fitted to held-out targets. All four windows for one target are processed jointly and stay together in every split.

The temporal model uses 29 node inputs, a shared two-layer GAT over four chronological window graphs, then embedding mean, standard deviation and last-minus-first pooling. It changes feature width, time representation, sampled edge aggregation and head width simultaneously. Improvement would support this combined candidate, not identify the causal contribution of any one feature. Both models use 300 epochs, seed 42 with identical pair-bound seed derivation, Adam 0.001, weight decay 0.0001, cosine schedule, unweighted BCE and threshold 0.5. No tuning follows scores.

Matched-pair AUROC is the fraction of 12 within-model positive/negative comparisons ranked correctly, with ties half credited. It is not pooled cross-fold AUROC. Only two independent negative targets exist. Scores/accuracy are descriptive, uncalibrated, and not external validation. Exposure/E2 point geometry does not establish ubiquitination or degradation.

Feature units: SASA and burial are Å². The summary column `mean_min_e2_A` averages each lysine's nearest E2-point distance over lysines, frames and poses; it is not the earlier graph scalar that took the minimum over all lysines first. Exposure is NZ SASA >5 Å². The joint occupancy uses each lysine's nearest point <80 Å plus simultaneous true-VHL/ligand/target heavy-atom bridging (<5 Å on both sides). Point reach fraction is the fraction of modeled E2 points within 80 Å of that individual lysine. These are not catalytic probabilities.
