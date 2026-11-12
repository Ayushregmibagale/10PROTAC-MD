# Reviewed reuse path for the PROTAC analyses

This package accepts compound-specific manifests. It never substitutes SJF8240 labels, paths, ligand atom counts or chain numbers for another compound. Existing completed experiments remain frozen. Other compounds can produce different results.

Use the project-level `RUN_ANALYSIS.py` entry point:

```
python RUN_ANALYSIS.py validate --study /path/to/study.json
python RUN_ANALYSIS.py batch --batch /path/to/batch.json --output /path/to/preflight.json
python RUN_ANALYSIS.py extract --manifest /path/to/trajectories.json --output /fresh/output
python RUN_ANALYSIS.py train --study /path/to/study.json --output /fresh/run/vhl --role vhl --device cpu
python RUN_ANALYSIS.py replay --study /path/to/study.json --output /same/run/vhl --role vhl
```

For GPU fitting, invoke the launcher with a CUDA-compatible Python and PyG runtime and request `--device cuda`; unavailability is an error, not an automatic CPU fallback. For independent CPU replay, use the base CPU Python with the GPU-only PYTHONPATH overlay removed. Logs and outputs belong to the caller's explicitly selected directory. Linux/WSL file locks prevent simultaneous runs using the same output. Training folds resume only with identical configuration, graph hashes, source and training runtime. Incomplete unreceipted weights are archived before refitting that fold. Extraction deliberately requires a fresh output and preserves partial failures; it does not pretend to resume an incomplete extraction.

`make_sjf8240_examples.py --project-root /path/to/SJF8240_reanalysis --output /path/to/examples` migrates only the existing SJF8240 inputs. It produces separate original-MD and lysine-temporal example manifests with 300 epochs/seed42. It does not fit them. Unknown SIK3/RPS6KB2 are prediction-only in the original example; the enriched example excludes them because audited enriched inputs are absent. This is a migration adapter, not a default cohort for other compounds.

## Required per-compound inputs

A study manifest supplies compound identity; explicit angstrom/ns units; model (`original_md` or `lysine_temporal`); evaluation (`leave_pair_out` or `lopo`); frozen optimizer/epochs/seed/loss/threshold; VHL/ELOB role accessions and hash-bound identity evidence; unique protein identity groups; target labels, label status and hash-bound JSON label evidence; and every role graph with its path/hash/compound/target/role. Confirmed-user labels are distinguished from confirmed experimental evidence. Hypothetical labels are rejected. Unknown labels must be null and are excluded from fitting and metrics. Missing graphs, nonfinite inputs, changed scopes, changed labels or duplicate target aliases stop execution. No absent scalar/geometry is filled with zero.

The trajectory manifest supplies hash-bound TPR/XTC files, a hash-bound reference-sequence table (`sequences` keyed by UniProt accession), the target-accession map, exact model sequences and CA selections for VHL, ELOB and target in every pose, ligand selection/heavy-atom count, exact expected time grid, duplicate-frame policy, contact threshold and edge threshold. Roles are checked against canonical sequences as well as the declared model sequence, so selecting ELOB while calling it VHL fails. Each sampled frame uses periodic distances and an actual simultaneous contact outer product. Poses receive equal weight and residue correspondence must match. Empty contact graphs remain empty; no edge is fabricated. Topology/trajectory identity and underlying MD/chemistry quality still require their own provenance review.

Supported geometry: VHL/ELOB systems only. CRBN or another ligase needs a separate identity/model specification; do not rename it VHL. Other compounds need their own binding/degradation labels and assay context. A ten-compound batch lists all expected compound IDs; omission of any listed compound is an error. Each study currently fits within one compound. Joint cross-compound training is not implemented; future evaluation must group repeated biological targets across compounds.

## Model and feature boundaries

The original MD architecture is preserved from corrected_md_only_v1. The enriched architecture is preserved from mechanistic_v3 (29 ordered node features, four chronological windows, mean/std/trend pooling). These are separate candidates. The configurable runner handles both graph formats and saves all fold weights/scores, membership, source/input hashes and CPU replay receipts. One graph or one four-window packet represents one target; poses/frames are not independent labelled targets. Its JSON-encoded held-out-set seed derivation is explicit and differs from the older name-delimited study seed; new runs are new versions, not replacements for historical numbers.

The generic extractor produces original-MD contact graphs. The existing enriched SASA/E2 extractor remains a frozen SJF8240-specific research pipeline. For another compound, enriched graph generation needs a separately reviewed configuration for E2 reference-frame calibration, whole-component imaging, lysine mapping, time windows and feature units/scales. The runner accepts validated enriched graphs; it does not claim to generate or validate a new E2 prior. Do not copy the SJF8240 point cloud, residue correspondence or chemistry assumptions blindly. Full E2–ubiquitin catalysis is not modeled.

Matched-pair AUROC compares the held-out positive/negative inside one fitted model. Pooled LOPO AUROC is explicitly named and flagged as a different estimator. Sensitivity, specificity and balanced accuracy use frozen threshold0.5. Pair counts do not increase the independent negative-protein count. No post-score threshold or seed optimization is performed.

## Validation scope

The generic extractor was replayed on all141 reviewed SJF8240 cases /28,341 frames and reproduced16 historical role graphs. Both GNN architectures receive short CPU plumbing tests with save/resume/source-drift checks, target grouping and separate CPU replay. Those one-epoch synthetic tests are not new degradation results. Other nine compounds have not been run or certified by this audit; they must pass their own complete preflight, extraction and model replay.
