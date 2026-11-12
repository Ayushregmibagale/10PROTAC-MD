# SJF8240 full-frame per-lysine delta_reach_50A validation

Read INTERPRETATION.md for the plain-language result and RESULTS.md for methods/all targets. This version recomputes the exact per-lysine50-A MD-frame occupancy contrast on201 frames per retained pose, with four fixed geometry/mask combinations and protein-level statistics. It is a scoped SJF8240 validation, not a universal E2 reference generator for other compounds.

No old result, raw structure, trajectory, label or model weight is modified. No new MD or GNN is fitted. Unknown SIK3/RPS6KB2 are excluded from validation. Source/config manifests are frozen, all141 case packets are source-hash bound, published legacy values are replayed, and eight cases are independently re-read with standard MDAnalysis unwrap. See FULL_FRAME_CHECK.docx for a step-by-step record.

Reproduction: first run `python scripts/check_inputs.py` to rehash every delivered packet and source input. Then use `python scripts/extract.py --workers 2`, `python scripts/analyze.py` and `python scripts/verify_replays.py` within this scoped version. The frozen extractor validates its cached source/packet hash but does not rehash all external files on the cached fast path; the explicit input audit is required before cached reuse. It checks all141 packets and every unique bound external file. For a fresh general study, use the project-level guarded workflow rather than these frozen research entry points.

The extract and reporting manifests bind their respective fixed sources. These scripts consume the existing audited SJF8240 role/time/contact evidence and are not a substitute for per-compound calibration.
