# Partial repaired-reference findings — 16 November 2026

Only MZ1 is complete and assessed; GSK215 is still running. Source/input/output integrity passes for the completed task. All50global starting orientations are byte-identical to the frozen baseline. No new Ernie jobs or MD.

MZ1:168scored,129operational-QCretained;1passes the unchanged coarse structural conjunction gate, bestRMSD4.4000A versus8.3451A before repair. Passing posecombined_38_7_0001 ranks16/129by rawtotal score, recovers60%of native residue contacts, and does not establish reliable ranking or convergence. One seed only; not a biological prediction.

Independent graph/stereochemistry review succeeds for all129retained poses. The passing pose has ligand–protein maximum vdWoverlap0.3795A (native0.5349A), E3–POIoverlap0.1832A(native0.3516A), and intraligand nonlocal overlap0.6307A(native0.4433A). These continuous diagnostic values are not a complete force-field/physical acceptance certificate. Across all retained poses,39have ligand–protein overlap more than0.5A above native;4have E3–POI overlap more than0.5A above native. OperationalQC alone remains insufficient.

Next: finish GSK215, run full independent integrity/recovery/steric review, archive complete raw outputs, then assess whether this repair merits further seeds and a bounded six-compound technical pilot. Do not release production based on one recovered pose.
