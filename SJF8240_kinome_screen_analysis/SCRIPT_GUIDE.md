# Script guide

Use `/home/ayushregmibagale/miniconda3/bin/python3` for these scripts. Source paths are recorded in `config/protocol.json` and `config/pools.json`; they are specific to this workstation. Dependencies include NumPy, pandas, SciPy, Biopython and MDAnalysis. The MD stage additionally imports the existing reviewed `protac_qc` and read-only trajectory helpers from the original analysis repository.

| Order | Script | Purpose | Main output |
|---|---|---|---|
| 1 | `01_inventory_kinome_inputs.py` | Inventory available clustered Results poses; hash-check duplicate copies; freeze identity, geometry and QC settings. | `results/inventory.tsv`, `config/pools.json`, `config/protocol.json` |
| 2 | `02_screen_static_bridging_poses.py --workers 8` | Verify protein identities, measure true-VHL/kinase contacts and pocket proximity, and flag severe overlaps. | Per-pool tables and completion receipts under `results/pools/` |
| 3 | `03_collect_top10_static_poses.py` | Require completion of all available pools, verify result hashes, rank passing poses and rehash selected source PDBs. | `results/top10_static_prospective_poses.tsv`, `results/top10_static_distinct_prospective_kinases.tsv`, all-pose ranking |
| 4 | `04_rank_existing_md_bridging.py` | Rank already completed MD by simultaneous true-VHL/kinase heavy-atom contacts across 201 saved frames. Reuse audited core contacts and recompute SIK3/RPS6KB2 contacts. | `results/top10_MD_prospective_poses.tsv`, MD target summary and validation receipts |
| 5 | `05_check_package.py` | Check packaged source hashes and syntax without rerunning docking, MD or analysis. | Terminal PASS or explicit failure |

Stage 4 is independent of stages 1–3. Do not launch duplicate jobs while the current static screen is active. Stage 1 writes configuration: do not rerun it during stage 2 or use it to overwrite an established run without archiving that run. Stage 3 deliberately refuses incomplete screens. A rejected pool remains rejected, not a non-bridging result. Cached static results verify code, configuration, result-table and reference-context hashes; a fresh complete rescan is needed if raw input poses change. Selected top poses are rehashed during collection.

## Geometry and interpretation

Legacy kinome selection called chain A VHL, although it contains ELOB plus kinase segments. The new analysis resolves protein identity and TER-separated segments. Many saved docking models contain only a VHL fragment. Full-VHL coordinates are therefore reconstructed by fitting the retained E3 reference through ELOC; this is a model-based check, explicitly distinct from directly observed coordinates. Both distances are retained.

Static eligibility uses simultaneous contacts below 5 Å, proximity to the reference-head pocket regions, at least three ligand atoms contacting each pocket region, reference-fit checks and an operational 1.5 Å severe-overlap screen. Eligible poses are ordered by the smaller of the two pocket contact-atom counts, then the worse pocket distance, then source path. Rosetta energies are annotations, not cross-protein degradation scores. Pocket proximity does not establish atom-mapped engagement of both cognate ligand heads.

The top-ten individual poses may contain repeated kinases; the distinct-kinase table selects one pose per accession. Prospective tables exclude the eight current labelled targets. SIK3 and RPS6KB2 have unknown experimental degradation status. Their exploratory labels are not validation evidence. Existing prospective MD is available for these two kinases only; static poses are not MD results.

## Validation status

Packaging checks establish source integrity and Python syntax only. Static screening is still running at initial packaging. MD ranking completed for 174 poses, including 33 prospective poses across SIK3 and RPS6KB2. Original datasets are unchanged. These scripts are not a prospectively validated degradation classifier. See per-pool receipts, rejection reasons and the existing validation directory for numerical evidence and limitations.

## Legacy scripts

The old `protac-md/scripts/01_setup/09_run_kinome_md_screen.py` is retained as historical provenance; it is not called here. Its chain-A interpretation, random pose sampling and accession mapping make it unsuitable as the current selection method. The numbered MD/GNN scripts in `SJF8240_renanalysis2` remain a separate downstream workflow; this folder supplies the missing broad-screen selection stage.


## Stronger head and clash review (completed)

Run stages 6–8 after static collection:

- `06_check_head_engagement.py --workers 6`: checks all 723 original QC survivors, maps the two chemical heads, compares contacts at 4/4.5/5 Å, audits element-specific overlaps and reference-head integrity. It does not perform structure repair or new MD.
- `07_verify_engagement_results.py`: independently replays previous top-ten raw-PDB kinase distances and shared linker overlaps, and checks result hashes.
- `08_write_engagement_report.py`: writes Markdown/Word findings, compact decisions and PNG/PDF comparison figures. Requires python-docx and matplotlib in addition to RDKit and the analysis dependencies.

Read `engagement_v2/FINDINGS.md` or `engagement_v2/FINDINGS.docx`. **None of the original top ten passes the stronger clash check.** Both head reference geometries are untrusted under their saved topology and all 723 retained static poses have a shared internal linker overlap. The original contact-count top ten is superseded for selection. Empty qualified top-ten tables are intentional; do not use them as evidence of absent kinase binding. New repair/regeneration work is required before a defensible MD shortlist can be approved.

The engagement review scripts, results, receipts and figures are physical files in `engagement_v2/`; the earlier bulk results/logs/validation links remain as documented. The stage-6 backend accepts `--analysis-root` and `--output` for a fresh archived comparison. The package check validates source integrity, not biological performance.


## Geometry repairs and chemical identity (14 November 2026)

The repair stage is in `repair_v3/`. Read `repair_v3/FINDINGS.md` or the Word document before selecting a pose. Three repaired historical-model poses—FER, TGFBR1 and PRKG1—pass the declared geometry/contact diagnostic. They are not supplier-matched SJF8240 binding predictions. The archived molecule has different connectivity and actual coordinate stereochemistry; old results have not been relabelled.

| Order | Script | Action |
|---|---|---|
| 9 | `09_audit_molecular_identity.py` | Record historical/supplier chemical graphs and save native-parent reference controls. |
| 10 | `10_rebuild_reference_heads.py --workers 4` | Regenerate VHL heads and test corrected kinase-head transfer for all198 shortlisted pools; rejects existing output rather than overwriting. |
| 11 | `11_repair_historical_poses.py --attempts 2 --workers 4` | First bounded local minimization of the original ten; preserves actual archived molecular graph/stereo and fixed proteins. |
| 12 | `12_refine_historical_poses.py --attempts 2 --workers 4` | Stronger overlap-penalty stage with unchanged external screening gates. |
| 13 | `13_verify_and_report_repairs.py` | Replay geometry, verify identity/fixed proteins and all198 VHL heads; regenerate full-E3 diagnostic PDBs, tables, Word/Markdown report and figures. |

Stages9–13 are already completed locally. Archive output directories before rerunning10–12: these stages deliberately refuse existing receipts. Use `KINOME_ANALYSIS_ROOT` only when referencing a different archived kinome analysis package. They are workstation-specific and require the reviewed reference sources in the original workflow. Neither minimize stage performs docking or MD.

All198 VHL reference heads pass geometry checks. All198 archived kinase-head transfers fail the retained pocket/domain criteria. The native MET self-alignment positive control passes. Therefore the broad supplier-SJF8240 top-ten selection remains unresolved and a suitable target-reference/candidate-generation stage is needed. Numerical success is not biological validation. Full structures remain local; compact reports, receipts and scripts may be versioned.
