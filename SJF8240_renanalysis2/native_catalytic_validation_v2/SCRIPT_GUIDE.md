# Ordered scripts and reproducibility

Use /home/ayushregmibagale/miniconda3/bin/python3. Dependencies:NumPy,pandas,SciPy,MDAnalysis,Biopython,scikit-learn,python-docx,matplotlib;the GNN runtime smoke test additionally uses PyTorch/PyG with the existing native-backend guard. Do not run with Python optimization flags. All scripts are already run and results verified in this version.

| Step | Purpose |
|---|---|
|00_verify_reference_inputs.py|Verify four public PDB/CCD inputs and their source-bound hashes.|
|01_prepare_and_validate_native.py|Prepare deposited/proxy references;correct native-target backbone/sidechain mapping;record every transferred atom and all control sites.|
|02_extract_md_native_references.py --workers 4|Score existing141 MD cases,201 frames each,usingthree references andtwo anchors. Uses read-only XTC and source/role/contact checks.|
|03_annotate_kinome_geometry.py|Geometry annotation of723 archived static poses. No experimental labels inferred. Eligibility remains blocked by original input QC.|
|04_compare_native_scores.py|Equal-pose,eight-target ranking and all declared cutoff/exposure/orientation/reference sensitivities.|
|05_prepare_matched_gnn_inputs.py|Build eight34-feature candidate graphs and eight zero-padded controls,preserving original29 features and every edge.|
|06_independent_md_verification.py|Independent raw-frame distance/exposure/clash replay and pairwise AUROC verification.|
|07_check_gnn_runtime.py|Guard incompatible optional native modules;check16 graphs and finite untrained CPUforward;zero optimizersteps.|
|08_write_complete_report.py|Write the initial all-workstream report. Final combined report comes fromstep13.|
|09_model_covalent_donor_tail.py|Fixed-backbone donor-tail geometry trial;rejected outputs remain explicitly failed.|
|10_check_covalent_tail_sterics.py|Check peptide-bond exclusions and deposited nonbonded neighbors forstep09;does not override its chemistry failure.|
|11_relax_terminal_oxygen.py|Bounded terminal-oxygen/Gly76 reconstruction independent oftarget labels.8RX0passes declared geometry;8R5Hdoesnot.|
|12_run_chemical_tail_sensitivity.py|Rebuild/check the chemical reference byte-identically and compare/replay completed secondary MD results. --extract reprocesses existingMD only.|
|13_finalize_all_findings.py|Independent global tail-neighbor replay;combine primary and secondary findings and write final Word/Markdown report.|
|14_check_package.py|Verify immutable source/configuration,two141-pose result sets,16 GNNinputs,723 kinomeannotation and numerical receipts.|

Step12's secondary extraction scripts,helpers,reference packet,protocol and full results are under chemical_tail_sensitivity/. The existing reference is byte-reproducible from the saved local tail model. Invalid initial tail attempts remain in validation/tail_attempt0;they are not current approved references.

Cached MD results require their frozen extraction source/configuration/reference manifest and packet hashes. A changed scientific method requires a new version,not refreshing the old manifest to accept stale results. Scripts00/14 audit current artifacts;steps01,03–13 regenerate their declared derivatives. Archive outputs before a deliberate revised experiment. No script authorizes a production simulation or trained biological release.

Original read-only dependencies are the protac-md trajectories,role audits and completed full-frame source bindings;mechanistic_v3 supplies the old temporal base graphs/model definition. The installed subfolder is self-contained for generated artifacts but is intentionally bound to those existing datasets on this workstation. The chemistry-correct supplier reference and historical simulation remain separate.

Main outputs:results/primary_target_scores.tsv,results/native_metrics.tsv,results/native_controls.tsv,results/kinome_pose_annotations.tsv,results/kinome_target_annotations.tsv,gnn_inputs/,chemical_tail_sensitivity/results/,validation/scientific_status.json,FINDINGS.docx. Raw geometric shortlist files are labelled DIAGNOSTIC_ONLY;none is a validated binder/degrader prediction.
