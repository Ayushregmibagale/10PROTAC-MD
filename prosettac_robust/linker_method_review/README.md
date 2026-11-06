# Experimental linker-method review

Status: EXPLORATORY_NOT_ADOPTED. Default workflow is unchanged.

`candidate_workflow` contains the frozen source variant. Center both head coordinate sets, disable ETKDG experimental torsion and basic-knowledge preferences, preserve all exact-head/stereochemical and full-complex checks, restore the original protein frame. `tests` contains three geometry controls. `paired_comparison.json` binds the comparison to 50 exactly identical local protein arrangements. `evidence` retains compact terminal receipts and focused controls. Full local run directories are deliberately Git-ignored.

The scored SJF8240/MET candidate control has 15 QC-valid poses from 50 local arrangements; the baseline has zero in those same arrangements. The GSK215 reference control passes its prespecified native-frame check. Neither result demonstrates prospective accuracy or biological degradation. Do not deploy this variant for the study until method choice, sampling controls and MD handoff are qualified and documented. The alternative retaining experimental torsion preferences requires a separate comparison.

MD control selection was deterministic filename order, independent of labels. Protein repair passed; minimization is in progress. Planned 20 ps production is only an integration diagnostic. Raw repaired/MD files remain private/local and are separate from released study outputs.
