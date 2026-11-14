# Experimental-label and chemistry review for the matched GNN ablation

The eight training/evaluation labels remain exactly the user's confirmed six-degrader/two-nondegrader map. This read-only review is not used to change labels,features,thresholds or training. It is separate from the new matched ablation.

The user-confirmed LF/PS conditions are RPE1,500nM,12-hour treatment. They differ from the published cell-line endpoint reviewed previously. The PS heavy channel is newly synthesized protein;neither PS channel alone establishes a degradation-rate constant. Authoritative LF processing/version and independent replication provenance remain unresolved. Named_new is a saved version,not proof that it is the final as-run pipeline.

## Version-explicit LF SJF8240-versus-DMSO view

| Target | User label | log2FC | BH-adjustedP | Interpretation |
|---|---|---:|---:|---|
| ABL1 | 0 | -0.084 | 0.6921 | No statistically established change in this table;not equivalence |
| AXL | 0 | 0.230 | 0.04139 | Observed increase |
| DDR2 | 1 | -2.126 | 1.055e-06 | Observed decrease |
| EPHA2 | 1 | -0.157 | 0.1338 | No statistically established change in this table;not equivalence |
| EPHB2 | 1 | -1.106 | 1.826e-05 | Observed decrease |
| MAPK14_full | 1 | -0.859 | 0.0008413 | Observed decrease |
| MET | 1 | -1.080 | 1.826e-05 | Observed decrease |
| RIPK2_full | 1 | -1.325 | 0.009906 | Observed decrease |

Five of the six user-positive labels show significant decreases in this named LF table. EPHA2 is inconclusive in this RPE1 table;this does not override the user's positive label or establish a non-degrader. ABL1 is inconclusive near zero;AXL increases versus DMSO. These observations are supportive/contextual evidence,not independent confirmation of each exact binary label or proof of direct ubiquitination. A nonsignificant result is not evidence of equivalence or no degradation.

The complete saved-version and contrast table is retained separately. In particular,AXL can decrease versus foretinib while increasing versus DMSO;the comparator matters. LF versions differ and are not averaged or cherry-picked to match the model. PS light/heavy channels are reported separately. No R pipeline was rerun,no multi-gigabyte raw spectra were reprocessed,andno new wet-lab result was generated. The earlier review's source summaries were checked against the named source files;LFvalues reproduce to1e-12.

The supplied folders also contain later PELSA datasets,but their generic processing READMEs do not establish binding for these eight targets. This review does not claim direct binding validation or analyze those newer assays into model labels. No messages or external submissions were sent.

## Chemistry and interpretation

The retained computational ligand has a historical connectivity/coordinate-stereochemistry provenance that differs from supplier-matched SJF8240. The new ablation preserves that actual historical model and all old results;it cannot validate a supplier-SJF8240 classifier merely by improving AUROC. Geometry repair does not change molecular connectivity. Historical and supplier-matched branches remain separate. Chemical identity audit:../native_catalytic_validation_v2/FINDINGS.md and../../SJF8240_kinome_screen_analysis/repair_v3/chemical_identity_audit.json.

The independent units are eight targets with only two negatives. Assay context,protein identity and chemical provenance limit biological interpretation even if the controlled input ablation improves ranking. Experimental sources were not used to tune the ablation or relabel EPHA2.
