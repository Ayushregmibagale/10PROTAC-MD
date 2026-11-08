# Historical methods draft reconciliation

The original `protac-md/docs/paper_methods.md` remains unchanged. Its claims must not be copied into a submission without reconciliation. This table separates retained evidence from unsupported draft wording.

| Topic | Supported description / correction |
|---|---|
| Panel | 17 historical evaluation rows (five saved positive, twelve negative), plus two prediction-only targets; not eight confirmed degraders among nineteen. See PROTEIN_PANEL_19.md. |
| Docking | Retained adapter performs random rigid-body global sampling with PyRosetta. PatchDock folder names do not establish use of the PatchDock executable. Uniform search budgets and as-run seed receipts are not established. |
| MD duration / temperature | Compiled production settings request 20 ns, 2 fs and 300 K. All 174 retained production checkpoints reach 20,000 ps. A separately retained extended AXL processed trajectory requires an explicit window. Not a uniform 100 ns / 310 K campaign. |
| Cutoffs | 113 production cases use 1.0 nm and 61 use 1.6 nm; a current template cannot replace historical compiled evidence. |
| Minimization | 174 compiled inputs use steepest descent, a 50,000-step cap and 1,000 kJ mol−1 nm−1 tolerance. 173 logs state force convergence; ABL1 pose001 stops at machine precision with Fmax 2877.7249. Do not claim all converged. |
| Equilibration | 99 retained compiled/checkpoint sets establish NVT 200 ps and two NPT stages of 500 ps each at 300 K. Missing intermediate files for 75 cases do not prove omitted stages. |
| Force field | All 174 retained topology include inventories name AMBER99SB-ILDN and TIP3P. Include names alone do not freeze external parameter-file versions. |
| Constraints | LINCS is supported. Do not claim all-bond constraints solely from the algorithm field; source templates request H-bond constraints. |
| Contact fraction | Strict <5 Å simultaneous ligand contacts; historical VHL selector is ELOB. Corrected true-VHL/heavy analysis is separate. Not an 8 Å thermodynamic stability probability. |
| Dynamic graphs | Historical construction uses ligand-heavy to ELOB/target Cα distances, marginal contact-frequency products and a 0.01 edge criterion; it is not a protein Cα–Cα joint occupancy graph. Direct-coordinate imaging limitations are documented. |
| Static graphs | Historical VHL node role is ELOB; node/edge schemas and current expanded packets must not be described as exact frozen original training inputs. |
| GNN | Main 17-row evaluation is the E2 scalar dual encoder; two prediction-only rows use a different writer without E2 scalars. Declared source defaults do not establish dated as-run hyperparameters. |
| Weights | Two retained MGDB pretraining checkpoints are readable and finite. They are not the missing per-fold trained states of the main E2 GNN. |
| Identity / labels | CDK17-labelled prepared structure matches CDK14/O94921. Assay missingness/discordance and historical-vs-supplier ligand connectivity remain disclosed; no silent relabeling or structure substitution. |

Primary retained evidence: `existing_result_audit/retained_stage_evidence_v1`, `retained_checkpoint_times_v1`, `compiled_md_protocol_v1`, `retained_gnn_weights_v1`, and INTEGRATION_REVIEW_RESULTS.md. This corrects the methods record without claiming a new biological validation or repeating the experiment.
