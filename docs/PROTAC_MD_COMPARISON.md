# Reference workflow comparison

The user requested comparison with `/home/ayushregmibagale/protac-md`.
That project was inspected read-only. `PROTAC_MD_REFERENCE_SNAPSHOT.json`
records exact hashes of the inspected reference files. Its results and inputs
have not been modified or imported as validation labels for this project.

| Topic | Reference project evidence | Decision for 10PROTAC-MD and reason |
|---|---|---|
| Force fields | The executable setup/screen scripts use AMBER99SB-ILDN protein, GAFF2 ligand and TIP3P water. The paper methods describe CHARMM36m/CGenFF, another version and temperature. | Use the actual AMBER/GAFF2/TIP3P workflow and write methods from executable records. Do not reproduce the inconsistent paper description. |
| Temperature and duration | Screen script uses 300 K, 20 ns initial MD and 100 ns follow-up. | Retain 300 K and explicit configurable 20 ns preparation; duration is a sampling plan, not evidence of convergence. No production trajectory was launched in remediation. |
| Ligand parameters | SJF8240-specific GAFF2/AM1-BCC parameter files exist. | Generate exact stereochemical supplier molecules separately for all five compounds using GAFF2/AM1-BCC. Verify names, graph connectivity, charge, masses and output hashes. A shared ligand topology cannot represent five different molecules. |
| Placement | Older topology builder uses sequential heavy-atom reassignment and a fallback for hydrogens. | Map deposited ligand to the full stereochemical chemical graph, transfer coordinates by the mapping, add/minimize H with every heavy atom fixed, and reject mismatches. |
| Missing residues | `pdb2gmx -missing` can bypass missing atoms. | Model internal loops explicitly from cached UniProt sequence with seeded PDBFixer, retain observed frame, check peptide bonds and ligand clashes, and record every added residue. No `-missing`. These loops remain modeled hypotheses. |
| Warnings | Existing setup permits blanket `-maxwarn 2`. | GROMACS preprocessing must succeed without warning overrides. Failures are retained and repaired explicitly. |
| Solvent selection | Recent screen dynamically identifies solvent. Older scripts use a numeric group. | Select `SOL` by name and create explicit mapped heavy-atom groups plus a whole-system thermostat/COM group. |
| Equilibration | EM, NVT, restrained NPT, ligand-release NPT, production. | Preserve this staged release; continue velocities from checkpoints after NVT, use C-rescale for equilibration and Parrinello–Rahman for production. No velocity regeneration between stages. |
| Resume | Existing production extension can use `-noappend` and concatenate/delete trajectory pieces. | Resume from a checkpoint with GROMACS append/checksum validation. Input hashes and immutable per-stage intent/completion records prevent reuse of changed systems. |
| Pose selection | The screen uses seeded random selection rather than only top energies. | Preserve the rationale: an MD subset must represent the ensemble and include independent replicas; energy ranking alone is not a probability distribution. Selection and random seeds must be recorded before simulation. |
| Features and ML | Existing models include MD/E2 features and were developed on a small panel. | Geometry-only output is not a complete model input. No missing-feature zero filling or automatic reuse of trained weights across the replacement engine/chemotypes. Cross-compound experimental validation is still required. |

The replacement docking engine is a separately identified constrained-chemistry
Rosetta protocol. It must not be described as an unmodified upstream pRosettaC
run. Native integration and small perturbation tests check chemical/engine
integration; they do not establish prospective kinase selectivity.

`linker_geometry_yield` counts accepted ternaries against *all* attempted
protein poses. Contacts calculated only among successful, head-constrained
ternaries are conditional geometry metrics and cannot be interpreted as a
degradation probability. A finished search with no feasible pose is an explicit
sampling outcome, not an experimental nondegrader label.
