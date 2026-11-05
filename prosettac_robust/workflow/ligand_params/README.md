# Verified ligand parameter states

Each compound/state folder contains atom mappings, execution provenance and a
parameter_validation.json receipt. The selected protein/ligand/water protocol is
AMBER99SB-ILDN/GAFF2/TIP3P at 300 K. Use 06_prepare_md.py to build and verify a
combined system; use 07_run_md_stage.py for staged, receipt-bound execution.

The chosen ligand topologies are the verified *_GMX.itp GAFF2 outputs. ACPYPE's
other generated variants, template MDPs and rungmx.sh remain archived artifacts;
they are not the combined-system protocol. File presence alone does not certify
parameters: the validation receipt and atom mapping must match.

neutral_reference denotes the supplier neutral graph, not a conclusion about the
populated physiological state. The GSK215 +1 state is a separate verified set.
The provisional pH assessment and sensitivity limitations are documented under
molecules/protonation_assessment_pH7p4/. No state populations are asserted.
