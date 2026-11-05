# Inherited warhead files

These files are retained historical inputs, not the repaired production inputs.
In particular, the inherited VS4718_warhead.sdf does not contain verified VS4718
chemistry and must not be used as the GSK215 head. The active workflow uses
compound-specific heads in reference_structures/validated/ or explicitly labeled
modeled_reference_candidates/, with complete supplier structures in molecules/.
Use the preparation and preflight scripts described in the project README.
The original artifacts remain intact for provenance.
