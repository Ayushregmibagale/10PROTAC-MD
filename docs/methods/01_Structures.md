# Structure provenance and domain selection

The inherited KLIFS/AF2 structure collection is not a released production dataset.
The archived merged manifest requires review of exact KLIFS ID, chain, alternate
state, domain identity, mutations, missing pocket residues and AF2 confidence.
UniProt kinase-keyword membership alone does not establish a protein-kinase fold.
Do not silently deduplicate tandem kinase domains by UniProt accession.

Verified preparation references are deposited mmCIF files with SHA256 provenance:
7PI4/7QB (GSK215 native FAK-VHL ternary complex), 1XKK/FMM (lapatinib parent
ligand), and 6BFN/DL1 (IRAK1 with JH-I-25 parent ligand). Parent-ligand references
are not identical to the conjugated head of the purchased PROTAC.

5JC7 is chicken MDA5 and must not be used as a FAK reference. 6EGA is an IRAK4
complex with a type II inhibitor and is not the documented IRAK1/JH-I-25 reference.
The originals remain archived; the new preparation code does not select them.

New placement requires explicit `protac,gene,uniprot_id,structure_path,chain,
res_start,res_end` TSV fields. Whole full-length AF2 chains and motif-derived
windows cannot automatically substitute for verified domain boundaries.
DFG state is a descriptor, not a universal ligand-compatibility rule.
