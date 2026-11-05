# Data retention and GitHub scope

Full raw output retention is the user's planning preference. The approximately
7.24 TB estimate is for 112 pairs at 200 global poses, 50 local poses and 100
conformers, extrapolated from a small recorded GSK215 diagnostic. It excludes MD,
retries, backups and reserve. No data deletion or remote campaign has occurred.

GitHub contains source, documentation, exact molecular identities, compact
reference/control structures, verified ligand parameters and a selected set of
small diagnostic receipts in validation_evidence/. The complete local project
retains audit_baselines/, audits/, benchmarks/, engine_audit/, domain_structures/,
prepared_screen/, legacy prosettac_inputs/ and structure_provenance/. Excluding
these paths from Git neither removes them nor authorizes later deletion.

The selected receipts preserve the paths/hashes from their actual historical
execution. They do not assert that every underlying large artifact is bundled
with the Git checkout. Before publication, deposit the finalized input inventory,
source metadata, raw outputs and immutable attempt receipts in the allocated
research data store and record its persistent accession. There is no invented
public data accession or assumed Ernie quota here.

The optional lossless pose archive restores complete original PDB bytes with
SHA256 verification. It does not delete the originals and is not a mandatory
replacement for full raw retention.

PyRosetta binaries/databases are never redistributed. Rosetta parameter helper
source is obtained directly from the official pinned repository; the project
records its file hashes and does not relicense it as LondonLab MIT code.

The separate 161-candidate pocket comparison would project approximately 10.41 TB
at those same budgets. It is not an expansion of the accepted 7.24-TB campaign
scope, and no jobs have been launched from it. Both input versions and exclusions
are retained for method comparison and cohort review.

The repaired runtime retains even zero-conformer markers and virtual compatibility
references. Earlier historical runs may lack those empty files because the
inherited upstream code removed them; their sampling ledgers remain intact.
The updated native-frame 100-conformer pilot projects ~3.37 TB for the original
112-pair model. It is a small conditional sample and does not reduce the accepted
~7.24-TB planning allowance or establish an allocated data store.
