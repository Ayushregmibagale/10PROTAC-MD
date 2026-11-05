# Repaired original pRosettaC baseline

The original protac-md project and the Ernie installation remain unchanged.
This isolated runtime preserves the existing global and local docking stages.
Exact source/patch provenance is in provenance.json, publication_repairs.patch
and runtime_manifest.json. LondonLab's MIT license is retained.

The adapter's former free-conformer rigid fit allowed the two bound ligand heads
to move by several angstroms. The repaired linker step constrains both complete
stereochemical heads during embedding and refinement; impossible geometries are
explicit sampling rejections. There is no free-conformer fallback. Official
Rosetta parameter tooling preserves explicit H, atom names, coordinates, and the
conformer library. Unknown residue types, worker failures and incomplete stages
are errors. Raw scores, including positive values, remain unshifted.

A native-frame GSK215 diagnostic completed the original local docking/linker/
scoring sequence with a pose passing the explicit 1.5 Å gross-overlap filter and
both head-core displacement checks. The receipt-compatible producer and the
bounded controller also recorded a genuine global no-valid-pose outcome.
These small integration controls do not demonstrate global recovery, sampling
convergence, binding affinity or campaign-wide accuracy.

The public repository omits Rosetta's separately licensed helper source and all
PyRosetta binaries/databases. Run fetch_rosetta_helpers.py to retrieve the exact
hash-verified helpers from the official Rosetta repository in an appropriately
licensed environment. The local working project already contains the verified
helpers used in the recorded controls. See vendor/SOURCE.json and vendor/NOTICE.md.

run.py is the receipt-compatible producer; execute.py invokes the isolated runtime.
Use scripts/01_setup/02_batch_prosettac.py --run --pilot with an explicit Python
interpreter and seed. Pilot budgets are enforced; full campaign mode requires
versioned scientific evidence, input/source hashes and a storage budget.
All raw diagnostics remain retained locally. No remote jobs are submitted here.

The linker search now declares 20 outer seed attempts and 20 internal ETKDG
embedding retries per attempt. Exhausting that budget is a sampling outcome.
All complete bound heads remain fixed, and ligand-rotamer orientation uses a ring
center and bonded heavy neighbors entirely inside the fixed VHL head. The
parameter generator verifies that this frame is invariant across conformers;
this preserves conformer selection without moving the heads during packing.

The controller freezes a campaign code/config package, copies a verified package
into each attempt, and runs both recorder and producer from that package. Later
working-tree edits cannot alter current or subsequent workers in that campaign.
Input drift and source-snapshot tampering fail verification. XML/flag settings are
covered as well as Python and JSON files. The original minimal source patch is
retained in publication_repairs.patch; added helpers are versioned separately,
and runtime_manifest.json records the assembled runtime.

Full raw retention includes empty failed-conformer markers and compatibility
reference files. Previously produced historical attempts are not rewritten.
