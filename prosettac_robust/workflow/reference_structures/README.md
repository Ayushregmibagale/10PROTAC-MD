# Reference structure scope

Use validated/ for the repaired experimental reference preparations. Parent
inhibitor structures are distinguished from exact conjugated PROTAC fragments.
Modeled modifications are separately stored in modeled_reference_candidates/.

PDB files directly in this directory are inherited historical inputs. Their
filenames do not establish the intended protein identity: 5JC7 is not the intended
FAK reference, and 6EGA must not be treated as an IRAK1 structure. The repaired
workflow does not use those legacy paths. Source identities, selections and
compound-specific mappings are documented in docs/REMEDIATION_LOG.md.
No inherited reference files have been deleted.
