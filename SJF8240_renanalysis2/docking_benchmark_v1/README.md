# Controlled docking benchmark and source audit

The native-recovery benchmark is not yet complete. The first completed stage is the comparison of published/source algorithms, pose budgets and error handling. See code_audit/FINDINGS.md. The upstream builder is publicly available and pinned to a174357135d599c7274d24779e0662e15461b0ef; its MIT license is preserved. The complete2022 modified docking/CRL4A driver was not found in that repository. A reconstructed method must be called an adaptation.

Numbered audit scripts:

1. code_audit/01_probe_upstream.py: run unmodified upstream example and invalid-input probes in isolated folders using the existing docking runtime. Does not run Rosetta docking. Records default example success and actual error-handling weaknesses.
2. code_audit/02_inventory_and_report.py: pin source hashes, compare published pose budgets with source defaults and four saved historical target summaries; report limitations. Summary counts are not verified unique model counts when raw directories are absent.

Our existing regression suite was run read-only with PYTHONDONTWRITEBYTECODE=1 and cacheprovider disabled:141passed. The preserved docking overlay lacks pytest; the suite was run in the analysis Python3.13 runtime, separately from upstream example probes in the docking Python3.10 overlay. This is software verification,not scientific recovery validation.

Prior known-structure recovery experiments exist for GSK215/7PI4: one seed failed the exploratory criterion and a subsequent adaptive seed passed. Those are supplementary,seed-dependent controls,not the proposed equal-input VHL–MZ1–BRD4 comparison. No new docking outcome is claimed here. No new SJF8240 panel,MD or GNN run has been launched.
