# User-supplied RPE1 experimental evidence — 2026-11-08

## Conditions and scope

User supplied E250916_LF and E260106_PS under `/home/ayushregmibagale/SJF8240_LF_and_PS` and confirmed RPE1, 500 nM, 12-hour treatment for both experiments; PS heavy marks synthesis/newly synthesized proteins. The separate labeling-pulse duration and foretinib dose were not independently established in this pass. These are user-confirmed conditions, not a notebook-derived protocol. Three sample IDs per condition are retained; IDs alone do not establish independent biological replication.

This evidence differs from the published MDA-MB-231 endpoint reviewed earlier. Keep the experiments separate. Original experimental files, historical labels, structures, GNN predictions and MD remain unchanged. No label release, differential-model refit or new simulation occurred.

## What the saved tables show

LF `_new` has canonical-accession results for all 19 targets, including six missing from the published abundance table. The five historical positive targets have negative log2FC with BH-adjusted P below 0.05 in this named saved version. The separately saved SJF8240-versus-foretinib LF contrast also has negative effects and BH P below 0.05 for all five. This is supportive retrospective abundance evidence, conditional on validating the experimental processing/version; it is not a new independent predictive test or proof of direct PROTAC-mediated degradation.

The following is a version-explicit view, not selection of an authoritative result file. All versions and original protein groups remain in the machine-readable exports. Missing PS rows are not negative results. A small or nonsignificant change is not equivalence/no degradation.

| Target | LF `_new` log2FC / BH P | PS light log2FC / BH P | PS heavy log2FC / BH P |
|---|---|---|---|
| MAP4K5 | -0.206 / 0.044 | -0.072 / 0.707 | -0.141 / 0.51 |
| SLK | -0.032 / 0.721 | 0.018 / 0.922 | 0.009 / 0.967 |
| DDR2 | -2.126 / 1.06e-06 | -1.330 / 0.0537 | -1.258 / 6.42e-05 |
| RIPK2_full | -1.325 / 0.00991 | -0.831 / 0.168 | -0.524 / 0.465 |
| MAPK14_full | -0.859 / 0.000841 | -0.352 / 0.431 | -0.135 / 0.517 |
| MET | -1.080 / 1.83e-05 | -0.545 / 0.374 | -0.200 / 0.465 |
| EPHB2 | -1.106 / 1.83e-05 | -0.088 / 0.823 | 0.079 / 0.89 |
| ABL1 | -0.084 / 0.692 | -0.135 / 0.848 | -0.129 / 0.662 |
| AXL | 0.230 / 0.0414 | 0.520 / 0.411 | 0.326 / 0.465 |
| EPHA2 | -0.157 / 0.134 | -0.054 / 0.635 | 0.023 / 0.886 |
| CDK17 | -0.041 / 0.708 | -0.107 / 0.618 | -0.095 / 0.644 |
| CDK2 | 0.025 / 0.913 | -0.008 / 0.972 | 0.029 / 0.87 |
| CDK4 | -0.035 / 0.768 | 0.006 / 0.987 | -0.021 / 0.94 |
| CDK6 | 0.046 / 0.783 | -0.011 / 0.975 | -0.005 / 0.991 |
| DAPK1 | -0.040 / 0.818 | Not available | 0.673 / 0.499 |
| IRAK1 | -0.018 / 0.902 | -0.008 / 0.987 | 0.640 / 0.629 |
| TLK2 | -0.138 / 0.212 | Not available | 1.112 / 0.51 |
| SIK3 | -0.018 / 0.946 | Not available | Not available |
| RPS6KB2 | -0.402 / 0.419 | Not available | Not available |

## Important interpretation and processing findings

- MAP4K5 in LF `_new` has log2FC −0.206 and BH P 0.044, a modest abundance decrease; calling it simply unaffected loses information. It does not meet the explicitly displayed exploratory magnitude cutoff of −0.3. This does not resolve the earlier published endpoint disagreement by relabeling.
- ABL1 LF abundance changes little (−0.084, BH P 0.692). This describes the current RPE1 endpoint; it does not invalidate another cell line/assay.
- AXL is higher versus DMSO in LF `_new` (+0.230, BH P 0.041), while lower versus foretinib (−0.360, BH P 0.017). Control-relative direction must be reported. A negative SJF/foretinib contrast alone is not evidence of reduction below vehicle.
- RPS6KB2 has only one precursor in each SJF8240 LF `_new` sample according to retained quantification. Its weak evidence cannot be promoted to a validated negative. SIK3 and RPS6KB2 remain prediction-only in the historical GNN evaluation despite now having retrospective experimental measurements available.
- PS light has matched result rows for 15/19 and PS heavy for 17/19. Heavy is the newly synthesized pool under the user-confirmed design. Neither heavy intensity nor light intensity alone establishes a kinetic degradation constant; pulse design, incorporation and joint analysis matter. Do not collapse channels or equate reduced heavy signal with a causal synthesis-rate change.
- CDK17/Q00537 is quantified, but the retained CDK17-labelled computational structure is CDK14/O94921. CDK14 is not recovered as a matching abundance result in the inspected named tables. These data do not justify transferring the CDK17 label to CDK14 or modifying the old structure.
- Multiple LF result files differ materially, including direction for some targets. `_new` in a filename is not an as-run receipt. Earlier MAPK14 files contain distinct Q16539 and Q16539-2 rows; they must not be averaged or silently collapsed.
- `E260106_PS/LIMMA_results_DMSO_vs_SJF8240_Foretinib.csv.gz` is byte-identical to the older LF results file. It is not independent PS replication. Use channel-declared PS files and hashes.
- The LF differential script says at least two samples per condition but filters total rows >=4, which does not enforce that claim. PS scripts include raw P<0.05 and BH FDR<0.20 plotting variants. Retained BH values are exported, with an explicit descriptive 0.05 / −0.3 annotation; no threshold is selected to optimize historical GNN scores.

## How the review was done

Read-only Python audits hash 57 compact source/result/script inputs, extract accession/gene matches and retain gene-only alternatives and protein-group ambiguity flags. Protein quantification exports retain each matched sample and precursor count; support tables count positive measurements and those with at least two precursors. No original R script was executed. Multi-gigabyte raw DIA-NN reports and spectral libraries were inventoried by metadata rather than reprocessed or copied into Git. Five relevant original R scripts are frozen in the isolated evidence packet.

The first summary attempt stopped on distinct MAPK14 isoforms; its source and failure receipt are preserved in `user_experimental_summary_v1`. Version 2 keeps protein-group identity in uniqueness keys, preserves isoform rows, and rejects ambiguous display selection. This prevents a silent canonical-accession merge.

Reproduction: `existing_result_audit/review_tools/audit_user_experimental_evidence.py` and `summarize_user_experimental_evidence.py`; immutable outputs in `user_experimental_evidence_v1` and `user_experimental_summary_v2`. Full source hashes, all result versions, sample mapping and targeted quantification are retained. They are descriptive audits, not a new training input release.

## Remaining decisions

The experimental evidence gap is substantially narrower: local LF data cover all nineteen expected targets. Full scientific sign-off still needs the authoritative LF processing/version and replicate/protocol provenance, peptide-specific support where weak/ambiguous, the CDK17/CDK14 computational correspondence, and historical GNN replay/independent validation. A label matrix should distinguish decreased, inconclusive, unquantified and increased, with assay conditions and control contrast frozen before model use. No docking or MD repetition is required to perform those reviews.
