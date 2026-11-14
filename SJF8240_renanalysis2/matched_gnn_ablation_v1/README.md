# Matched catalytic-feature GNN ablation

Completed exploratory24-fold CPU training comparison,not a supplier-SJF8240 classifier release. Read FINDINGS.md/FINDINGS.docx and EXPERIMENTAL_EVIDENCE_REVIEW.md.

Primary matched-pair ranking AUROC:zero-padded control0.750,native-feature model0.833. Two pair rankings improve,one worsens,nine are unchanged. Both models have sensitivity1.000,specificity0.000,balancedaccuracy0.500 atfixed0.5. Ranking improves byone netcomparison;held-out non-degrader classification remains unresolved.

Both models use34-input temporal GAT architecture,identical folds,seeds andinitialweights,300epochs andtheoriginal loss/optimizer settings.24 checkpoints and48 held-out scores have independent CPUreplay. No threshold,seed,epoch orendpoint selected fromtheoutcome. Allold models anddata remain unchanged.

Weights are physically saved locally in results/<variant>/<pair>.pt;source-bound receipts andpredictiontables accompany them. Software checks do not establish biological validation. There are onlyeight independent labelledtargets andtwo negatives. UnknownSIK3/RPS6KB2 are excluded. Parallel read-only experimental review didnotchange labels.
