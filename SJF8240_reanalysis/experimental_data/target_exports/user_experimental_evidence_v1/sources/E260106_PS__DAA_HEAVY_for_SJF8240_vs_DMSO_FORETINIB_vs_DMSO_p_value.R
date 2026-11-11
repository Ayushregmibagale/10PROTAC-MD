#### DAA 2 ####
# Proteomics Analysis - LIMMA (DMSO / Foretinib / SJF8240) - HEAVY CHANNEL
# + Selective degraded / increased plot with custom colors
#
# Outputs:
#  1) LIMMA_HEAVY_results_all_contrasts.csv.gz
#  2) LIMMA_HEAVY_SJF8240_vs_Foretinib.csv.gz
#  3) Selective_Degraded_Increased_HEAVY_labeled.png
#  4) Selective_Degraded_Increased_table_HEAVY_scored.csv.gz

library(data.table)
library(ggplot2)
library(limma)
library(ggrepel)

cat("🚀 LIMMA Analysis - HEAVY CHANNEL (DMSO / Foretinib / SJF8240)\n")
cat("=============================================================\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# 🎛️  THRESHOLDS — edit here, changes propagate everywhere
#     NOTE: using RAW p-value (P.Value), NOT FDR (adj_P_Val)
# ══════════════════════════════════════════════════════════════════════════════
pval_cut     <- 0.05   # raw p-value threshold
fc_cut       <- 0.3    # |log2FC| threshold

cat("🎛️  Thresholds: raw p-value <", pval_cut, "| |logFC| >=", fc_cut, "\n\n")

# ══════════════════════════════════════════════════════════════════════════════
#### 📂 PART 1: LIMMA DIFFERENTIAL EXPRESSION ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════
cat("📊 PART 1: LIMMA DIFFERENTIAL EXPRESSION ANALYSIS - HEAVY CHANNEL\n")
cat("=================================================================\n")

cat("📂 Loading protein quantitation data...\n")
DT <- fread("protein_quant.csv.gz")
cat("   Rows loaded:", nrow(DT), "| Columns:", ncol(DT), "\n")

req_cols <- c("Protein_group", "Sample", "Condition", "Replicate",
              "LFQ_H", "N_precursors_H", "Genes", "Protein_description")
missing <- setdiff(req_cols, names(DT))
if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))
cat("✅ All required columns present\n")

DT <- DT[Condition %in% c("DMSO", "Foretinib", "SJF8240")]
cat("✅ Conditions retained:", paste(unique(DT$Condition), collapse = ", "), "\n")
cat("   Samples per condition:\n")
print(DT[, .(N_rows = .N, N_replicates = uniqueN(Replicate)), by = Condition])

DT[, Genes   := fifelse(is.na(Genes) | !nzchar(Genes), "UNMAPPED_GENE", Genes)]
DT[, Feature := paste0(Genes, "|", Protein_group)]

DT <- DT[!is.na(LFQ_H) & LFQ_H > 0]
cat("✅ Filtered to valid LFQ_H > 0\n")
cat("   Features:", uniqueN(DT$Feature), "\n")

DT[, log2_quant_H := log2(LFQ_H)]

DT <- DT[!is.na(N_precursors_H) & N_precursors_H >= 2]
cat("✅ Filtered N_precursors_H >= 2\n")
cat("   Features:", uniqueN(DT$Feature), "\n")

keepers <- DT[, .(n_reps = uniqueN(Replicate)), by = .(Feature, Condition)][
  n_reps >= 2
][, .(n_conds = uniqueN(Condition)), by = Feature][
  n_conds >= 2, Feature
]
DT <- DT[Feature %in% keepers]
cat("✅ Filtered: >= 2 replicates in >= 2 conditions\n")
cat("   Features:", uniqueN(DT$Feature), "\n")

target_median <- DT[, .(med = median(log2_quant_H, na.rm = TRUE)), by = Sample][, median(med)]
DT[, log2_quant_H_norm := log2_quant_H - median(log2_quant_H, na.rm = TRUE) + target_median,
   by = Sample]
cat("✅ Sample-wise median centering applied\n")
cat("   Per-sample median check (should all be ~", round(target_median, 3), "):\n")
print(DT[, .(median_norm = round(median(log2_quant_H_norm, na.rm = TRUE), 4)),
         by = Sample][order(Sample)])

wDT <- dcast(DT, Feature ~ Sample, value.var = "log2_quant_H_norm", fill = NA_real_)
cat("✅ Wide format:", nrow(wDT), "features x", ncol(wDT) - 1, "samples\n")

sample_map <- unique(DT[, .(Sample, Condition, Replicate)])
sample_map[, Condition := factor(Condition, levels = c("DMSO", "Foretinib", "SJF8240"))]

cat("\n🔍 Sample map:\n")
print(sample_map[order(Condition, Replicate, Sample)])

design_mat <- model.matrix(~ 0 + Condition, data = sample_map)
colnames(design_mat) <- gsub("^Condition", "", colnames(design_mat))
rownames(design_mat) <- sample_map$Sample

cat("\nDesign matrix:\n")
print(design_mat)

common_samples <- intersect(rownames(design_mat), colnames(wDT))
cat("\n✅ Common samples between design and data:", length(common_samples), "\n")
if (length(common_samples) < 2) stop("❌ Not enough common samples")

design_mat  <- design_mat[common_samples, , drop = FALSE]
wDT_subset  <- wDT[, c("Feature", common_samples), with = FALSE]
features    <- wDT_subset$Feature
expr_matrix <- as.matrix(wDT_subset[, -"Feature", with = FALSE])
rownames(expr_matrix) <- features
expr_matrix <- expr_matrix[rowSums(!is.na(expr_matrix)) > 0, , drop = FALSE]

cat("✅ Expression matrix:", nrow(expr_matrix), "rows x", ncol(expr_matrix), "cols\n")
cat("   Missing value rate:", round(mean(is.na(expr_matrix)) * 100, 1), "%\n")

if (!identical(rownames(design_mat), colnames(expr_matrix))) {
  design_mat <- design_mat[colnames(expr_matrix), , drop = FALSE]
}

fit <- lmFit(expr_matrix, design_mat)
cat("✅ Linear models fitted\n")

contrast.matrix <- makeContrasts(
  Foretinib_vs_DMSO    = Foretinib - DMSO,
  SJF8240_vs_DMSO      = SJF8240   - DMSO,
  SJF8240_vs_Foretinib = SJF8240   - Foretinib,
  levels = design_mat
)
cat("✅ Contrast matrix:\n")
print(contrast.matrix)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
cat("✅ eBayes moderation applied\n")
cat("   Prior df:", round(fit2$df.prior, 2),
    "| s2.prior:", round(fit2$s2.prior, 4), "\n")

# ── Calculate SE ──────────────────────────────────────────────────────────────
contrast_names <- colnames(contrast.matrix)

results_list <- lapply(contrast_names, function(contrast) {
  tt <- as.data.table(topTable(fit2, coef = contrast, number = Inf, sort.by = "none"))
  tt[, Contrast     := contrast]
  tt[, Feature      := rownames(expr_matrix)]
  tt[, SE_model     := sqrt(fit2$s2.post) * sqrt(fit2$cov.coefficients[contrast, contrast])]
  tt[, SE_empirical := fifelse(is.finite(t) & t != 0, abs(logFC / t), NA_real_)]
  tt[, SE           := pmax(SE_model, SE_empirical, na.rm = TRUE)]
  tt[is.na(SE) | SE <= 0 | !is.finite(SE), SE := 0.1]
  tt[, .(Feature, logFC, AveExpr, SE, SE_model, SE_empirical, t, P.Value, adj.P.Val, Contrast)]
})

results_H <- rbindlist(results_list)
ann       <- unique(DT[, .(Feature, Genes, Protein_group, Protein_description)])
results_H <- merge(results_H, ann, by = "Feature", all.x = TRUE)
setnames(results_H, "adj.P.Val", "adj_P_Val")

# ── Per-contrast diagnostics ──────────────────────────────────────────────────
cat("\n📊 Per-contrast summary:\n")
cat(strrep("\u2500", 70), "\n")
for (ct in contrast_names) {
  sub <- results_H[Contrast == ct]
  cat(sprintf("\n  \u25b6 %s\n", ct))
  cat(sprintf("    Features tested   : %d\n", nrow(sub)))
  cat(sprintf("    logFC range       : [%.3f, %.3f]\n",
              min(sub$logFC, na.rm = TRUE), max(sub$logFC, na.rm = TRUE)))
  cat(sprintf("    P.Value range     : [%.2e, %.3f]\n",
              min(sub$P.Value, na.rm = TRUE), max(sub$P.Value, na.rm = TRUE)))
  cat(sprintf("    adj_P_Val range   : [%.2e, %.3f]\n",
              min(sub$adj_P_Val, na.rm = TRUE), max(sub$adj_P_Val, na.rm = TRUE)))
  cat(sprintf("    p < 0.001         : %d\n", sub[P.Value < 0.001, .N]))
  cat(sprintf("    p < 0.01          : %d\n", sub[P.Value < 0.010, .N]))
  cat(sprintf("    p < 0.05          : %d\n", sub[P.Value < 0.050, .N]))
  cat(sprintf("    p < 0.05 & |FC|>=%.1f : %d\n", fc_cut,
              sub[P.Value < 0.05 & abs(logFC) >= fc_cut, .N]))
  cat(sprintf("    FDR < 0.05        : %d\n", sub[adj_P_Val < 0.05, .N]))
  cat(sprintf("    FDR < 0.20        : %d\n", sub[adj_P_Val < 0.20, .N]))
  cat(sprintf("    Top 10 by P.Value:\n"))
  print(sub[order(P.Value)][1:min(10, .N),
                            .(Genes,
                              logFC     = round(logFC,     3),
                              P.Value   = signif(P.Value,  3),
                              adj_P_Val = signif(adj_P_Val, 3))])
}
cat(strrep("\u2500", 70), "\n")

fwrite(results_H, "LIMMA_HEAVY_results_all_contrasts.csv.gz")
cat("\n\u2705 Saved: LIMMA_HEAVY_results_all_contrasts.csv.gz\n")

fwrite(results_H[Contrast == "SJF8240_vs_Foretinib"],
       "LIMMA_HEAVY_SJF8240_vs_Foretinib.csv.gz")
cat("\u2705 Saved: LIMMA_HEAVY_SJF8240_vs_Foretinib.csv.gz\n")

# ══════════════════════════════════════════════════════════════════════════════
#### 📈 PART 2: SELECTIVE DEGRADED / INCREASED PLOT
# ══════════════════════════════════════════════════════════════════════════════
cat("\n\U0001f4c8 PART 2: Selective degraded / increased plot\n")
cat("===============================================\n")
cat("\U0001f3db\ufe0f  Using RAW p-value <", pval_cut, "| |logFC| >=", fc_cut, "\n")

# A: other drug "not affected" if (p>=pval_cut) OR (|logFC|<fc_cut)
# B: stricter: other drug must be (p>=pval_cut) AND (|logFC|<fc_cut)
other_not_affected_mode <- "A"

wide_dt <- dcast(
  results_H[Contrast %in% c("Foretinib_vs_DMSO", "SJF8240_vs_DMSO")],
  Feature + Genes + Protein_group + Protein_description ~ Contrast,
  value.var = c("logFC", "P.Value", "adj_P_Val")
)

setnames(wide_dt,
         old = c("logFC_Foretinib_vs_DMSO", "logFC_SJF8240_vs_DMSO",
                 "P.Value_Foretinib_vs_DMSO",   "P.Value_SJF8240_vs_DMSO",
                 "adj_P_Val_Foretinib_vs_DMSO", "adj_P_Val_SJF8240_vs_DMSO"),
         new = c("logFC_Fore", "logFC_SJF",
                 "P_Fore",     "P_SJF",
                 "FDR_Fore",   "FDR_SJF")
)

# ── Significance flags using RAW p-value ─────────────────────────────────────
wide_dt[, Fore_down := !is.na(P_Fore) & P_Fore < pval_cut & logFC_Fore <= -fc_cut]
wide_dt[, Fore_up   := !is.na(P_Fore) & P_Fore < pval_cut & logFC_Fore >=  fc_cut]
wide_dt[, SJF_down  := !is.na(P_SJF)  & P_SJF  < pval_cut & logFC_SJF  <= -fc_cut]
wide_dt[, SJF_up    := !is.na(P_SJF)  & P_SJF  < pval_cut & logFC_SJF  >=  fc_cut]

if (other_not_affected_mode == "A") {
  wide_dt[, Fore_not_affected := is.na(P_Fore) | P_Fore >= pval_cut | abs(logFC_Fore) < fc_cut]
  wide_dt[, SJF_not_affected  := is.na(P_SJF)  | P_SJF  >= pval_cut | abs(logFC_SJF)  < fc_cut]
} else {
  wide_dt[, Fore_not_affected := is.na(P_Fore) | (P_Fore >= pval_cut & abs(logFC_Fore) < fc_cut)]
  wide_dt[, SJF_not_affected  := is.na(P_SJF)  | (P_SJF  >= pval_cut & abs(logFC_SJF)  < fc_cut)]
}

# ── Diagnostics BEFORE categorisation ────────────────────────────────────────
cat("\n\U0001f52c Filter diagnostics (raw p <", pval_cut, "):\n")
cat(sprintf("   SJF_down  (p<%.2f & logFC<=-%.1f) : %d\n", pval_cut, fc_cut, wide_dt[SJF_down  == TRUE, .N]))
cat(sprintf("   SJF_up    (p<%.2f & logFC>= %.1f) : %d\n", pval_cut, fc_cut, wide_dt[SJF_up    == TRUE, .N]))
cat(sprintf("   Fore_down (p<%.2f & logFC<=-%.1f) : %d\n", pval_cut, fc_cut, wide_dt[Fore_down == TRUE, .N]))
cat(sprintf("   Fore_up   (p<%.2f & logFC>= %.1f) : %d\n", pval_cut, fc_cut, wide_dt[Fore_up   == TRUE, .N]))
cat(sprintf("   Fore_not_affected                 : %d\n", wide_dt[Fore_not_affected == TRUE, .N]))
cat(sprintf("   SJF_not_affected                  : %d\n", wide_dt[SJF_not_affected  == TRUE, .N]))
cat(sprintf("   SJF_down  & Fore_not_affected     : %d  <- SJF selective degraded\n",
            wide_dt[SJF_down == TRUE & Fore_not_affected == TRUE, .N]))
cat(sprintf("   SJF_up    & Fore_not_affected     : %d  <- SJF increased\n",
            wide_dt[SJF_up   == TRUE & Fore_not_affected == TRUE, .N]))

# Threshold sweep for SJF8240
cat("\n\U0001f4c9 SJF8240 vs DMSO — hits across p-value thresholds:\n")
for (pv in c(0.001, 0.01, 0.05, 0.10, 0.20)) {
  n_down <- wide_dt[!is.na(P_SJF) & P_SJF < pv & logFC_SJF <= -fc_cut, .N]
  n_up   <- wide_dt[!is.na(P_SJF) & P_SJF < pv & logFC_SJF >=  fc_cut, .N]
  cat(sprintf("   p < %.3f : down=%d  up=%d\n", pv, n_down, n_up))
}
cat("   (no FC filter) p < 0.05 :", wide_dt[!is.na(P_SJF) & P_SJF < 0.05, .N], "\n")

# ── Category assignment: "both" first, selective overwrites ──────────────────
wide_dt[, Category := "Other / not selective"]
wide_dt[SJF_down  & Fore_down,       Category := "Degraded by both"]
wide_dt[SJF_up    & Fore_up,         Category := "Increased by both"]
wide_dt[SJF_down  & Fore_not_affected, Category := "SJF selective degraded"]
wide_dt[SJF_up    & Fore_not_affected, Category := "SJF increased"]
wide_dt[Fore_down & SJF_not_affected,  Category := "Foretinib selective degraded"]
wide_dt[Fore_up   & SJF_not_affected,  Category := "Foretinib increased"]

cat("\n\U0001f4ca Category counts:\n")
print(wide_dt[, .N, by = Category][order(-N)])

# ── Scores (using raw p for scoring too) ─────────────────────────────────────
wide_dt[, P_SJF_safe  := pmax(P_SJF,  1e-300)]
wide_dt[, P_Fore_safe := pmax(P_Fore, 1e-300)]
wide_dt[, Score_SJF  := sign(logFC_SJF)  * abs(logFC_SJF)  * (-log10(P_SJF_safe))]
wide_dt[, Score_Fore := sign(logFC_Fore) * abs(logFC_Fore) * (-log10(P_Fore_safe))]
wide_dt[, SelectivityScore := Score_SJF - Score_Fore]

selective_cats <- c("SJF selective degraded", "SJF increased",
                    "Foretinib selective degraded", "Foretinib increased")
wide_dt[, SelectiveSignificant := Category %in% selective_cats]
wide_dt[, Label := fifelse(!is.na(Genes) & nzchar(Genes) & Genes != "UNMAPPED_GENE",
                           Genes, Feature)]

# ── Print hit lists for all non-grey categories ───────────────────────────────
all_cats <- c("SJF selective degraded", "SJF increased",
              "Foretinib selective degraded", "Foretinib increased",
              "Degraded by both", "Increased by both")

for (cat_name in all_cats) {
  cat(sprintf("\n=== %s ===\n", cat_name))
  hits <- wide_dt[Category == cat_name,
                  .(Genes,
                    logFC_SJF  = round(logFC_SJF,  3),
                    P_SJF      = signif(P_SJF,     3),
                    FDR_SJF    = signif(FDR_SJF,   3),
                    logFC_Fore = round(logFC_Fore,  3),
                    P_Fore     = signif(P_Fore,     3),
                    FDR_Fore   = signif(FDR_Fore,   3))][order(logFC_SJF)]
  if (nrow(hits) == 0) cat("  (none at current thresholds)\n") else print(hits)
}

# ── Custom colors ─────────────────────────────────────────────────────────────
custom_cols <- c(
  "SJF selective degraded"       = "#1f77b4",   # blue
  "SJF increased"                = "#ff7f0e",   # orange
  "Foretinib selective degraded" = "#2ca02c",   # green
  "Foretinib increased"          = "#d62728",   # red
  "Degraded by both"             = "#9467bd",   # purple
  "Increased by both"            = "#8c564b",   # brown
  "Other / not selective"        = "grey70"
)
wide_dt[, Category := factor(Category, levels = names(custom_cols))]

# ── Labels: Top 10 per category by |score| ───────────────────────────────────
top10_per_cat <- function(dt, cat_name) {
  sub <- dt[Category == cat_name & !is.na(Label)]
  if (nrow(sub) == 0) return(sub)
  # rank by absolute combined score (strongest effect + most significant first)
  sub[, abs_score := abs(Score_SJF) + abs(Score_Fore)]
  head(sub[order(-abs_score)], 10)
}

label_dt <- rbindlist(lapply(
  names(custom_cols)[names(custom_cols) != "grey70"],
  top10_per_cat,
  dt = wide_dt
), use.names = TRUE, fill = TRUE)
label_dt <- unique(label_dt, by = "Feature")

cat("\n\U0001f3f7\ufe0f  Gene labels to draw (top 10 per category):", nrow(label_dt), "\n")
cat("   Breakdown:\n")
print(label_dt[, .N, by = Category][order(-N)])

# ── Quadrant labels ───────────────────────────────────────────────────────────
x_range <- range(wide_dt$logFC_Fore, na.rm = TRUE)
y_range <- range(wide_dt$logFC_SJF,  na.rm = TRUE)
x_pad   <- diff(x_range) * 0.04
y_pad   <- diff(y_range) * 0.04

quadrant_labels <- data.frame(
  x     = c(x_range[1] + x_pad, x_range[2] - x_pad,
            x_range[1] + x_pad, x_range[2] - x_pad),
  y     = c(y_range[2] - y_pad, y_range[2] - y_pad,
            y_range[1] + y_pad, y_range[1] + y_pad),
  label = c(
    "Foretinib \u2193  SJF8240 \u2191\n(SJF increased)",
    "Foretinib \u2191  SJF8240 \u2191\n(Both increased)",
    "Foretinib \u2193  SJF8240 \u2193\n(Both degraded)",
    "Foretinib \u2191  SJF8240 \u2193\n(SJF selective degraded)"
  ),
  hjust = c(0, 1, 0, 1)
)

# ── Plot ──────────────────────────────────────────────────────────────────────
tryCatch({
  p <- ggplot(wide_dt, aes(x = logFC_Fore, y = logFC_SJF)) +
    geom_hline(yintercept = c(-fc_cut, 0, fc_cut), linetype = "solid", alpha = 0.35) +
    geom_vline(xintercept = c(-fc_cut, 0, fc_cut), linetype = "solid", alpha = 0.35) +
    # Quadrant annotation
    geom_text(
      data = quadrant_labels,
      aes(x = x, y = y, label = label, hjust = hjust),
      vjust = 1, size = 3.2, colour = "grey40",
      lineheight = 0.85, inherit.aes = FALSE
    ) +
    # Points — grey "other" drawn first, coloured on top
    geom_point(data = wide_dt[Category == "Other / not selective"],
               aes(color = Category), alpha = 0.5, size = 1.2) +
    geom_point(data = wide_dt[Category != "Other / not selective"],
               aes(color = Category), alpha = 0.85, size = 1.8) +
    scale_color_manual(values = custom_cols, drop = FALSE) +
    # Labels for ALL non-grey genes
    ggrepel::geom_text_repel(
      data               = label_dt,
      aes(label          = Label, color = Category),
      size               = 2.8,
      fontface           = "bold",
      max.overlaps       = Inf,
      box.padding        = 0.3,
      point.padding      = 0.2,
      min.segment.length = 0,
      segment.size       = 0.3,
      segment.alpha      = 0.6,
      show.legend        = FALSE
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right", panel.grid.minor = element_blank()) +
    labs(
      title    = "Selective degraded / increased proteins (HEAVY)",
      subtitle = paste0(
        "x = Foretinib vs DMSO logFC, y = SJF8240 vs DMSO logFC",
        " | raw p<", pval_cut, " & |logFC|>=", fc_cut,
        " | All significant genes labelled"
      ),
      x     = "Foretinib vs DMSO (log2 FC)",
      y     = "SJF8240 vs DMSO (log2 FC)",
      color = "Class"
    )
  
  ggsave("Selective_Degraded_Increased_HEAVY_labeled.png", p,
         width = 14, height = 10, dpi = 300, bg = "white")
  cat("\u2705 Plot saved: Selective_Degraded_Increased_HEAVY_labeled.png\n")
}, error = function(e) {
  cat("\u26a0\ufe0f Could not create plot:", e$message, "\n")
})

fwrite(wide_dt, "Selective_Degraded_Increased_table_HEAVY_scored.csv.gz")
cat("\u2705 Table saved: Selective_Degraded_Increased_table_HEAVY_scored.csv.gz\n")

cat("\n\U0001f3af DONE\n")
cat("======\n")