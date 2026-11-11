#### DAA 2 ####
# Proteomics Analysis - LIMMA (DMSO / Foretinib / SJF8240) - LIGHT CHANNEL
# + Selective degraded / increased plot with custom colors
#
# Outputs:
#  1) LIMMA_LIGHT_results_all_contrasts.csv.gz
#  2) LIMMA_LIGHT_SJF8240_vs_Foretinib.csv.gz
#  3) Selective_Degraded_Increased_LIGHT_labeled.png
#  4) Selective_Degraded_Increased_table_LIGHT_scored.csv.gz

library(data.table)
library(ggplot2)
library(limma)
library(ggrepel)

cat("🚀 LIMMA Analysis - LIGHT CHANNEL (DMSO / Foretinib / SJF8240)\n")
cat("=============================================================\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# 🎛️  THRESHOLDS — edit here, changes propagate everywhere
# ══════════════════════════════════════════════════════════════════════════════
alpha         <- 0.20   # FDR threshold for significance (relaxed from 0.05)
alpha_strict  <- 0.05   # kept for reference / annotation
fc_cut        <- 0.3    # |log2FC| threshold
topN_SJF      <- 50     # max gene labels per SJF class in plot

cat("🎛️  Thresholds: FDR <", alpha, "| |logFC| >=", fc_cut, "\n\n")

# ══════════════════════════════════════════════════════════════════════════════
#### 📂 PART 1: LIMMA DIFFERENTIAL EXPRESSION ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════
cat("📊 PART 1: LIMMA DIFFERENTIAL EXPRESSION ANALYSIS - LIGHT CHANNEL\n")
cat("=================================================================\n")

cat("📂 Loading protein quantitation data...\n")
DT <- fread("protein_quant.csv.gz")
cat("   Rows loaded:", nrow(DT), "| Columns:", ncol(DT), "\n")

req_cols <- c("Protein_group", "Sample", "Condition", "Replicate",
              "LFQ_L", "N_precursors_L", "Genes", "Protein_description")
missing <- setdiff(req_cols, names(DT))
if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))
cat("✅ All required columns present\n")

DT <- DT[Condition %in% c("DMSO", "Foretinib", "SJF8240")]
cat("✅ Conditions retained:", paste(unique(DT$Condition), collapse = ", "), "\n")
cat("   Samples per condition:\n")
print(DT[, .(N_rows = .N, N_replicates = uniqueN(Replicate)), by = Condition])

DT[, Genes := fifelse(is.na(Genes) | !nzchar(Genes), "UNMAPPED_GENE", Genes)]
DT[, Feature := paste0(Genes, "|", Protein_group)]

DT <- DT[!is.na(LFQ_L) & LFQ_L > 0]
cat("✅ Filtered to valid LFQ_L > 0\n")
cat("   Features:", uniqueN(DT$Feature), "\n")

DT[, log2_quant_L := log2(LFQ_L)]

DT <- DT[!is.na(N_precursors_L) & N_precursors_L >= 2]
cat("✅ Filtered N_precursors_L >= 2\n")
cat("   Features:", uniqueN(DT$Feature), "\n")

keepers <- DT[, .(n_reps = uniqueN(Replicate)), by = .(Feature, Condition)][
  n_reps >= 2
][, .(n_conds = uniqueN(Condition)), by = Feature][
  n_conds >= 2, Feature
]
DT <- DT[Feature %in% keepers]
cat("✅ Filtered: >= 2 replicates in >= 2 conditions\n")
cat("   Features:", uniqueN(DT$Feature), "\n")

target_median <- DT[, .(med = median(log2_quant_L, na.rm = TRUE)), by = Sample][, median(med)]
DT[, log2_quant_L_norm := log2_quant_L - median(log2_quant_L, na.rm = TRUE) + target_median, by = Sample]
cat("✅ Sample-wise median centering applied\n")
cat("   Per-sample median check (should all be ~", round(target_median, 3), "):\n")
print(DT[, .(median_norm = round(median(log2_quant_L_norm, na.rm = TRUE), 4)), by = Sample][order(Sample)])

wDT <- dcast(DT, Feature ~ Sample, value.var = "log2_quant_L_norm", fill = NA_real_)
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

results_L <- rbindlist(results_list)
ann       <- unique(DT[, .(Feature, Genes, Protein_group, Protein_description)])
results_L <- merge(results_L, ann, by = "Feature", all.x = TRUE)
setnames(results_L, "adj.P.Val", "adj_P_Val")

# ── Per-contrast diagnostics ──────────────────────────────────────────────────
cat("\n📊 Per-contrast summary:\n")
cat(strrep("─", 70), "\n")
for (ct in contrast_names) {
  sub <- results_L[Contrast == ct]
  cat(sprintf("\n  ▶ %s\n", ct))
  cat(sprintf("    Features tested : %d\n", nrow(sub)))
  cat(sprintf("    logFC range     : [%.3f, %.3f]\n",
              min(sub$logFC, na.rm = TRUE), max(sub$logFC, na.rm = TRUE)))
  cat(sprintf("    P.Value range   : [%.2e, %.3f]\n",
              min(sub$P.Value, na.rm = TRUE), max(sub$P.Value, na.rm = TRUE)))
  cat(sprintf("    adj_P_Val range : [%.2e, %.3f]\n",
              min(sub$adj_P_Val, na.rm = TRUE), max(sub$adj_P_Val, na.rm = TRUE)))
  cat(sprintf("    FDR < 0.05  : %d\n", sub[adj_P_Val < 0.05, .N]))
  cat(sprintf("    FDR < 0.10  : %d\n", sub[adj_P_Val < 0.10, .N]))
  cat(sprintf("    FDR < 0.20  : %d\n", sub[adj_P_Val < 0.20, .N]))
  cat(sprintf("    FDR < 0.20 & |logFC|>=%.1f : %d\n", fc_cut,
              sub[adj_P_Val < 0.20 & abs(logFC) >= fc_cut, .N]))
  cat(sprintf("    Top 10 by P.Value:\n"))
  print(sub[order(P.Value)][1:min(10, .N),
                            .(Genes, logFC = round(logFC, 3),
                              P.Value   = signif(P.Value,   3),
                              adj_P_Val = signif(adj_P_Val, 3))])
}
cat(strrep("─", 70), "\n")

fwrite(results_L, "LIMMA_LIGHT_results_all_contrasts.csv.gz")
cat("\n✅ Saved: LIMMA_LIGHT_results_all_contrasts.csv.gz\n")

fwrite(results_L[Contrast == "SJF8240_vs_Foretinib"],
       "LIMMA_LIGHT_SJF8240_vs_Foretinib.csv.gz")
cat("✅ Saved: LIMMA_LIGHT_SJF8240_vs_Foretinib.csv.gz\n")

# ══════════════════════════════════════════════════════════════════════════════
#### 📈 PART 2: SELECTIVE DEGRADED / INCREASED PLOT
# ══════════════════════════════════════════════════════════════════════════════
cat("\n📈 PART 2: Selective degraded / increased plot\n")
cat("===============================================\n")
cat("🎛️  Using FDR threshold:", alpha, "| |logFC| >=", fc_cut, "\n")

# A: other drug "not affected" if (FDR>=alpha) OR (|logFC|<fc_cut)
# B: stricter: other drug must be (FDR>=alpha) AND (|logFC|<fc_cut)
other_not_affected_mode <- "A"

wide_dt <- dcast(
  results_L[Contrast %in% c("Foretinib_vs_DMSO", "SJF8240_vs_DMSO")],
  Feature + Genes + Protein_group + Protein_description ~ Contrast,
  value.var = c("logFC", "adj_P_Val", "P.Value")
)

setnames(wide_dt,
         old = c("logFC_Foretinib_vs_DMSO", "logFC_SJF8240_vs_DMSO",
                 "adj_P_Val_Foretinib_vs_DMSO", "adj_P_Val_SJF8240_vs_DMSO",
                 "P.Value_Foretinib_vs_DMSO",   "P.Value_SJF8240_vs_DMSO"),
         new = c("logFC_Fore", "logFC_SJF",
                 "FDR_Fore",   "FDR_SJF",
                 "P_Fore",     "P_SJF")
)

# ── Significance flags ────────────────────────────────────────────────────────
wide_dt[, Fore_down := !is.na(FDR_Fore) & FDR_Fore < alpha & logFC_Fore <= -fc_cut]
wide_dt[, Fore_up   := !is.na(FDR_Fore) & FDR_Fore < alpha & logFC_Fore >=  fc_cut]
wide_dt[, SJF_down  := !is.na(FDR_SJF)  & FDR_SJF  < alpha & logFC_SJF  <= -fc_cut]
wide_dt[, SJF_up    := !is.na(FDR_SJF)  & FDR_SJF  < alpha & logFC_SJF  >=  fc_cut]

if (other_not_affected_mode == "A") {
  wide_dt[, Fore_not_affected := is.na(FDR_Fore) | FDR_Fore >= alpha | abs(logFC_Fore) < fc_cut]
  wide_dt[, SJF_not_affected  := is.na(FDR_SJF)  | FDR_SJF  >= alpha | abs(logFC_SJF)  < fc_cut]
} else {
  wide_dt[, Fore_not_affected := is.na(FDR_Fore) | (FDR_Fore >= alpha & abs(logFC_Fore) < fc_cut)]
  wide_dt[, SJF_not_affected  := is.na(FDR_SJF)  | (FDR_SJF  >= alpha & abs(logFC_SJF)  < fc_cut)]
}

# ── Diagnostics BEFORE categorisation ────────────────────────────────────────
cat("\n🔬 Filter diagnostics (FDR <", alpha, "):\n")
cat(sprintf("   SJF_down  (FDR<alpha & logFC<=-%.1f) : %d\n", fc_cut, wide_dt[SJF_down  == TRUE, .N]))
cat(sprintf("   SJF_up    (FDR<alpha & logFC>= %.1f) : %d\n", fc_cut, wide_dt[SJF_up    == TRUE, .N]))
cat(sprintf("   Fore_down (FDR<alpha & logFC<=-%.1f) : %d\n", fc_cut, wide_dt[Fore_down == TRUE, .N]))
cat(sprintf("   Fore_up   (FDR<alpha & logFC>= %.1f) : %d\n", fc_cut, wide_dt[Fore_up   == TRUE, .N]))
cat(sprintf("   Fore_not_affected                    : %d\n", wide_dt[Fore_not_affected == TRUE, .N]))
cat(sprintf("   SJF_not_affected                     : %d\n", wide_dt[SJF_not_affected  == TRUE, .N]))
cat(sprintf("   SJF_down  & Fore_not_affected        : %d  <- SJF selective degraded\n",
            wide_dt[SJF_down == TRUE & Fore_not_affected == TRUE, .N]))
cat(sprintf("   SJF_up    & Fore_not_affected        : %d  <- SJF increased\n",
            wide_dt[SJF_up   == TRUE & Fore_not_affected == TRUE, .N]))

# Looser sweep across thresholds
cat("\n📉 SJF8240 vs DMSO — hits at relaxed thresholds:\n")
for (a in c(0.05, 0.10, 0.20, 0.30, 0.50)) {
  n_down <- wide_dt[!is.na(FDR_SJF) & FDR_SJF < a & logFC_SJF <= -fc_cut, .N]
  n_up   <- wide_dt[!is.na(FDR_SJF) & FDR_SJF < a & logFC_SJF >=  fc_cut, .N]
  cat(sprintf("   FDR < %.2f : down=%d  up=%d\n", a, n_down, n_up))
}
cat("   (no FC filter) hits at FDR<0.20:",
    wide_dt[!is.na(FDR_SJF) & FDR_SJF < 0.20, .N], "\n")

# ── Category assignment: "both" first, selective overwrites ──────────────────
wide_dt[, Category := "Other / not selective"]
wide_dt[SJF_down  & Fore_down, Category := "Degraded by both"]
wide_dt[SJF_up    & Fore_up,   Category := "Increased by both"]
wide_dt[SJF_down  & Fore_not_affected, Category := "SJF selective degraded"]
wide_dt[SJF_up    & Fore_not_affected, Category := "SJF increased"]
wide_dt[Fore_down & SJF_not_affected,  Category := "Foretinib selective degraded"]
wide_dt[Fore_up   & SJF_not_affected,  Category := "Foretinib increased"]

cat("\n📊 Category counts:\n")
print(wide_dt[, .N, by = Category][order(-N)])

# ── Scores ────────────────────────────────────────────────────────────────────
wide_dt[, FDR_SJF_safe  := pmax(FDR_SJF,  1e-300)]
wide_dt[, FDR_Fore_safe := pmax(FDR_Fore, 1e-300)]
wide_dt[, Score_SJF  := sign(logFC_SJF)  * abs(logFC_SJF)  * (-log10(FDR_SJF_safe))]
wide_dt[, Score_Fore := sign(logFC_Fore) * abs(logFC_Fore) * (-log10(FDR_Fore_safe))]
wide_dt[, SelectivityScore := Score_SJF - Score_Fore]

selective_cats <- c("SJF selective degraded", "SJF increased",
                    "Foretinib selective degraded", "Foretinib increased")
wide_dt[, SelectiveSignificant := Category %in% selective_cats]
wide_dt[, Label := fifelse(!is.na(Genes) & nzchar(Genes) & Genes != "UNMAPPED_GENE",
                           Genes, Feature)]

# ── Print SJF selective hit lists ─────────────────────────────────────────────
cat("\n=== SJF selective degraded ===\n")
sjf_deg <- wide_dt[Category == "SJF selective degraded",
                   .(Genes,
                     logFC_SJF  = round(logFC_SJF,  3),
                     FDR_SJF    = signif(FDR_SJF,   3),
                     logFC_Fore = round(logFC_Fore,  3),
                     FDR_Fore   = signif(FDR_Fore,   3))][order(logFC_SJF)]
if (nrow(sjf_deg) == 0) cat("  (none at current thresholds)\n") else print(sjf_deg)

cat("\n=== SJF increased ===\n")
sjf_inc <- wide_dt[Category == "SJF increased",
                   .(Genes,
                     logFC_SJF  = round(logFC_SJF,  3),
                     FDR_SJF    = signif(FDR_SJF,   3),
                     logFC_Fore = round(logFC_Fore,  3),
                     FDR_Fore   = signif(FDR_Fore,   3))][order(-logFC_SJF)]
if (nrow(sjf_inc) == 0) cat("  (none at current thresholds)\n") else print(sjf_inc)

# ── Custom colors ─────────────────────────────────────────────────────────────
custom_cols <- c(
  "SJF selective degraded"       = "#1f77b4",
  "SJF increased"                = "#ff7f0e",
  "Foretinib selective degraded" = "#2ca02c",
  "Foretinib increased"          = "#d62728",
  "Degraded by both"             = "#9467bd",
  "Increased by both"            = "#8c564b",
  "Other / not selective"        = "grey70"
)
wide_dt[, Category := factor(Category, levels = names(custom_cols))]

# ── Labels: SJF only ──────────────────────────────────────────────────────────
build_label_dt <- function(dt, cat, score_col, ascending = TRUE, n) {
  sub <- dt[Category == cat]
  if (nrow(sub) == 0) return(sub)
  if (ascending) sub <- sub[order(get(score_col))] else sub <- sub[order(-get(score_col))]
  head(sub, n)
}

label_dt <- rbindlist(list(
  build_label_dt(wide_dt, "SJF selective degraded", "Score_SJF", ascending = TRUE,  n = topN_SJF),
  build_label_dt(wide_dt, "SJF increased",          "Score_SJF", ascending = FALSE, n = topN_SJF)
), use.names = TRUE, fill = TRUE)
label_dt <- unique(label_dt, by = "Feature")
cat("\n🏷️  Labels to draw:", nrow(label_dt), "\n")

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
    geom_text(
      data = quadrant_labels,
      aes(x = x, y = y, label = label, hjust = hjust),
      vjust = 1, size = 3.2, colour = "grey40",
      lineheight = 0.85, inherit.aes = FALSE
    ) +
    geom_point(aes(color = Category), alpha = 0.75, size = 1.3) +
    scale_color_manual(values = custom_cols, drop = FALSE) +
    ggrepel::geom_text_repel(
      data      = label_dt,
      aes(label = Label, color = Category),
      size      = 3.2,
      fontface  = "bold",
      max.overlaps     = Inf,
      box.padding      = 0.35,
      point.padding    = 0.25,
      min.segment.length = 0,
      show.legend      = FALSE
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right", panel.grid.minor = element_blank()) +
    labs(
      title    = "Selective degraded / increased proteins (LIGHT)",
      subtitle = paste0(
        "x = Foretinib vs DMSO logFC, y = SJF8240 vs DMSO logFC",
        " | FDR<", alpha, " & |logFC|>=", fc_cut,
        " | Bold labels = SJF8240 selective"
      ),
      x     = "Foretinib vs DMSO (log2 FC)",
      y     = "SJF8240 vs DMSO (log2 FC)",
      color = "Class"
    )
  
  ggsave("Selective_Degraded_Increased_LIGHT_labeled.png", p,
         width = 12, height = 9, dpi = 300, bg = "white")
  cat("✅ Plot saved: Selective_Degraded_Increased_LIGHT_labeled.png\n")
}, error = function(e) {
  cat("⚠️ Could not create plot:", e$message, "\n")
})

fwrite(wide_dt, "Selective_Degraded_Increased_table_LIGHT_scored.csv.gz")
cat("✅ Table saved: Selective_Degraded_Increased_table_LIGHT_scored.csv.gz\n")

cat("\n🎯 DONE\n")
cat("======\n")