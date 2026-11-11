#### 📦 Load Libraries ####
library(data.table)
library(ggplot2)
library(ggupset)
library(patchwork)
library(imputeLCMD)
library(limma)
library(ggrepel)

# Make script execution reproducible
set.seed(123)

# Set a plotting theme for the whole document
theme_set(theme(
  text = element_text(size = 8),
  panel.border = element_rect(linewidth = 0.25, colour = "black", fill = NA),
  panel.background = element_blank(),
  panel.grid = element_blank(),
  axis.text = element_text(colour = "black"),
  strip.background = element_blank(),
  strip.text = element_text(size = 8, vjust = 1)
))

#### 🧬 Load and Prepare Data ####
DT <- fread("protein_quant_new.csv.gz")

# Keep only proteins with LFQ intensity values
DT <- DT[!is.na(LFQ)]

# Log2-transform LFQ intensities
DT[, log2_quant := log2(LFQ)]

# Print selected columns
DT[1:5, .(Sample, Protein_group, log2_quant)]

#### 📊 PCA Analysis ####
# Run a PCA to check for major trends
PCA_input <- dcast(DT, Protein_group ~ Sample, value.var = "log2_quant")
PCA_input <- PCA_input[complete.cases(PCA_input)]
PCA_input <- t(data.frame(PCA_input, row.names = "Protein_group"))
DT_PCA <- prcomp(PCA_input, scale. = TRUE)
var_explained <- summary(DT_PCA)$importance
DT_PCA <- data.table(DT_PCA$x, keep.rownames = "Sample")
DT_PCA <- merge(DT_PCA, unique(DT[, .(Sample, Condition)]))

# Create PCA plot
pPCA <- ggplot(DT_PCA, aes(x = PC1, y = PC2, colour = Condition)) +
  geom_point(size = 2) +
  scale_colour_brewer(palette = "Set1") +
  xlab(paste("PC1 [", round(var_explained[2, 1] * 100, 0), "%]", sep = "")) +
  ylab(paste("PC2 [", round(var_explained[2, 2] * 100, 0), "%]", sep = "")) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.title = element_blank())

# Create barplot with number of remaining proteins per sample
pBar <- ggplot(DT[, .N, .(Sample, Condition)], aes(x = Condition, y = N, fill = Sample, label = N)) +
  geom_bar(stat = "identity", position = "dodge2") +
  geom_text(position = position_dodge2(width = 0.9), vjust = 0.5, hjust = 1.2, size = 3, angle = 90) +
  scale_fill_brewer() +
  labs(x = NULL, y = "Protein count", fill = NULL)

# Arrange the two plots into a single figure
pPCA + pBar + plot_annotation(tag_levels = "A")

#### 🔍 Protein Filtering ####
# Calculate protein stats before filtering
protein_stats <- DT[, .(
  n_samples = .N,
  n_conditions = uniqueN(Condition)
), Protein_group]

# Keep proteins detected in at least 2 samples per condition
keepers <- protein_stats[n_samples >= 4, Protein_group]

# Filter table
DT <- DT[Protein_group %in% keepers]

cat("Proteins retained after filtering:", uniqueN(DT$Protein_group), "\n")

#### 📐 UpSet Plot ####
pUpSet <- ggplot(
  DT[, .(Sample = list(Sample)), by = Protein_group],
  aes(x = Sample)
) +
  geom_bar() +
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5, size = 3) +
  scale_x_upset(order_by = "degree", reverse = TRUE,
                sets = DT[order(Condition), unique(Sample)]) +
  labs(x = NULL, y = "Protein count") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

pUpSet

#### 📏 Median Normalization ####
# Identify the "median of the medians"
target_median <- DT[, median(log2_quant), Sample][, median(V1)]

# Scale protein quantities to obtain this target_median for all samples
DT[, log2_quant_norm := log2_quant - median(log2_quant) + target_median, by = Sample]

# Plot before and after normalization
pNorm <- ggplot(
  data = melt(DT, measure.vars = c("log2_quant", "log2_quant_norm")),
  aes(x = Condition, fill = Sample, y = value)
) +
  facet_wrap("variable", labeller = labeller(variable = c("log2_quant" = "before", "log2_quant_norm" = "after"))) +
  geom_boxplot(outliers = FALSE) +
  scale_fill_brewer() +
  labs(x = NULL, y = "log2 protein quantities")

pNorm

#### 🧮 LIMMA Analysis - FIXED ####
# Extract sample to condition mapping
sample_map <- unique(DT[, .(Sample, Condition)])

# Create design matrix
design <- model.matrix(~ 0 + Condition, data = sample_map)
colnames(design) <- gsub("Condition", "", colnames(design))
rownames(design) <- sample_map$Sample

# Cast data into wide format
wDT <- dcast(DT, Protein_group ~ Sample, value.var = "log2_quant_norm", fill = NA)

# Set column order in wDT to match design matrix
setcolorder(wDT, neworder = c("Protein_group", sample_map$Sample))

# Check missing values
cat("Missing values:", sum(is.na(wDT[, -"Protein_group"])), "\n")

# Convert to matrix with Protein_group as rownames (CRITICAL FIX)
wDT_matrix <- as.matrix(wDT[, -"Protein_group"])
rownames(wDT_matrix) <- wDT$Protein_group

# Fit the linear model with the matrix
fit <- lmFit(wDT_matrix, design)

# Define contrasts (DMSO as baseline)
cont.matrix <- makeContrasts(
  SJF8240_vs_DMSO = SJF8240 - DMSO,
  Foretinib_vs_DMSO = Foretinib - DMSO,
  levels = design
)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

#### 📋 Extract Results - FIXED ####
extract_res <- function(name) {
  # Get the topTable result with row names preserved
  tt <- topTable(fit2, coef = name, number = Inf, adjust = "BH")
  
  # Convert to data.table with row names as Protein_group
  tbl <- as.data.table(tt, keep.rownames = "Protein_group")
  tbl[, Contrast := name]
  
  return(tbl)
}

results <- rbindlist(lapply(colnames(cont.matrix), extract_res))

# Add gene and protein names
results <- merge(
  results,
  unique(DT[, .(Protein_group, Genes, Protein_description)]),
  by = "Protein_group",
  all.x = TRUE
)

# Check the merge worked
cat("\n=== Checking merge ===\n")
cat("Total results:", nrow(results), "\n")
cat("Results with gene names:", sum(!is.na(results$Genes)), "\n")
cat("First few rows with genes:\n")
print(results[!is.na(Genes)][1:5, .(Protein_group, Genes, logFC, P.Value)])

# Set column order
setcolorder(results, c("Protein_group", "Genes", "Protein_description", "Contrast", 
                       "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"))

# Save results
fwrite(results, "LIMMA_results_DMSO_vs_SJF8240_Foretinib_new.csv.gz")

# Print summary
cat("\n=== Results Summary ===\n")
results[, .(
  Total = .N,
  Sig_0.05 = sum(adj.P.Val < 0.05, na.rm = TRUE),
  Sig_0.1 = sum(adj.P.Val < 0.1, na.rm = TRUE)
), by = Contrast]

#### 🌋 Volcano Plots with Gene Labels ####
# Cut-offs
fc_cutoff <- 1
p_cutoff  <- 0.05

# Significance + highlights
results[, significant := abs(logFC) > fc_cutoff & P.Value < p_cutoff]
results[, highlight   := Genes %in% c("MET", "MAPK14")]

# Safety: ensure Genes is character
results[, Genes := as.character(Genes)]

# Check how many proteins are significant
cat("\n=== Significant Proteins Check ===\n")
results[, .(
  Total = .N,
  Significant = sum(significant, na.rm = TRUE),
  With_Gene_Names = sum(!is.na(Genes)),
  Sig_with_Genes = sum(significant & !is.na(Genes), na.rm = TRUE)
), by = Contrast]

# Show top significant proteins
cat("\n=== Top Significant Proteins ===\n")
print(results[significant == TRUE][order(P.Value)][1:10, .(Contrast, Genes, logFC, P.Value, adj.P.Val)])

# Check MET and MAPK14
cat("\n=== MET and MAPK14 ===\n")
print(results[Genes %in% c("MET", "MAPK14"), .(Protein_group, Genes, Contrast, logFC, P.Value)])

# Function to create volcano plot
make_volcano <- function(contrast_name, title_label) {
  df <- results[Contrast == contrast_name]
  
  cat("\n=== ", title_label, " ===\n")
  cat("Total proteins:", nrow(df), "\n")
  cat("Significant proteins:", sum(df$significant, na.rm = TRUE), "\n")
  cat("Sig with gene names:", sum(df$significant & !is.na(df$Genes), na.rm = TRUE), "\n")
  
  ggplot(df, aes(x = logFC, y = -log10(P.Value), colour = significant)) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dotted", colour = "grey50", linewidth = 0.25) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dotted", colour = "grey50", linewidth = 0.25) +
    geom_point(alpha = 0.7, size = 1.5) +
    # Label significant proteins (not MET/MAPK14) - inherits red color
    geom_text_repel(
      data = df[significant == TRUE & highlight == FALSE & !is.na(Genes)],
      aes(label = Genes),
      size = 2.5,
      fontface = "bold",
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.size = 0.3,
      min.segment.length = 0,
      seed = 42
    ) +
    # Highlight MET and MAPK14 in blue
    geom_point(data = df[highlight == TRUE], colour = "deepskyblue2", size = 3) +
    geom_label(
      data = df[highlight == TRUE & !is.na(Genes)],
      aes(label = Genes),
      colour = "deepskyblue2",
      size = 3,
      fontface = "bold"
    ) +
    scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "red3")) +
    labs(x = paste(title_label, "\nlog2 fold-change"),
         y = "-log10 P-value") +
    theme(legend.position = "none")
}

# Generate & save plots
p1 <- make_volcano("SJF8240_vs_DMSO", "SJF8240 vs DMSO")
p2 <- make_volcano("Foretinib_vs_DMSO", "Foretinib vs DMSO")

ggsave("Volcano_SJF8240_vs_DMSO.png", p1, width = 7, height = 7, dpi = 300)
ggsave("Volcano_Foretinib_vs_DMSO.png", p2, width = 7, height = 7, dpi = 300)

combined <- p1 + p2 + plot_annotation(tag_levels = "A")
ggsave("Volcano_Combined.png", combined, width = 14, height = 7, dpi = 300)

combined