############################################################################################################# #
# Transform SILAC-labelled precursor-level output from DIA-NN into protein quantities
# Date: Updated Feb 2025
############################################################################################################# #

# Load required libraries
library(data.table); library(iq); library(ggplot2); library(patchwork); library(scales)

#### Load the report and parameterise the script ####

# Load the precursor report (relevant columns only)
report <- fread("report.tsv",  # Default is report.tsv.gz (second-pass/MBR). Alternatively use report-first-pass.tsv.gz
                select = c("Run", "Protein.Group", "Protein.Ids", "First.Protein.Description", "Genes", "Stripped.Sequence", "Precursor.Id", "Proteotypic", 
                           "Precursor.Quantity",    # Note: we use unnormalised quantities for SILAC data
                           "Precursor.Translated", "Channel.Q.Value",  # These are SILAC-specific columns
                           "Q.Value", "Global.Q.Value", "PG.Q.Value", "Global.PG.Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value", "Precursor.Charge", "RT", "Mass.Evidence", "RT.Start", "RT.Stop"))

# Is this a first- or second-pass report? Required for FDR filtering.
workflow_type <- "second"  # Either "second" (default) or "first"

# Show full sample names
report[, .N, Run]


library(data.table)

# ============================
# Replace dashes with underscores
# ============================
report[, Run := gsub("-", "_", Run)]

# ============================
# Remove trailing timestamp if present (e.g. _20260107005249)
# ============================
report[, Run_clean := sub("_[0-9]{14}$", "", Run)]

# ============================
# Extract Condition: DMSO, Foretinib, SJF8240
#   - map PBS -> DMSO
# ============================
report[, Condition := fifelse(grepl("_Foretinib_", Run_clean, ignore.case = TRUE), "Foretinib",
                              fifelse(grepl("_SJF8240_",   Run_clean, ignore.case = TRUE), "SJF8240",
                                      fifelse(grepl("_PBS_",      Run_clean, ignore.case = TRUE), "DMSO",
                                              fifelse(grepl("_DMSO_",     Run_clean, ignore.case = TRUE), "DMSO",
                                                      NA_character_))))]

# ============================
# Extract Replicate from end of Run_clean: ..._1 / ..._2 / ..._3
# ============================
report[, Replicate := fifelse(grepl("_(\\d+)$", Run_clean),
                              paste0("R", sub(".*_(\\d+)$", "\\1", Run_clean)),
                              NA_character_)]

# ============================
# Construct Sample name: Condition_R#
# ============================
report[, Sample := fifelse(!is.na(Condition) & !is.na(Replicate),
                           paste0(Condition, "_", Replicate),
                           NA_character_)]

# ============================
# Order samples by Condition then Replicate
# ============================
desired_condition_order <- c("DMSO", "Foretinib", "SJF8240")
desired_replicate_order <- c("R1", "R2", "R3")

report[, Condition := factor(Condition, levels = desired_condition_order)]
report[, Replicate := factor(Replicate, levels = desired_replicate_order)]

report <- report[order(Condition, Replicate)]
report[, Sample := factor(Sample, levels = unique(Sample))]

# ============================
# Preview
# ============================
head(report[, .(Run, Run_clean, Condition, Replicate, Sample)], 10)
tail(report[, .(Run, Run_clean, Condition, Replicate, Sample)], 10)


# NB: Code below here can usually be executed without further adjustment

#### Data processing section ####

# ----------------------------
# Helper: safe filter only if column exists
# ----------------------------
apply_filters <- function(dt, filters) {
  for (f in filters) {
    col <- f$col
    if (col %in% names(dt)) {
      dt <- dt[ get(col) <= f$thr ]
    }
  }
  dt
}

# ----------------------------
# Apply FDR filters dynamically based on workflow type
# ----------------------------
if (workflow_type == "second") {
  filters <- list(
    list(col = "Q.Value",          thr = 0.01),
    list(col = "PG.Q.Value",       thr = 0.05),
    list(col = "Lib.Q.Value",      thr = 0.01),
    list(col = "Lib.PG.Q.Value",   thr = 0.01),
    list(col = "Channel.Q.Value",  thr = 0.01)
  )
  report <- apply_filters(report, filters)
  message("Second-pass FDR filtering applied")
  
} else if (workflow_type == "first") {
  filters <- list(
    list(col = "Q.Value",          thr = 0.01),
    list(col = "PG.Q.Value",       thr = 0.05),
    list(col = "Global.Q.Value",   thr = 0.01),
    list(col = "Global.PG.Q.Value",thr = 0.01),
    list(col = "Channel.Q.Value",  thr = 0.01)
  )
  report <- apply_filters(report, filters)
  message("Single-pass FDR filtering applied")
} else {
  stop("workflow_type must be 'first' or 'second'")
}

# ----------------------------
# Convert zeroes to NAs for precursor quantities
# ----------------------------
if ("Precursor.Quantity" %in% names(report)) {
  report[ Precursor.Quantity == 0, Precursor.Quantity := NA_real_ ]
} else stop("Missing column: Precursor.Quantity")

if ("Precursor.Translated" %in% names(report)) {
  report[ Precursor.Translated == 0, Precursor.Translated := NA_real_ ]
} else {
  # If you don't have translated quantities, create as NA (or set equal to Quantity if that's intended)
  report[, Precursor.Translated := NA_real_ ]
}

# Remove precursors that were not quantified at all (NAs are ok for translated values)
report <- report[ !is.na(Precursor.Quantity) ]

# ----------------------------
# Normalise Precursor.Quantity to adjust for unequal sample loading
# ----------------------------
target_sum <- report[, .(sumQ = sum(Precursor.Quantity, na.rm = TRUE)), by = Run][, median(sumQ)]

report[, Precursor.Quantity :=
         Precursor.Quantity / sum(Precursor.Quantity, na.rm = TRUE) * target_sum,
       by = Run]

# FIXED: by = Run (your version had a stray ", Run" which is wrong)
report[, Precursor.Translated :=
         Precursor.Translated / sum(Precursor.Translated, na.rm = TRUE) * target_sum,
       by = Run]

# ----------------------------
# Process SILAC labels (based on Precursor.Id containing SILAC labels)
# ----------------------------
if (!("Precursor.Id" %in% names(report))) stop("Missing column: Precursor.Id")

# Control: should be zero (no precursors with both labels)
report[ Precursor.Id %like% "SILAC-.-L" & Precursor.Id %like% "SILAC-.-H", .N ]

report[ Precursor.Id %like% "SILAC-.-L", label := "L" ]
report[ Precursor.Id %like% "SILAC-.-H", label := "H" ]
report[, Precursor.Id.nolabels := gsub("SILAC-.-.", "SILAC", Precursor.Id)]

# ----------------------------
# Intra-sample comparison: H/L ratios using translated quantities
# ----------------------------
SILAC <- dcast(
  report,
  formula   = Precursor.Id.nolabels + Protein.Group + Run ~ label,
  value.var = "Precursor.Translated"
)

# Avoid log2(H/L) when missing or zero
SILAC <- SILAC[ !is.na(H) & !is.na(L) & H > 0 & L > 0 ]
SILAC[, precursor_log2_HL_ratio := log2(H / L)]

SILAC <- SILAC[
  ,
  .(log2_HL_ratio = median(precursor_log2_HL_ratio, na.rm = TRUE),
    N_ratios      = .N),
  by = .(Run, Protein.Group)
]

# ----------------------------
# Cross-sample comparisons: MaxLFQ on non-translated quantities
# ----------------------------
LFQ_fun <- function(dt, LFQ_colname) {
  tmp <- dt[, .(
    protein_list = Protein.Group,
    sample_list  = Run,
    id           = Precursor.Id,
    quant        = log2(Precursor.Quantity)
  )]
  
  tmp <- fast_MaxLFQ(tmp)
  
  tmp <- data.table(
    tmp$estimate,
    Precursor_group = tmp$annotation,
    keep.rownames = "Protein.Group"
  )
  
  tmp <- melt(
    tmp,
    id.vars = c("Protein.Group", "Precursor_group"),
    variable.name = "Run",
    value.name = LFQ_colname
  )
  
  tmp <- tmp[ Precursor_group == "" & !is.na(get(LFQ_colname)) ][, -"Precursor_group"]
  
  tmp[, (LFQ_colname) := 2^get(LFQ_colname)]
  
  return(tmp)  # FIXED: your function didn’t return anything
}

LFQ_T <- LFQ_fun(report, "LFQ")
LFQ_H <- LFQ_fun(report[label == "H"], "LFQ_H")
LFQ_L <- LFQ_fun(report[label == "L"], "LFQ_L")

# ----------------------------
# Crude "total signal": sum precursor quantities per protein+run
# ----------------------------
int_T <- report[,             .(Intensity   = sum(Precursor.Quantity, na.rm = TRUE)), by = .(Protein.Group, Run)]
int_H <- report[label == "H", .(Intensity_H = sum(Precursor.Quantity, na.rm = TRUE)), by = .(Protein.Group, Run)]
int_L <- report[label == "L", .(Intensity_L = sum(Precursor.Quantity, na.rm = TRUE)), by = .(Protein.Group, Run)]

# ----------------------------
# Precursor counts per protein+run
# ----------------------------
if (!("Proteotypic" %in% names(report))) report[, Proteotypic := 0L]

counts_H <- report[label == "H",
                   .(N_precursors_H = uniqueN(Precursor.Id),
                     N_precursors_proteotypic_H = sum(Proteotypic, na.rm = TRUE)),
                   by = .(Protein.Group, Run)
]

counts_L <- report[label == "L",
                   .(N_precursors_L = uniqueN(Precursor.Id),
                     N_precursors_proteotypic_L = sum(Proteotypic, na.rm = TRUE)),
                   by = .(Protein.Group, Run)
]

# ----------------------------
# Protein annotations (NO Timepoint in your new design)
# ----------------------------
ann_cols <- c("Protein.Group", "Run", "Sample", "Condition", "Replicate",
              "First.Protein.Description", "Genes")
ann_cols <- ann_cols[ann_cols %in% names(report)]

annotations <- unique(report[, ..ann_cols])

# ----------------------------
# Merge everything
# ----------------------------
merge_function <- function(x, y) merge(x, y, by = c("Protein.Group", "Run"), all = TRUE)

proteins <- Reduce(merge_function, list(
  LFQ_T, LFQ_H, LFQ_L,
  int_T, int_H, int_L,
  counts_H, counts_L,
  SILAC,
  annotations
))

# Remove protein-sample pairs without values
proteins <- proteins[ !is.na(Intensity) ]

# Tidy names
if ("First.Protein.Description" %in% names(proteins)) {
  setnames(proteins,
           old = c("Protein.Group", "First.Protein.Description"),
           new = c("Protein_group", "Protein_description"))
} else {
  setnames(proteins, old = "Protein.Group", new = "Protein_group")
}

# Order columns (by name; robust)
preferred_first <- intersect(
  c("Run", "Sample", "Condition", "Replicate", "Protein_group", "Protein_description", "Genes",
    "LFQ", "LFQ_H", "LFQ_L",
    "Intensity", "Intensity_H", "Intensity_L",
    "log2_HL_ratio", "N_ratios",
    "N_precursors_H", "N_precursors_L",
    "N_precursors_proteotypic_H", "N_precursors_proteotypic_L"),
  names(proteins)
)
setcolorder(proteins, c(preferred_first, setdiff(names(proteins), preferred_first)))

# Apply sample order from report (if Sample exists)
if ("Sample" %in% names(proteins) && "Sample" %in% names(report)) {
  proteins[, Sample := factor(Sample, levels = levels(report$Sample))]
}

# Write out protein result table (optional)
#fwrite(proteins, "protein_quant.csv.gz")
fwrite(report, "reportfile.tsv")




#### QC plotting section ####

# ----------------------------
# Set Sample plotting order: Condition then Replicate
# ----------------------------
desired_condition_order <- c("DMSO", "Foretinib", "SJF8240")
desired_replicate_order <- c("R1", "R2", "R3")

# Ensure factors exist (if not already set upstream)
report[,   Condition := factor(Condition, levels = desired_condition_order)]
report[,   Replicate := factor(Replicate, levels = desired_replicate_order)]
proteins[, Condition := factor(Condition, levels = desired_condition_order)]
proteins[, Replicate := factor(Replicate, levels = desired_replicate_order)]

sample_levels <- proteins[order(Condition, Replicate), unique(Sample)]
report[,   Sample := factor(Sample, levels = sample_levels)]
proteins[, Sample := factor(Sample, levels = sample_levels)]

# ----------------------------
# Theme (keep yours)
# ----------------------------
theme_set(
  theme(
    line = element_line(linewidth = 0.25),
    text = element_text(size = 7),
    rect = element_rect(colour = "black", linewidth = 0.25, fill = NA),
    panel.border = element_rect(linewidth = 0.25),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(colour = "black", size = 6),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key.width = unit(4, "mm"),
    legend.key = element_rect(colour = NA),
    strip.background = element_blank(),
    strip.text = element_text(vjust = 1)
  )
)

# ----------------------------
# Colours
# ----------------------------
HL_cols <- c(heavy = "#37abc8ff", H = "#37abc8ff", light = "#ca00a3ff", L = "#ca00a3ff", total = "black")
scale_HL_labels <- c(L = "light", H = "heavy")

# Stable condition colours (avoid re-order surprises)
cond_levels <- levels(proteins$Condition)
condition_color_mapping <- hue_pal()(length(cond_levels))
names(condition_color_mapping) <- cond_levels

# ----------------------------
# (1) Intensity distributions
# ----------------------------
plot_int_DT <- rbindlist(list(
  report[!is.na(Precursor.Quantity),
         .(Sample, `log2 quantity` = log2(pmax(Precursor.Quantity, 1)), Type = "Precursor quantity")],
  
  proteins[!is.na(LFQ),
           .(Sample, `log2 quantity` = log2(pmax(LFQ, 1)), Type = "Protein MaxLFQ")],
  
  proteins[!is.na(Intensity),
           .(Sample, `log2 quantity` = log2(pmax(Intensity, 1)), Type = "Protein Intensity")]
), use.names = TRUE, fill = TRUE)

plot_int_DT[, Type := factor(Type, levels = c("Precursor quantity", "Protein MaxLFQ", "Protein Intensity"))]

pInt <- ggplot(plot_int_DT, aes(x = Sample, y = `log2 quantity`)) +
  facet_wrap(~Type, scales = "free_x") +
  geom_boxplot(outliers = FALSE, linewidth = 0.25) +
  theme(axis.title.y = element_blank()) +
  coord_flip()

plot(pInt)

# ----------------------------
# (2) Cumulative intensity distributions (Light vs Heavy)
# ----------------------------
cum_cols <- intersect(c("LFQ_H", "LFQ_L", "Intensity_H", "Intensity_L"), names(proteins))

plot_int <- proteins[, lapply(.SD, function(x) sum(x, na.rm = TRUE)), by = Sample, .SDcols = cum_cols]
plot_int <- melt(plot_int, id.vars = "Sample", value.name = "Cumulative quantity", variable.name = "label")

plot_int[, Type := fifelse(label %like% "LFQ", "Protein MaxLFQ", "Protein intensity")]
plot_int[, label := gsub(".+_", "", label)]
plot_int[, Type := factor(Type, levels = c("Protein MaxLFQ", "Protein intensity"))]

pCumInt <- ggplot(plot_int, aes(x = Sample, y = `Cumulative quantity`, fill = label)) +
  facet_wrap(~Type, scales = "free_x") +
  geom_bar(stat = "identity") +
  ylab("Cumulative protein intensity") +
  scale_fill_manual(values = HL_cols, guide = guide_legend(title = "SILAC\nlabel"), labels = scale_HL_labels) +
  theme(axis.title.y = element_blank()) +
  coord_flip()

plot(pCumInt)

# ----------------------------
# (3) Counts: precursors, peptides, protein groups
# ----------------------------
count_cols <- intersect(c("Precursor.Id", "Stripped.Sequence", "Protein.Group"), names(report))

# If label missing (non-SILAC), make a dummy label so plot still runs
if (!("label" %in% names(report))) report[, label := "total"]

counts <- report[, lapply(.SD, uniqueN), by = .(Sample, label), .SDcols = count_cols]
counts <- melt(counts, id.vars = c("Sample", "label"), value.name = "IDs")

facet_labels <- as_labeller(c(
  Precursor.Id = "Precursors",
  Stripped.Sequence = "Peptides",
  Protein.Group = "Protein groups"
))

pCount <- ggplot(counts, aes(x = Sample, y = IDs/1000, fill = label,
                             label = format(IDs, big.mark = ",", scientific = FALSE))) +
  facet_wrap(~variable, scales = "free_x", labeller = facet_labels) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(size = 2, hjust = 1.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = HL_cols, guide = guide_legend(title = "SILAC\nlabel"), labels = scale_HL_labels) +
  ylab("Number of IDs [x1,000]") +
  theme(axis.title.y = element_blank()) +
  coord_flip()

plot(pCount)

# ----------------------------
# (4) Data completeness vs precursor cut-offs
# ----------------------------
f_prot_count <- function(dt, var_name) {
  tmp <- dt[, .(N_samples = .N), by = Protein_group][, .(N_proteins = .N), by = N_samples]
  tmp <- tmp[order(-N_samples)]
  tmp[, cumulative_protein_N := cumsum(N_proteins)]
  tmp[, Precursors := var_name]
  tmp
}

data_completeness <- rbindlist(list(
  f_prot_count(proteins, "all"),
  f_prot_count(proteins[!is.na(N_precursors_L) & N_precursors_L >= 2], "min 2 light"),
  f_prot_count(proteins[!is.na(N_precursors_H) & N_precursors_H >= 2], "min 2 heavy"),
  f_prot_count(proteins[!is.na(N_precursors_L) & !is.na(N_precursors_H) & N_precursors_L >= 2 & N_precursors_H >= 2],
               "min 2 light &\nmin 2 heavy")
), use.names = TRUE, fill = TRUE)

pNAs <- ggplot(data_completeness, aes(x = N_samples, y = cumulative_protein_N/1000, colour = Precursors)) +
  geom_point(size = 0.8) +
  geom_line(linewidth = 0.25) +
  xlab("Number of samples") +
  ylab("Number of proteins [x1,000]") +
  scale_x_continuous(breaks = seq(1, max(data_completeness$N_samples), 1)) +
  scale_y_continuous(limits = c(0, max(data_completeness$cumulative_protein_N, na.rm = TRUE)/1000)) +
  scale_color_brewer(palette = "Dark2") +
  theme(panel.grid.major = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
        panel.grid.minor = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"))

plot(pNAs)

# ----------------------------
# (5a) Precursor CVs across replicates within each Condition + label
#     IMPORTANT: compute per precursor across Samples (replicates), not pooled rows
# ----------------------------
CV_function <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x) / mean(x) * 100
}

# Build per-sample precursor totals first (median across rows per sample)
precursor_by_sample <- report[
  !is.na(Precursor.Quantity) & !is.na(label) & !is.na(Condition),
  .(qty = median(Precursor.Quantity, na.rm = TRUE)),
  by = .(Precursor.Id.nolabels, label, Condition, Sample)
]

precursor_CVs <- precursor_by_sample[
  ,
  .(PrecursorCV = CV_function(qty), N_samples = .N),
  by = .(Precursor.Id.nolabels, label, Condition)
][N_samples >= 3 & is.finite(PrecursorCV)]

precursor_CVs[, precursor_rank := frank(PrecursorCV), by = .(Condition, label)]
precursor_CVs[, label := factor(label, levels = c("L", "H"))]

pCV1 <- ggplot(precursor_CVs, aes(x = precursor_rank/1000, y = PrecursorCV, colour = Condition)) +
  facet_wrap(~label, labeller = as_labeller(c(L = "Light", H = "Heavy"))) +
  geom_line(linewidth = 0.25) +
  labs(x = "Number of precursors [x1,000]", y = "precursor CV [%]") +
  coord_cartesian(ylim = c(0, 50)) +
  scale_color_manual(values = condition_color_mapping) +
  theme(panel.grid.major = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
        panel.grid.minor = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
        legend.position = "none")

plot(pCV1)

# ----------------------------
# (5b) Protein CVs across replicates within each Condition + label
# ----------------------------
protein_long <- melt(
  proteins,
  id.vars = c("Protein_group", "Condition", "Sample"),
  measure.vars = intersect(c("LFQ_L", "LFQ_H"), names(proteins)),
  variable.name = "label",
  value.name = "LFQ_intensity",
  na.rm = TRUE
)
protein_long[, label := gsub("LFQ_", "", label)]

protein_CVs <- protein_long[
  ,
  .(ProteinCV = CV_function(LFQ_intensity), N_samples = uniqueN(Sample)),
  by = .(Protein_group, label, Condition)
][N_samples >= 3 & is.finite(ProteinCV)]

protein_CVs[, protein_rank := frank(ProteinCV), by = .(Condition, label)]
protein_CVs[, label := factor(label, levels = c("L", "H"))]

pCV2 <- ggplot(protein_CVs, aes(x = protein_rank/1000, y = ProteinCV, colour = Condition)) +
  facet_wrap(~label, labeller = as_labeller(c(L = "Light", H = "Heavy"))) +
  geom_line(linewidth = 0.25) +
  labs(x = "Number of proteins [x1,000]", y = "protein CV [%]") +
  coord_cartesian(ylim = c(0, 50)) +
  scale_color_manual(values = condition_color_mapping) +
  theme(panel.grid.major = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
        panel.grid.minor = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
        legend.position = "none")

plot(pCV2)

# ----------------------------
# PCA (Light and Heavy) - NO Timepoint, use Replicate as shape
# ----------------------------
PCA_input_L <- dcast(proteins, Protein_group ~ Sample, value.var = "LFQ_L")
PCA_input_L <- PCA_input_L[complete.cases(PCA_input_L)]
PCA_input_L <- t(data.frame(PCA_input_L, row.names = "Protein_group"))
DT_PCA_L <- prcomp(PCA_input_L, scale. = TRUE)
var_explained_L <- summary(DT_PCA_L)$importance
DT_PCA_L <- data.table(DT_PCA_L$x, keep.rownames = "Sample")
DT_PCA_L <- merge(DT_PCA_L, unique(proteins[, .(Sample, Condition, Replicate)]), by = "Sample", all.x = TRUE)

PCA_input_H <- dcast(proteins, Protein_group ~ Sample, value.var = "LFQ_H")
PCA_input_H <- PCA_input_H[complete.cases(PCA_input_H)]
PCA_input_H <- t(data.frame(PCA_input_H, row.names = "Protein_group"))
DT_PCA_H <- prcomp(PCA_input_H, scale. = TRUE)
var_explained_H <- summary(DT_PCA_H)$importance
DT_PCA_H <- data.table(DT_PCA_H$x, keep.rownames = "Sample")
DT_PCA_H <- merge(DT_PCA_H, unique(proteins[, .(Sample, Condition, Replicate)]), by = "Sample", all.x = TRUE)

pPCA_L <- ggplot(DT_PCA_L, aes(x = PC1, y = PC2, colour = Condition, shape = Replicate)) +
  ggtitle("Light proteins") +
  geom_point(size = 2.5) +
  xlab(paste0("PC1 [", round(var_explained_L["Proportion of Variance", "PC1"] * 100, 0), "%]")) +
  ylab(paste0("PC2 [", round(var_explained_L["Proportion of Variance", "PC2"] * 100, 0), "%]")) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_blank(), legend.position = "right")
plot(pPCA_L)

pPCA_H <- ggplot(DT_PCA_H, aes(x = PC1, y = PC2, colour = Condition, shape = Replicate)) +
  ggtitle("Heavy proteins") +
  geom_point(size = 2.5) +
  xlab(paste0("PC1 [", round(var_explained_H["Proportion of Variance", "PC1"] * 100, 0), "%]")) +
  ylab(paste0("PC2 [", round(var_explained_H["Proportion of Variance", "PC2"] * 100, 0), "%]")) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_blank(), legend.position = "bottom")
plot(pPCA_H)

# ----------------------------
# Combine plots (patchwork)
# ----------------------------
combined_plot <- (
  pInt /
    pCumInt /
    pCount /
    patchwork::wrap_elements(pNAs + plot_spacer()) /
    patchwork::wrap_elements(pCV1 + pCV2) /
    patchwork::wrap_elements(pPCA_L + pPCA_H)
) +
  plot_layout(ncol = 1, heights = c(1, 1, 1.8, 1, 1, 1)) +
  plot_annotation(theme = theme(plot.margin = margin(20, 20, 20, 20)))

plot(combined_plot)

ggsave(
  "DIANNprotein_quant_plot_LARGE.png",
  combined_plot,
  width = 35, height = 80, units = "cm",
  dpi = 300, bg = "white"
)
