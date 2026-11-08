############################################################################################################# #
# Transform label-free, precursor-level output from DIA-NN into protein quantities
# Date: Updated Feb 2025
############################################################################################################# #

# Load required libraries
library(data.table); library(iq); library(ggplot2); library(patchwork)

#### Load the report and parameterise the script ####

# Load the precursor report (relevant columns only)
report <- fread("report (copy 1).tsv",  # Default is report.tsv.gz (second-pass/MBR). Alternatively use report-first-pass.tsv.gz
                select = c("Run", "Protein.Group", "Protein.Ids", "First.Protein.Description", "Genes", "Stripped.Sequence", "Precursor.Id", "Proteotypic", 
                           "Precursor.Normalised",   # Note: we use DIA-NN's normalised data for LF data
                           "Q.Value", "Global.Q.Value", "PG.Q.Value", "Global.PG.Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value"))

# Is this a first- or second-pass report? Required for FDR filtering.
workflow_type <- "second"  # Either "second" (default) or "first"

# Show full sample names
report[, .N, Run]

# # Run provides the full sample names. Use regex patterns to extract three annotation columns:
#   # Sample: Simplified but still unique(!) sample names, e.g. R1_SS48, R1_SS6, R2_SS48
#   # Condition: Name of experimental condition/ timepoint, e.g. SS6, SS48
#   # Replicates: Name or replicate, e.g. R1, R2, R3
# report[,    Sample := gsub(".*_(R\\d_SS\\d+)_.*", "\\1", Run)] 
# report[, Condition := gsub("R\\d_", "", Sample)] 
# report[, Replicate := gsub("_.+", "", Sample)] 
# 
# # Check if it worked as expected. You are aiming for "sample = condition + replicate = unique, simplified version of Run" 
# report[, .N, .(Run, Sample, Condition, Replicate)]

# Extract replicate, condition, and simplified sample from Run (minimal)
report[, Replicate := sub(".*_(\\d+)$", "R\\1", Run)]
report[, Condition := sub(".*_([^_]+)_\\d+$", "\\1", Run)]
report[,    Sample := paste(Replicate, Condition, sep = "_")]

# Check
report[, .N, .(Run, Sample, Condition, Replicate)]



# Enable this to use only proteotypic peptides (disabled by default, but number of proteotypic peptides per protein is always reported)
# report <- report[ Proteotypic == 1 ]

# NB: Code below here can usually be executed without further adjustment, 
# except for the dimensions of the plot (final line) in some cases.


#### Data processing section ####

# Apply FDR filters dynamically based on workflow type
if (workflow_type == "second") {
  report <- report[ Q.Value <= 0.01           &   # Set run-specific precursor FDR
                      PG.Q.Value <= 0.05        &   # Set run-specific protein group FDR (recommended at 5% by DIA-NN Github)
                      Lib.Q.Value <= 0.01       &   # Set global precursor FDR for library entry, only use for MBR-style two-pass reports
                      Lib.PG.Q.Value <= 0.01    ]   # Set global protein group FDR for library entry, only use for MBR-style two-pass reports
  print("Second-pass FDR filtering applied")
  
} else if (workflow_type == "first") {
  report <- report[ Q.Value <= 0.01           &   # Set run-specific precursor FDR
                      PG.Q.Value <= 0.05        &   # Set run-specific protein group FDR (recommended at 5% by DIA-NN Github)
                      Global.Q.Value <= 0.01    &   # Set global precursor FDR (use instead of 'Lib.Q.Value' for single-pass searches)
                      Global.PG.Q.Value <= 0.01 ]   # Set global protein group FDR (use this instead of 'Lib.PG.Q.Value' for single-pass searches)
  print("Single-pass FDR filtering applied") 
}

# Remove precursors that were not quantified
report <- report[ !is.na(Precursor.Normalised) & Precursor.Normalised != 0 ]

# Define function to calculate LFQs using iq's fast_MaxLFQ. It returns a protein-sample matrix in long format
LFQ_fun <- function(dt, LFQ_colname){                     
  tmp <- dt[, .(protein_list = Protein.Group, sample_list = Run, id = Precursor.Id, quant = log2( Precursor.Normalised ))] # Transform input table as required by iq package
  tmp <- fast_MaxLFQ( tmp )                                                                                          # Apply the fast_MaxLFQ function
  tmp <- data.table( tmp$estimate, Precursor_group = tmp$annotation, keep.rownames = "Protein.Group" )               # Transform back into data.table with LFQ value and precursor grouping for each protein...
  tmp <- melt(tmp, id.vars = c("Protein.Group", "Precursor_group"), variable.name = "Run", value.name = LFQ_colname) # ... and melt into long-format with appropriate colname
  tmp <- tmp[ Precursor_group == "" & !is.na( get(LFQ_colname) ) ][, -"Precursor_group" ]                            # Return output without proteins that could not be quantified with connected precursors, remove precursor grouping column
  tmp <- tmp[, (LFQ_colname) := 2^get(LFQ_colname)]                                                                  # Un-log the data
}

# Calculate MaxLFQs
LFQ <- LFQ_fun( report, "LFQ")

# For a crude estimate of "total signal", calculate the sum of all precursors per sample
int <- report[, .(Intensity = sum( Precursor.Normalised )), .(Protein.Group, Run) ]

# For annotation, count the precursors per protein and sample
counts_1 <- report[, .(N_precursors = uniqueN(Precursor.Id)),        .(Protein.Group, Run) ]
counts_2 <- report[, .(N_precursors_proteotypic = sum(Proteotypic)), .(Protein.Group, Run) ]

# Get protein annotations from report file
annotations <- unique( report[, .(Protein.Group, Run, Sample, Condition, Replicate, First.Protein.Description, Genes)] )

# Combine all these data into a single output table
merge_function <- function(x, y){ merge(x, y, by = c("Protein.Group", "Run"), all = TRUE)}    # Define a function for full outer joins
proteins <- Reduce(merge_function, list(LFQ, int, counts_1, counts_2, annotations))                               

# Remove protein-sample pairs without values
proteins <- proteins[ !is.na(Intensity) ]

# Tidy up column names and order
setnames(proteins, old = c("Protein.Group", "First.Protein.Description"), new = c("Protein_group", "Protein_description"))
setcolorder(proteins, neworder = c(2, 1, 10, 11))

# Define sample order
sample_levels <- c(
  paste0("R", 1:3, "_Control"),
  paste0("R", 1:3, "_DMSO"),
  paste0("R", 1:3, "_SJF8240"),
  paste0("R", 1:3, "_Foretinib")
)

# Apply ordered factor levels to Sample column
report[, Sample := factor(Sample, levels = sample_levels)]
proteins[, Sample := factor(Sample, levels = sample_levels)]

# Write out protein result table
#fwrite(proteins, "protein_quant_new.csv.gz")


#### QC plotting section ####

# Sort samples by desired plotting order
report[, Sample := factor(Sample, levels = sort(unique(Sample))) ]
proteins[, Sample := factor(Sample, levels = sort(unique(Sample))) ]

# Define overall plotting theme
theme_set( theme( line = element_line(linewidth = 0.25), text = element_text(size = 7), rect = element_rect(colour = "black", linewidth = 0.25, fill = NA),
                  panel.border = element_rect(linewidth = 0.25), panel.background = element_blank(), panel.grid = element_blank(),
                  axis.text = element_text(colour = "black", size = 6),
                  legend.background = element_blank(), legend.box.background = element_blank(), legend.key.width = unit(4, "mm"), legend.key = element_rect(colour = NA),
                  strip.background = element_blank(), strip.text = element_text(vjust = 1)))

# (1) Combine data to plot intensity distributions
plot_int_DT <- rbindlist( list(
  report[, .(Sample, `log2 quantity` = log2(Precursor.Normalised), Type = "Precursor quantity")],
  proteins[, .(Sample, `log2 quantity` = log2(LFQ)                 , Type = "Protein MaxLFQ")],
  proteins[, .(Sample, `log2 quantity` = log2(Intensity)           , Type = "Protein Intensity")]
))
plot_int_DT[, Type := factor(Type, levels = c("Precursor quantity", "Protein MaxLFQ", "Protein Intensity"))]   # Set plotting order

# Plot intensity distributions
pInt <- ggplot(plot_int_DT, aes(x = Sample, y = `log2 quantity`))+
  facet_wrap("Type", scales = "free_x")+
  geom_boxplot(outliers = FALSE, linewidth = 0.25)+
  theme(axis.title.y = element_blank())+
  coord_flip()
plot(pInt)

# (2) Count precursors, peptides and proteins to plot their numbers
count_cols <- c("Precursor.Id", "Stripped.Sequence", "Protein.Group")     # Define the columns with the unique features to be counted
counts <- report[ , lapply(.SD, uniqueN), Sample, .SDcols = count_cols ]  # For each sample, count the number of precursors, peptides and proteins

# Plot precursor, peptide and protein counts
facet_labels <- as_labeller( c(Precursor.Id = "Precursors", Stripped.Sequence = "Peptides", Protein.Group = "Protein groups")) # Create plot labellers for facets

pCount <- ggplot( data = melt(counts, id.vars = "Sample", value.name = "IDs"),
                  aes(x = Sample, y = IDs/1000, label = format(IDs, big.mark = ",", scientific = FALSE)))+
  facet_wrap("variable", scales = "free_x", labeller = facet_labels)+
  geom_bar(alpha = 0.5, stat = "identity")+
  geom_text(size = 2.2, hjust = 1.2)+
  ylab("Number of IDs [x1,000]")+
  theme(axis.title.y = element_blank())+
  coord_flip()
plot(pCount)

# (3) Plot data completeness for different minimum precursor cut-offs

# Define function to count the cumulative number of proteins across samples
f_prot_count <- function(dt, var_name){ 
  tmp <- dt[, .(N_samples = .N), Protein_group][, .(N_proteins = .N), N_samples]
  tmp[ order(-N_samples), cumulative_protein_N := cumsum(N_proteins)]
  data.table(tmp, Precursors = var_name)
}

# Count proteins
data_completeness <- rbind( f_prot_count(proteins                                 , "all" ),
                            f_prot_count(proteins[ N_precursors >= 2 ]            , "min 2"),
                            f_prot_count(proteins[ N_precursors_proteotypic >= 2 ], "min 2 pro.typic"))

# Plot the result
pNAs <- ggplot(data_completeness, aes(x = N_samples, y = cumulative_protein_N/1000, colour = Precursors))+
  geom_point(size = 0.8)+
  geom_line(linewidth = 0.25)+
  xlab("Number of samples")+
  ylab("Number of proteins [x1,000]")+
  scale_x_continuous( breaks = seq(1, 1000, 1))+
  scale_y_continuous( limits = c(0, max(data_completeness$cumulative_protein_N )/1000))+
  theme(panel.grid.major = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
        panel.grid.minor = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"))
plot(pNAs)

# (4a) Plot precursor CV values
CV_function <- function(x){ sd(x) / mean(x) * 100 }  # A custom function to calculate CV values and express them as %
precursor_CVs <- report[, .(PrecursorCV = CV_function(Precursor.Normalised), N = .N), .(Precursor.Id, Condition)] # Calculate precursor CVs across replicates
precursor_CVs <- precursor_CVs[ N >= 3 ]  # Restrict to precursors detected in at least 3 replicates
precursor_CVs[, precursor_rank := frank(PrecursorCV), Condition] # Rank precursors by CV value

pCV1 <- ggplot(precursor_CVs, aes(x = precursor_rank/1000, y = PrecursorCV, colour = Condition))+
  geom_line()+
  labs(x = "Number of precursors [x1,000]", y = "precursor CV [%]")+
  coord_cartesian(ylim = c(0,50))+
  theme(panel.grid.major = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
        panel.grid.minor = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"))

plot(pCV1)

# (4b) Plot protein CV values
protein_CVs <- proteins[, .(LFQ_CV = CV_function(LFQ), N = .N), .(Protein_group, Condition)]  # Calculate protein CVs
protein_CVs <- protein_CVs[ N >= 3 ]  # Restrict to proteins detected in at least 3 replicates
protein_CVs[, protein_rank := frank(LFQ_CV), Condition ] # Rank proteins by their CV value

pCV2 <- ggplot(protein_CVs, aes(x = protein_rank/1000, y = LFQ_CV, colour = Condition))+
  geom_line()+
  labs(x = "Number of proteins [x1,000]", y = "protein CV [%]")+
  coord_cartesian(ylim = c(0,50))+
  theme(panel.grid.major = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
        panel.grid.minor = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"))
plot(pCV2)

# Run a PCA to check for major trends
PCA_input <- dcast(proteins, Protein_group ~ Sample, value.var = "LFQ") # Cast this table into a wide data format for PCA analysis
PCA_input <- PCA_input[ complete.cases(PCA_input) ]                     # Can't do PCA with missing values
PCA_input <- t( data.frame( PCA_input, row.names = "Protein_group" ))   # Convert to dataframe and transpose to compare samples
DT_PCA <- prcomp( PCA_input , scale. = TRUE )                           # Run PCA
var_explained <- summary(DT_PCA)$importance                             # Extract variance explained
DT_PCA <- data.table(DT_PCA$x, keep.rownames = "Sample")                # Prepare for plotting
DT_PCA <- merge(DT_PCA, unique( proteins[, .(Sample, Condition, Replicate)]), by = "Sample")  # Annotate condition and replicates

# Create PCA plot
pPCA <- ggplot(DT_PCA, aes(x = PC1, y = PC2, colour = Condition, shape = Replicate))+
  geom_point(size = 2, alpha = 0.8 )+
  xlab(paste("PC1 [", round( var_explained[ rownames(var_explained) == "Proportion of Variance" , "PC1"] * 100, 0), "%]", sep = ""))+
  ylab(paste("PC2 [", round( var_explained[ rownames(var_explained) == "Proportion of Variance" , "PC2"] * 100, 0), "%]", sep = ""))+
  theme( axis.text = element_blank(), axis.ticks = element_blank(), legend.title = element_blank())
plot(pPCA)

# Combine plots using patchwork
custom_plot_layout <- "AAAA
                       BBBB
                       DDEE
                       CCFF"
combined_plot <- pInt + pCount + free(pNAs, type = "label") + free(pCV1, type = "label") + pCV2 + free(pPCA, type = "label") + plot_layout(design = custom_plot_layout)
combined_plot  # Print the plot on screen
ggsave("protein_quant_plot_new.pdf", combined_plot, width = 18, height = 20, units = "cm", bg = "white")


