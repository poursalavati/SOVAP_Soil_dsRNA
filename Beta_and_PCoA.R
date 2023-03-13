## Beta diversity refers to the differences in species composition among different samples or communities.
## It is often used to study the variation in biodiversity across different environmental conditions or spatial scales.
## In microbial ecology, beta diversity is often calculated using taxonomic data from different samples.

## Principal Coordinate Analysis (PCoA) is a multivariate statistical technique that can be used to visualize and explore 
## the patterns of beta diversity among samples. PCoA transforms the original data into a lower dimensional space, 
## where each axis (principal coordinate) captures a proportion of the total variation in the data. 
## By plotting the samples in this reduced dimensional space, one can visualize the similarity or dissimilarity between
## different samples based on their taxonomic composition.

## in order to perform PCoA on beta diversity data, it is common to first calculate a Bray-Curtis dissimilarity matrix 
## based on the abundance data of each taxon across all samples. The resulting distance matrix is then used as input for 
## the PCoA algorithm, which projects the data into a lower-dimensional space while preserving the distance relationships 
## between the samples. The resulting PCoA plot can provide insights into the underlying patterns of variation in the data.

## To prepare data for PCoA analysis, one needs to merge all taxonomic data from different samples into a single table.
## The resulting table should have the same taxonomic categories (i.e., columns) for each sample, and the rows should
## represent the relative abundance or presence/absence of different taxa. In the example header given, the taxonomic
## categories are labeled as "Taxon", and the samples are labeled by their sample IDs (e.g., CCC.dsRNA_1, CCC.dsRNA_2, etc.).
## Each cell in the table represents the relative abundance or presence/absence of a given taxon in a given sample. Once the
## table is prepared, it can be used as input for PCoA or other beta diversity analyses:

## Taxon,CCC.dsRNA_1,CCC.dsRNA_2,CCC.dsRNA_3,CCC.dsRNA_4,RPT.dsRNA_1,RPT.dsRNA_2,RPT.dsRNA_3,RPT.dsRNA_4,VANA_1,VANA_2,VANA_3,VANA_4

## An easy way to merge and prepare the table:

# Load required libraries
library(dplyr)
library(purrr)

# Get a list of file names matching the pattern "*_Abundance.tsv" using list.files() function
file_list <- list.files(pattern = "*_Abundance.tsv")

# Use lapply() function to read each file in the file list using read.table() function, 
# which returns a list of data frames.
# header = TRUE means the first row is a header row and sep = "\t" indicates tab-separated file.
# Then use reduce() function to merge all the data frames in the list by "Taxon" column 
# using a full join (i.e., keep all rows from all tables) and return a single data frame.
df_merged <- file_list %>% 
  lapply(read.table, header = TRUE, sep = "\t") %>% 
  reduce(full_join, by = "Taxon")

# Write the merged data frame to a tab-separated file "merged_data.tsv" 
# using write.table() function with row.names = FALSE to exclude row names from output.
write.table(df_merged, file = "merged_data.tsv", sep = "\t", row.names = FALSE)

## Calculate beta diversity and plot PCoA data

# Load libraries
library(vegan)
library(ggplot2)

# Read the data
df <- read.csv("merged_data.tsv", header = TRUE, row.names = 1)

# Create a list of the taxonomic levels and prefixes.
tax_levels <- c("Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_prefixes <- c("r__", "k__", "p__", "c__", "o__", "f__", "g__", "s__")

generate_plot <- function(tax_level, tax_prefix) {
# Subset the data for the taxonomic level
  subset_data <- df[grep(tax_prefix, rownames(df)), ]
# Remove taxonomic prefix from row names
  rownames(subset_data) <- sub(paste0("^", tax_prefix), "", rownames(subset_data))
# Transpose the data frame
  subset_data_transposed <- t(subset_data)
# Calculate Bray-Curtis dissimilarity matrix
  d <- vegdist(subset_data_transposed, method = "bray")
# Perform PCoA
  pcoa <- cmdscale(d, k = min(2, ncol(subset_data_transposed) - 1), eig = TRUE)
# Create a data frame with PCoA scores and sample IDs
  pcoa_df <- data.frame(pcoa$points, Sample = row.names(subset_data_transposed))
# Split Sample column into two columns: Sample Type (C, R, or V) and Replicate Number (1-4)
  pcoa_df <- pcoa_df %>%
    separate(Sample, into = c("Viromics_Method", "Replicate"), sep = "_") %>%
    mutate(Replicate = as.factor(Replicate))
# Set color palette for each sample type
  palette <- c("CCC.dsRNA" = "#4CAF50", "RPT.dsRNA" = "#2196F3", "VANA" = "#FFC107")
# Plot the PCoA
  plot <- ggplot(pcoa_df, aes(x = `X1`, y = `X2`, color = Viromics_Method)) +
    geom_point(size = 5, stroke = 2, alpha = .5, fill = "white") +
    scale_color_manual(values = palette) +
    ggtitle(paste0("PCoA plot - ", tax_level)) +
    xlab(paste0("PC1 (", round(pcoa$eig[1]/sum(pcoa$eig)*100, 1), "%)")) +
    ylab(paste0("PC2 (", round(pcoa$eig[2]/sum(pcoa$eig)*100, 1), "%)")) +
    theme(legend.position = "bottom")
	ggsave(paste(tax_level, "PCoA_plot.svg"), plot, width = 8, height = 8, units = "in")
  return(plot)
}

plots <- list()

for (i in 1:length(tax_levels)) {
  plot <- generate_plot(tax_levels[i], tax_prefixes[i])
  plots[[i]] <- plot
}
