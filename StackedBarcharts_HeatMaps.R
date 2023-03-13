## The code loads several packages including tidyverse, data.table, reshape2, dplyr, purrr and ggplot2. 
## Then it sets the folder path and reads in all files with the pattern "*_sum.tsv" and joins them 
## based on "Taxon" column, saving the merged data to a new file.

## The code then creates two options of the generate_plot function that generate a stacked bar chart for each taxonomic level. 
## The first option generates a basic stacked bar chart for each taxonomic level, while the second option scales the y-axis to 100 percent
## on each plot. The code then creates a list of plots for each taxonomic level using a loop and another function generate_heatmap 
## that generates a heatmap for each taxonomic level.

# Load required packages
library(tidyverse)
library(data.table)
library(reshape2)
library(dplyr)
library(purrr)
library(ggplot2)

# Change to the working directory and get a list of all files with the pattern "*_sum.tsv"
file_list <- list.files(pattern = "*_sum.tsv")

# Read in each file and join them based on "Taxon" column
df_merged <- file_list %>% 
  lapply(read.table, header = TRUE, sep = "\t") %>% 
  reduce(full_join, by = "Taxon")

# Write the merged data to a new file
write.table(df_merged, file = "merged_data.tsv", sep = "\t", row.names = FALSE)

#Create a list of the taxonomic levels and prefixes.
tax_levels <- c("Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_prefixes <- c("r__", "k__", "p__", "c__", "o__", "f__", "g__", "s__")

#OPT1: Create a function to generate a stacked bar chart for each taxonomic level.
generate_plot <- function(tax_level, tax_prefix) {
  # Subset the data for the taxonomic level
  subset_data <- subset(data, grepl(tax_prefix, Taxon))
  
  # Melt the data using the `melt` function from the `reshape2` library
  melted_data <- melt(subset_data, id.vars = "Taxon", variable.name = "Sample", value.name = "TPM")
  
  # Create the stacked bar chart using the `ggplot` function
  plot <- ggplot(melted_data, aes(x = Sample, y = TPM, fill = Taxon)) + 
    geom_bar(stat = "identity") +
    labs(title = paste(tax_level, "TPM values for each sample"), x = "Sample", y = "TPM") +
    theme_bw()
  
  return(plot)
}

#OPT2: To scale the y-axis to 100 percent on each plot, you can add the scale_y_continuous function to the generate_plot function as follows:
generate_plot <- function(tax_level, tax_prefix) {
  # Subset the data for the taxonomic level
  subset_data <- subset(data, grepl(tax_prefix, Taxon))
  subset_data$Taxon <- gsub(paste0("^", tax_prefix), "", subset_data$Taxon)
  # Melt the data using the `melt` function from the `reshape2` library
  melted_data <- melt(subset_data, id.vars = "Taxon", variable.name = "Sample", value.name = "TPM")
  
  # Create the stacked bar chart using the `ggplot` function
  plot <- ggplot(melted_data, aes(x = Sample, y = TPM/sum(TPM), fill = Taxon)) + 
    geom_bar(position = "fill", stat = "identity", color = "gray", linewidth=0.3) +
    labs(title = paste(tax_level, "TPM values for each sample"), x = "Sample", y = "Percentage") +
    theme_bw() + scale_fill_manual(values = colours) +
    scale_y_continuous(labels = scales::percent)
    ggsave(paste(tax_level, "stacked_barplot.svg"), plot, width = 8, height = 8, units = "in")
  
  return(plot)
}

#Generate a list of the plots for each taxonomic level using a loop.
plots <- list()

for (i in 1:length(tax_levels)) {
  plot <- generate_plot(tax_levels[i], tax_prefixes[i])
  plots[[i]] <- plot
}

# Define a function to generate a heatmap
generate_heatmap <- function(tax_level, tax_prefix) {
  # Subset the data for the taxonomic level
  subset_data <- subset(data, grepl(tax_prefix, Taxon))
  # Remove the prefix from the taxon names
  subset_data$Taxon <- gsub(paste0("^", tax_prefix), "", subset_data$Taxon)
    # Melt the data using the `melt` function from the `reshape2` library
  melted_data <- melt(subset_data, id.vars = "Taxon", variable.name = "Sample", value.name = "TPM")
    # Calculate the sum of TPM values for each sample
  sum_tpm <- aggregate(melted_data$TPM, list(melted_data$Sample), sum)
  names(sum_tpm) <- c("Sample", "sum_tpm")
    # Merge the sum of TPM values with the melted data
  melted_data <- merge(melted_data, sum_tpm)
    # Convert TPM to percentage scale
  melted_data$TPM_percent <- melted_data$TPM / melted_data$sum_tpm * 100
    # Create the heatmap using the `ggplot` function
  plot <- ggplot(melted_data, aes(x = Sample, y = Taxon, fill = TPM_percent)) + 
    geom_tile(colour = "gray") +
    geom_text(aes(label = paste(round(TPM_percent, 1))), size = 2.5, colour = "black") +
    labs(title = paste(tax_level, "heatmap"), x = "Sample", y = "Taxon") +
    scale_fill_gradient2(low = "#1B1E09", mid = "#ffffff", high = "#C7E63E",
                         na.value = "#FFFFFF",
                         name = "TPM (%)",
                         limits = c(0, 100),
                         labels = scales::comma_format()) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 0))
    ggsave(paste(tax_level, "heatmap_plot.svg"), plot, width = 8, height = 8, units = "in") 
  return(plot)
}
plots <- list()

for (i in 1:length(tax_levels)) {
  plot <- generate_heatmap(tax_levels[i], tax_prefixes[i])
  plots[[i]] <- plot
}
