#This R code is used to manipulate and prepare necessary data and tables to plot taxonomy results in stacked bar charts and heatmaps. 

#The code sets the folder path where the taxonomy and abundance data are stored. 
#The list.files function is used to get a list of all file names with the extensions ".taxo" and ".tsv" in the specified folder.
#The code then loops over each pair of files and reads them in using read.table function.
#The first file contains the abundance data in TPM format, while the second file contains the taxonomy information for each ID. 
#The colClasses parameter is set to "character" to ensure that the taxonomy columns are read as character data. 
#The na.strings parameter is set to "NA" to ensure that any missing values are read in as "NA".
#The fill parameter is set to TRUE to ensure that the shorter rows are padded with "NA" values.
#The code then manipulates the data by renaming the abundance columns to more informative names and extracting only the relevant columns from the taxonomy data. 
#The taxonomy columns are then split using the tstrsplit function, and the resulting columns are added to the taxonomy data. 
#The merged_files are obtained by merging the abundance data and the taxonomy data by the ID column.
#Finally, the code writes the output to a file with the extension "_Abundance.tsv". 
#The output file contains the taxonomy information and abundance data in a tab-separated format that can be used to graph the taxonomy data.

# Load required packages
library(tidyverse)
library(data.table)
library(reshape2)
library(dplyr)
library(purrr)
library(ggplot2)


# Set the folder path
folder_path <- "path/to/taxonomy/and/abundance/data"
# Get a list of all file names with the extensions "*.taxo" and "*.tsv"
file_names <- list.files(folder_path, pattern = ".*\\.(taxo|tsv)$", full.names = TRUE)
# Loop over each pair of files

for (i in seq_along(file_names)) {
  # Read in the files
  tpm <- read.table(file_names[i], header = FALSE, sep = "\t")
  mydata <- read.table(gsub("\\.tsv$", ".taxo", file_names[i]), header = FALSE, sep = "\t", colClasses = "character", na.strings = "NA", fill = TRUE)
  
  # Manipulate the data
  colnames(tpm) <- c("ID", "Count", "CPM", "TPM", "FPKM", "Len")
  mydata <- mydata[, c(1, 2, 16)]
  newcols <- tstrsplit(as.character(mydata[, 3]), ";", fixed = TRUE)
  mydata <- cbind(mydata[, c(1, 2)], newcols)
  colnames(mydata) <- c("HIT", "ID", "Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  merged_files <- merge(mydata, tpm, by="ID")
  
  # Write the output to a file
  write.table(merged_files, gsub("\\.(taxo|tsv)$", "_Abundance.tsv", file_names[i]), sep="\t", row.names=FALSE, quote=FALSE)
}

#This R code reads in a table of taxonomic information and abundance data and calculates the sums of TPM and CPM for each taxonomic level 
#(Realm, Kingdom, Phylum, Class, Order, Family, Genus, and Species). 
#The results are combined into a single table (all_sums) and then written to a CSV file.
#The read.table() function is used to read in the input data. 
#The group_by() and summarize() functions from the dplyr package are used to calculate the sums for each taxonomic level.
#The bind_rows() function is used to combine the results into a single table. Finally, the write.csv() function is used to write the output to a CSV file.

#Insert prepared table
data <- read.table('*_Abundance.tsv', header = TRUE, stringsAsFactors = FALSE, sep="\t")

# calculate sums of TPM and CPM for each taxonomic level
realm_sum <- data %>%  group_by(Realm) %>%  summarize(sum_TPM = sum(TPM), sum_CPM = sum(CPM))
kingdom_sum <- data %>%  group_by(Kingdom) %>%  summarize(sum_TPM = sum(TPM), sum_CPM = sum(CPM))
phylum_sum <- data %>%  group_by(Phylum) %>%  summarize(sum_TPM = sum(TPM), sum_CPM = sum(CPM))
class_sum <- data %>%  group_by(Class) %>%  summarize(sum_TPM = sum(TPM), sum_CPM = sum(CPM))
order_sum <- data %>%  group_by(Order) %>%  summarize(sum_TPM = sum(TPM), sum_CPM = sum(CPM))
family_sum <- data %>%  group_by(Family) %>%  summarize(sum_TPM = sum(TPM), sum_CPM = sum(CPM))
genus_sum <- data %>%  group_by(Genus) %>%  summarize(sum_TPM = sum(TPM), sum_CPM = sum(CPM))
species_sum <- data %>%  group_by(Species) %>%  summarize(sum_TPM = sum(TPM), sum_CPM = sum(CPM))

# combine all sums into a single table
all_sums <- bind_rows(
  realm_sum %>% mutate(level = "Realm", Taxon = Realm),
  kingdom_sum %>% mutate(level = "Kingdom", Taxon = Kingdom),
  phylum_sum %>% mutate(level = "Phylum", Taxon = Phylum),
  class_sum %>% mutate(level = "Class", Taxon = Class),
  order_sum %>% mutate(level = "Order", Taxon = Order),
  family_sum %>% mutate(level = "Family", Taxon = Family),
  genus_sum %>% mutate(level = "Genus", Taxon = Genus),
  species_sum %>% mutate(level = "Species", Taxon = Species)
) %>% select(Taxon, level, sum_TPM, sum_CPM)

# write the output to a CSV file
write.csv(all_sums, "*_Abundance_sums.csv", row.names = FALSE)


# Set the folder path and get a list of all files with the pattern "*_sum.tsv"
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
