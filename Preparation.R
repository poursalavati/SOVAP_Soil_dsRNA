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
  colnames(tpm) <- c("ID", "Count", "CPM", "TPM", "FPKM", "Length")
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
