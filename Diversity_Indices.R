# This R code reads in four (or more) CSV files containing abundance data for different samples, merges them into one file,
# and calculates richness and diversity indices (shannon and simpson) for each taxonomic level using the vegan package.
# The code defines the taxonomy prefixes and their corresponding level names, reads in the data from the merged CSV file, and initializes
# an empty data frame to store the calculated richness and diversity values. The code then loops through each taxonomy level,
# subsets the data to only include rows with the current taxonomy prefix, calculates the richness and diversity indices for each
# sample separately, and adds the resulting values to the data frame.
# Finally, the code prints the resulting data frame and saves it as a CSV file named "merged_files.csv".

# Load the vegan package for diversity calculation
library(vegan)

# Read in four CSV files containing abundance data
file1 <- read.csv("1_Abundance_sums.csv", header=TRUE, stringsAsFactors=FALSE)
file2 <- read.csv("2_Abundance_sums.csv", header=TRUE, stringsAsFactors=FALSE)
file3 <- read.csv("3_Abundance_sums.csv", header=TRUE, stringsAsFactors=FALSE)
file4 <- read.csv("4_Abundance_sums.csv", header=TRUE, stringsAsFactors=FALSE)

# Merge the four files based on the Taxon column, retaining all rows from all files
merged_files <- merge(merge(merge(file1, file2, by="Taxon", all=TRUE), file3, by="Taxon", all=TRUE), file4, by="Taxon", all=TRUE)

# Write the merged file to a new CSV file
write.csv(merged_files, "merged_files.csv", row.names=FALSE)

# Define the taxonomy prefixes and their corresponding level names
taxa_prefixes <- c("r__", "k__", "p__", "c__", "o__", "f__", "g__", "s__")
taxa_levels <- c("Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Read in the merged data from the CSV file
data <- read.csv("merged_files.csv", header=TRUE, stringsAsFactors=FALSE)

# Initialize an empty data frame to store the richness and diversity values
result_df <- data.frame(Taxon = character(),
                        Level = character(),
                        X1_Richness = numeric(),
                        X2_Richness = numeric(),
                        X3_Richness = numeric(),
                        X4_Richness = numeric(),
                        X1_Shannon = numeric(),
                        X2_Shannon = numeric(),
                        X3_Shannon = numeric(),
                        X4_Shannon = numeric(),
                        X1_Simpson = numeric(),
                        X2_Simpson = numeric(),
                        X3_Simpson = numeric(),
                        X4_Simpson = numeric(),
                        stringsAsFactors = FALSE)

# Loop through each taxonomy level
for (i in 1:length(taxa_prefixes)) {
  # Subset the data to only include rows with the current taxonomy prefix
  prefix_data <- data[grep(paste0("^", taxa_prefixes[i]), data$Taxon), ]
  
  # Calculate the richness for each TPM column separately
  X1_richness <- sum(prefix_data$X1_TPM > 0)
  X2_richness <- sum(prefix_data$X2_TPM > 0)
  X3_richness <- sum(prefix_data$X3_TPM > 0)
  X4_richness <- sum(prefix_data$X4_TPM > 0)
  
  # Calculate the Shannon and Simpson diversity indices for each TPM column separately
  X1_shannon <- diversity(prefix_data$X1_TPM, index = "shannon")
  X2_shannon <- diversity(prefix_data$X2_TPM, index = "shannon")
  X3_shannon <- diversity(prefix_data$X3_TPM, index = "shannon")
  X4_shannon <- diversity(prefix_data$X4_TPM, index = "shannon")
  
  X1_simpson <- diversity(prefix_data$X1_TPM, index = "simpson")
  X2_simpson <- diversity(prefix_data$X2_TPM, index = "simpson")
  X3_simpson <- diversity(prefix_data$X3_TPM, index = "simpson")
  X4_simpson <- diversity(prefix_data$X4_TPM, index = "simpson")

# Add the richness and diversity values to the data frame
  result_df <- rbind(result_df, data.frame(Taxon = taxa_prefixes[i],
                                           Level = taxa_levels[i],
                                           X1_Richness = X1_richness,
                                           X2_Richness = X2_richness,
                                           X3_Richness = X3_richness,
				                                   X4_Richness = X4_richness,
                                           X1_Shannon = X1_shannon,
                                           X2_Shannon = X2_shannon,
                                           X3_Shannon = X3_shannon,
                                           X4_Shannon = X4_shannon,
                                           X1_Simpson = X1_simpson,
                                           X2_Simpson = X2_simpson,
                                           X3_Simpson = X3_simpson,
                                           X4_Simpson = X4_simpson))
}

# Print the resulting data frame
print(result_df)



### If you have replicates for each sample or related samples that you want to group together, you can use the rest of the code to create a violin plot.
## Make sure your dataframe have this header for instance:
## Level	CCC_dsRNA1	CCC_dsRNA2	CCC_dsRNA3	CCC_dsRNA4	RPT_dsRNA1	RPT_dsRNA2	RPT_dsRNA3	RPT_dsRNA4	VANA1	VANA2	VANA3	VANA4

# Load required libraries
library(ggplot2)
library(tidyr)

# Read in the data from a CSV file named "M-rich.csv"
df <- read.csv("Diversity.csv")

# Reshape the data from wide to long format using pivot_longer() function
df_long <- pivot_longer(df, cols = -Level, names_to = c("Sample", "Replicate"), names_sep = "(?<=[A-Za-z])(?=[0-9])")

# Set the order of the levels for plotting
level_order <- c("Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
df_long$Level <- factor(df_long$Level, levels = level_order)

# Create the plot using ggplot2
ggplot(df_long, aes(x = Sample, y = value, fill = Sample)) +
geom_violin() + # Add violin plot layer
geom_dotplot(aes(group = interaction(Sample, Replicate)), binaxis = "y", stackdir = "center", dotsize = 0.8) + # Add dot plot layer
geom_boxplot(aes(group = Sample), width = 0.1, alpha = 0.5, color = "black") + # Add box plot layer
facet_wrap(~Level, scales = "free", ncol = 4) + # Create faceted plot based on level
scale_fill_manual(values = c("#009E73", "#00739E", "#E69F00")) + # Set custom fill colors
labs(x = "Sample", y = "Value") + # Set axis labels
theme_bw() # Set plot theme

# Save the plot as an SVG file named "myplot.svg" with dimensions of 15 inches by 8 inches
ggsave("myplot.svg", width = 15, height = 8)

