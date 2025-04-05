# Usage: Rscript b_diversity_heatmap.R <path/to/beta_diversity_matrix.tsv>

# Load required libraries
library(ggplot2)      # Visualization
library(reshape2)     # Data transformation (melt function)
library(viridis)      # Color scale
library(dplyr)        # Data manipulation

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript b_diversity_heatmap.R <path/to/beta_diversity_matrix.tsv>")
file_path <- args[1]
if (!file.exists(file_path)) stop("Error: Beta diversity matrix file not found!")

# Read the file and separate metadata from data
file_content <- readLines(file_path)

metadata_lines <- grep("^#", file_content, value = TRUE)  # Extract metadata lines
df <- read.table(text = file_content[!grepl("^#", file_content)], # Extract data lines (excluding metadata)
                 header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, row.names = 1) %>% 
  {diag(.) <- NA  # Set diagonal values to NA
  .[, ncol(.):1]}  # Reverse column order

# Create a lookup table for correct sample name matching
sample_lookup <- setNames(
  sub("_L001", "", sub(".*/", "", sub("\\.bracken.*", "", metadata_lines))),  # Extract and clean sample names
  sub("^#(\\d+).*", "\\1", metadata_lines))  # Extract numeric indices

# Convert dataframe to a numeric matrix
df_matrix <- as.matrix(suppressWarnings(apply(df, 2, as.numeric)))  
rownames(df_matrix) <- ifelse(rownames(df) %in% names(sample_lookup), sample_lookup[rownames(df)], rownames(df))
colnames(df_matrix) <- ifelse(colnames(df_matrix) %in% names(sample_lookup), sample_lookup[colnames(df_matrix)], colnames(df_matrix))

# Convert matrix to long format for ggplot
b_diversity_melted <- melt(df_matrix, na.rm = TRUE)  # Remove NA values

# Define output file path for the heatmap
heatmap_path <- file.path(dirname(file_path), "beta_diversity_heatmap.png")

# Save heatmap as a high-resolution PNG
size = 200
png(heatmap_path, width = 9*size, height = 8*size, res = 1.5*size)

# Generate the heatmap
ggplot(b_diversity_melted, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +  # Use white grid lines
  scale_fill_viridis_c(option = "D", name = "β-diversity") +  
  geom_text(
    aes(
      label = sprintf("%.2f", value),
      color = ifelse(value > quantile(value, 0.75, na.rm = TRUE), "black", "white")),  # Adjust text color for readability
    size = 3, show.legend = FALSE) +
  scale_color_identity() +  # Use specified text colors
  theme_minimal() +
  theme(
    aspect.ratio=1,
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels 
    panel.grid = element_blank(),  # Remove background grid
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +  # Center and bold title
  labs(x = NULL, y = NULL, title = "β-diversity Heatmap")  # Remove axis titles
dev.off()

# Print confirmation message
cat("Heatmap saved to:", heatmap_path, "\n")