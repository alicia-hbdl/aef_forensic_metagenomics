# Usage: Rscript jaccard_similarity_heatmap.R <path/to/jaccard_results.txt>

# Load required libraries
library(ggplot2)      # Visualization
library(reshape2)     # Data transformation (melt function)
library(viridis)      # Color scale
library(stringr)      # String manipulation (mutate function)
library(dplyr)        # Data manipulation

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript jaccard_similarity_heatmap.R <path/to/jaccard_results.txt>")
file_path <- args[1]
if (!file.exists(file_path)) stop("Error: Jaccard result file not found!")

# Read Jaccard similarity data and clean sample names
jaccard_data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(
    Sample1 = str_remove_all(Sample1, "_L001"),  # Remove "_L001" from sample names
    Sample2 = str_remove_all(Sample2, "_L001"))

# Extract unique sample names
samples <- sort(unique(c(jaccard_data$Sample1, jaccard_data$Sample2)), decreasing=T)

# Initialize a symmetric Jaccard similarity matrix with NA values
jaccard_matrix <- matrix(NA, nrow = length(samples), ncol = length(samples), dimnames = list(samples, samples))

# Populate the matrix with Jaccard similarity values
for (i in 1:nrow(jaccard_data)) {
  s1 <- jaccard_data$Sample1[i]
  s2 <- jaccard_data$Sample2[i]
  jaccard_matrix[s1, s2] <- jaccard_data$Jaccard[i]  # Assign value to (s1, s2)
  jaccard_matrix[s2, s1] <- jaccard_data$Jaccard[i]}  # Ensure symmetry

# Mask the upper triangle and set diagonal to NA
for (i in 1:nrow(jaccard_matrix)) {
  for (j in i:ncol(jaccard_matrix)) {
    jaccard_matrix[i, j] <- NA}}  # Mask upper triangle

# Flip rows for better visualization
jaccard_matrix <- jaccard_matrix[nrow(jaccard_matrix):1, ]

# Convert matrix to long format for ggplot
jaccard_melted <- melt(jaccard_matrix, na.rm = TRUE)  # Remove NA values

# Define output file path for the heatmap
heatmap_path <- file.path(dirname(file_path), "jaccard_similarity_heatmap.png")

# Save heatmap as a high-resolution PNG
size = 200
png(heatmap_path, width = 9*size, height = 8*size, res = 1.5*size)

# Generate the heatmap
ggplot(jaccard_melted, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") + 
  scale_fill_viridis_c(option = "A", name = "JSI") +  
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
  labs(x = NULL, y = NULL, title = "Jaccard Similarity Heatmap")  # Remove axis titles
dev.off()

# Print confirmation message
cat("Heatmap saved at:", heatmap_path, "\n")
