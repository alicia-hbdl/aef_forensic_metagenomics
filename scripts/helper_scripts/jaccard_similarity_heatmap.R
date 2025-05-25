#!/usr/bin/env Rscript

# This script generates a Jaccard similarity heatmap from a tab-separated input file.
# The file must contain three columns: Sample1, Sample2, and Jaccard (as output by `bedtools jaccard`).

# Usage: Rscript jaccard_similarity_heatmap.R <path/to/jaccard_results.txt>

# Load required libraries silently
suppressPackageStartupMessages({
  library(tidyverse)   # Data manipulation and visualization (dplyr, ggplot2, etc.)
  library(viridis)     # Colorblind-friendly color scales
})

# Parse and validate input file path
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript jaccard_similarity_heatmap.R <path/to/jaccard_results.txt>")
file_path <- args[1]
if (!file.exists(file_path)) stop("Error: Jaccard result file not found!")

# Read and clean Jaccard data
jaccard_data <- read_tsv(file_path, show_col_types = FALSE) %>%
  mutate(across(c(Sample1, Sample2), ~ str_remove_all(.x, "_L001|_human"))) # Remove unwanted suffixes from sample names

# Extract sorted list of unique sample names
samples <- sort(unique(c(jaccard_data$Sample1, jaccard_data$Sample2)), decreasing = TRUE)

# Create full symmetric data
full_data <- jaccard_data %>%
  bind_rows(jaccard_data %>% rename(Sample1 = Sample2, Sample2 = Sample1)) %>% # Add mirrored pairs
  complete(Sample1 = samples, Sample2 = samples) %>% # Ensure all combinations exist
  filter(Sample1 != Sample2 & match(Sample1, samples) > match(Sample2, samples))  # Lower triangle only

# Determine upper limit for color scale (minimum 0.01)
scale_max <- max(0.01, max(jaccard_data$Jaccard, na.rm = TRUE))

# Build the heatmap
heatmap <- ggplot(full_data, aes(fct_rev(Sample2), Sample1, fill = Jaccard)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "D", name = "JSI", limits = c(0, scale_max)) +  # Use color scale fixed at 0–scale_max
  labs(x = NULL, y = NULL, title = "Jaccard Similarity Heatmap") +
  theme_minimal() +
  theme(aspect.ratio=1, plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   
        plot.background = element_rect(fill = "white", color = NA),
        legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8))

# Set plot size: 0.5 inch per sample, minimum 3
size <- max(3, (0.5 * length(samples)))
size
# Save the heatmap as a PNG in the same directory as the input file
ggsave(file.path(dirname(file_path), "jaccard_similarity_heatmap.png"), heatmap, width = size, height = size, dpi = 300)