#!/usr/bin/env Rscript

# Generate a β-diversity heatmap from a TSV matrix with optional metadata headers (#).

# Usage: Rscript b_diversity_heatmap.R <path/to/beta_diversity_matrix.tsv>

# Load necessary libraries
suppressPackageStartupMessages({
  library(tidyverse)    # for data analysis, manipulation, and visualization (loads dplyr, tidyr, ggplot2, readr, stringr, etc.)
  library(reshape2)     # to restructure and aggregate data
  library(viridis)      # for generating the color maps
})

# Parse and validate input
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript b_diversity_heatmap.R <matrix.tsv>")
file_path <- args[1]
if (!file.exists(file_path)) stop("Error: File not found!")

# Parse and validate input
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript b_diversity_heatmap.R <matrix.tsv>")
file_path <- args[1]
if (!file.exists(file_path)) stop("Error: File not found!")

# Read full file lines (to separate metadata from data)
lines <- readLines(file_path)

# Extract metadata lines (starting with "#") for sample ID lookup
metadata <- grep("^#", lines, value = TRUE)

# Read β-diversity matrix as a data frame
data <- read.table(header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, 
                   text = lines[!grepl("^#", lines)], # Skips metadata lines
                   row.names = 1, # Keeps row names
                   na.strings = "x.xxx") # Treats 'x.xxx' strings as NA

# Set diagonal to NA (self-comparisons) and reverse column order 
diag(data) <- NA
data <- data[, rev(seq_len(ncol(data)))]

# Build a lookup table to map numeric IDs to cleaned sample names from metadata
lookup <- setNames(
  sub("_L001", "", sub(".*/", "", sub("\\.bracken.*", "", metadata))),
  sub("^#(\\d+).*", "\\1", metadata)
)

# Convert data to numeric matrix and apply cleaned sample names as dimnames
matrix <- apply(data, 2, as.numeric) %>% `dimnames<-`(list(
  ifelse(rownames(data) %in% names(lookup), lookup[rownames(data)], rownames(data)),
  ifelse(colnames(data) %in% names(lookup), lookup[colnames(data)], colnames(data))
))

# Melt the matrix into long format for ggplot
df_long <- melt(matrix, na.rm = TRUE)

# Create heatmap plot object
heatmap <- ggplot(df_long, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "D", name = "β-diversity") + 
  labs(x = NULL, y = NULL, title = "Jaccard Similarity Heatmap") +
  theme_minimal() +
  theme(aspect.ratio=1, plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   
        plot.background = element_rect(fill = "white", color = NA),
        legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8))

# Set plot size: 0.5 inch per sample, minimum 3
size <- max(3, (0.5 * ncol(data)))

# Save heatmap 
ggsave(file.path(dirname(file_path), "beta_diversity_heatmap.png"), heatmap, width = size, height = size, dpi = 300)
