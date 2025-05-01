#!/usr/bin/env Rscript

# This script generates a phylogenetic tree based on BLAST results, highlighting key species.

# Usage: Rscript human_aligned_tree.R <path/to/combined_blast_results.txt>

# Load necessary libraries
suppressPackageStartupMessages({
  library(tidyverse)    # for data analysis, manipulation, and visualization (loads dplyr, tidyr, ggplot2, readr, stringr, etc.)
  library(ape)          # for reading, writing, plotting, and manipulating phylogenetic trees
  library(taxize)       # to interact with web APIs for taxonomic tasks
  library(treeio)       # to import and store phylogenetic tree with associated data
  library(ggtree)       # for visualization and annotation of phylogenetic trees
})

# Set NCBI API key for taxonomic queries
options(ENTREZ_KEY = "5a8133264ac32a3f11c0f1e666a90d96c908")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("❗ Usage: Rscript human_aligned_tree.R <path/to/combined_blast_results.txt>")
file_path <- args[1]
if (!file.exists(file_path)) stop("❌ Error: BLAST results file not found!")

# Read BLAST results and extract unique taxonomic IDs
taxids <- read_tsv(file_path, col_names = FALSE, show_col_types = FALSE) %>%
  pull(X2) %>%  # Extract taxonomic IDs from the second column
  unique() %>%  # Remove duplicates
  na.omit()     # Remove missing values

# Retrieve full taxonomy classification for taxonomic IDs
taxonomy_data <- classification(taxids, db = "ncbi")

# Convert taxonomy data into a phylogeny tree
tree <- class2tree(taxonomy_data)$phylo

# Define key species to highlight in the tree
highlight_species <- c(
  "Bacillus subtilis", "Enterococcus faecalis", "Escherichia coli", 
  "Limosilactobacillus fermentum", "Listeria monocytogenes", 
  "Pseudomonas aeruginosa", "Salmonella enterica", 
  "Staphylococcus aureus", "Saccharomyces cerevisiae", 
  "Cryptococcus neoformans"
)

# Generate a phylogenetic tree plot
tree_plot <- ggtree(tree, layout = "rectangular", size = 0.3) +
  geom_tiplab(aes(label = label, fontface = "bold",
                  color = ifelse(label %in% highlight_species, "red", "black")),
              align = TRUE, size = 2, family = "mono") +  # Highlight specified species in red
  scale_color_identity() +
  geom_label2(aes(label = label, subset = !isTip), 
              hjust = 1.05, fill = "white", size = 3.5, label.size = 0, family = "mono") +
  theme_tree() +
  theme(axis.ticks = element_blank()) + 
  xlim(-mean(tree$edge.length)/2, max(tree$edge.length)^2/16)

# Save the tree plot
output_path <- file.path(dirname(file_path), "combined_blast_tree.png")
ggsave(output_path, plot = tree_plot, width = 15, height = 22.5, dpi = 300)

# Print confirmation message
cat("✅ Phylogenetic tree saved to:", output_path, "\n")
