#!/usr/bin/env Rscript

# This script generate a phylogenetic tree from BLAST results, highlighting taxa from human-aligned reads vs. ground truth.

# Usage: Rscript human_aligned_tree.R <path/to/combined_blast_results.txt>

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)    # for data analysis, manipulation, and visualization (loads dplyr, tidyr, ggplot2, readr, stringr, etc.)
  library(ape)          # for reading, writing, plotting, and manipulating phylogenetic trees
  library(taxize)       # to interact with web APIs for taxonomic tasks
  library(treeio)       # to import and store phylogenetic tree with associated data
  library(ggtree)       # for visualization and annotation of phylogenetic trees
  library(argparse)
  library(viridis)
})

# Set NCBI API key for taxize
options(ENTREZ_KEY = "5a8133264ac32a3f11c0f1e666a90d96c908")

# Argument parsing
parser <- ArgumentParser(description = "Karyotype plot script with BLAST results and ground truth")
parser$add_argument("-t", "--ground-truth", required = TRUE, help = "Path to the ground truth file")
parser$add_argument("-b", "--blast-results", required = TRUE, help = "Path to the BLAST results file")
args <- parser$parse_args()
ground_truth <- args$ground_truth
blast_result <- args$blast_results

# Input validation
if (!file.exists(ground_truth)) stop("❌ Ground truth file not found at: ", ground_truth)
if (!file.exists(blast_result)) stop("❌ BLAST results file not found at: ", blast_result)

# Load BLAST results and extract unique taxonomic IDs
taxids <- read_tsv(blast_result, col_names = FALSE, show_col_types = FALSE) %>%
  pull(X2) %>%  # Extract taxonomic IDs from the second column
  unique() %>%  # Remove duplicates
  na.omit()     # Remove missing values

# Retrieve full taxonomy for each taxid
taxonomy_data <- classification(taxids, db = "ncbi")

# Convert taxonomy data to a phylogenetic tree
tree <- class2tree(taxonomy_data)$phylo

# Load ground truth species
highlight_species <- read_csv(ground_truth, show_col_types = FALSE) %>%
  pull(species) %>%
  unique()

# Set dynamic scaling based on number of tips
scale = length(tree$tip.label)

# Create the phylogenetic tree plot
tree_plot <- ggtree(tree, layout = "rectangular", size = 0.3, branch.length="none") +
  geom_tiplab(aes(label = label), fontface = ifelse(tree$tip.label == "Homo sapiens", "bold", "plain"),
              align = FALSE, size = 3, offset = scale / 40) + # Add tip labels, bolding "Homo sapiens"
  geom_tippoint(aes(color = factor(case_when(label == "Homo sapiens" ~ "Human",
                                 label %in% highlight_species ~ "Target species",
                                 TRUE ~ "Other species"), 
                       levels = c("Human", "Target species", "Other species"))),
  size = 3) + # Add colored tip points by category
  scale_color_manual(name = "Read Status", 
                     values = c("Human" =viridis(3, option = "D")[2], 
                                "Target species" = viridis(3, option = "D")[3], 
                                "Other species" = viridis(3, option = "D")[1])) +
  labs(x = NULL, y = NULL, title = "Taxonomic Tree of Taxa Identified by BLASTing Human-Aligned Reads") +
  theme_tree() +
  xlim(0, scale*1.3) +  # Expand space for labels
  theme(axis.ticks = element_blank(),
        plot.title = element_text(size = 10, face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   
        plot.background = element_rect(fill = "white", color = NA), 
        legend.key = element_rect(fill = NA, color = NA),
        legend.position = "inside", legend.position.inside =c(0.01, 0.99), legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("gray", 0.4), color = NA),
        legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) 

# Save the tree as a PNG image
ggsave(file.path(dirname(blast_result), "combined_blast_tree.png"), tree_plot, width =6, height = 5, dpi = 300)
