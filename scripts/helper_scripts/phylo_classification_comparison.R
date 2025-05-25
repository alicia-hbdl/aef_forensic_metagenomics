#!/usr/bin/env Rscript

# This script generates a phylogenetic tree and heatmaps to visualize taxonomic differences at species, genus, and phylum levels.
# The input requires a combined Bracken report and a ground truth species file, both in CSV format.

# Usage: Rscript script.R -r/--reports <path/to/combined_reports.csv> [-t/--ground-truth <path/to/ground_truth.csv>]

# Load necessary libraries
suppressPackageStartupMessages({
  library(tidyverse)    # for data analysis, manipulation, and visualization (loads dplyr, tidyr, ggplot2, readr, stringr, etc.)
  library(patchwork)    # for combining multiple plots
  library(ape)          # for reading, writing, plotting, and manipulating phylogenetic trees
  library(taxize)       # to interact with web APIs for taxonomic tasks
  library(ggtree)       # for visualization and annotation of phylogenetic trees
  library(treeio)       # to import and store phylogenetic tree with associated data
  library(optparse)     # to write "#!" shebang scripts that accept short and long flag/options
  library(pheatmap)     # to implement heatmaps with control over dimensions and appearance  
  library(ggtext)       # to render of complex formatted plot labels
  library(scales)       # to map data to aesthetics
  library(cowplot)
})

# Set NCBI API key for taxonomic queries
options(ENTREZ_KEY = "5a8133264ac32a3f11c0f1e666a90d96c908")

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list=list(
  make_option(c("-r", "--reports"), type="character"),
  make_option(c("-t", "--ground-truth"), type="character", default = NULL)  # make ground-truth optional
)))

# Check if the reports file is provided
if (is.null(opt$reports)) {
  stop("Usage: Rscript script.R -r/--reports <path/to/combined_reports.csv> [-t/--ground-truth <path/to/ground_truth.csv>]")
}


# Load and clean combined Bracken reports
if (file.exists(opt$reports)) {
  read_counts <- read_csv(opt$reports, show_col_types = FALSE) %>%
    mutate(
      species = recode(species,
                       "Bacillus intestinalis" = "Bacillus spizizenii",
                       "Cryptococcus gattii VGI" = "Cryptococcus gattii",
                       "Cryptococcus gattii VGII" = "Cryptococcus deuterogattii"
      )
    ) %>%
    replace(is.na(.), 0) %>%
    filter(rowMeans(select(., -species)) > 10)
  
  total_reads <- colSums(select(read_counts, -species))  # Fix missing total_reads
} else {
  stop("❌ Combined reports file not found.")
}

# Load ground truth species if it exists, or assign an empty table
if (file.exists(opt$`ground-truth`)) {
  ground_truth <- read_csv(opt$`ground-truth`, show_col_types = FALSE)
  highlight_species <- ground_truth$species
  
  # Find species present in read_counts but missing in ground_truth
  missing_species <- setdiff(read_counts$species, highlight_species)
  
  # Add missing species to ground_truth with 0 abundance
  if (length(missing_species) > 0) {
    ground_truth <- bind_rows(
      ground_truth,
      tibble(species = missing_species, abundance = 0)
    )
  }
} else {
  # If no ground truth file, create an empty data frame with correct column names
  ground_truth <- tibble(species = character(), abundance = numeric())
}

# Retrieve taxonomic hierarchy
#tax_ids <- get_uid(read_counts$species)  # Get taxonomic IDs
#taxonomy_data <- classification(tax_ids)  # Fetch taxonomy data

# Directly get taxonomy data (skip tax_id fetching step)
taxonomy_data <- classification(read_counts$species, db = "ncbi", http_version = 2L)  # Fetch taxonomy data directly

# Construct phylogenetic tree
tree <- class2tree(taxonomy_data)$phylo  # Convert classification to tree

# Set dynamic scaling based on number of tips
scale = length(tree$tip.label)

# Plot phylogenetic tree
tree_plot <- ggtree(tree, layout = "rectangular", size = 0.3, branch.length="none") +
  geom_tippoint(aes(color = factor(case_when(label %in% highlight_species ~ "Target species",TRUE ~ "Other species"))), size = 3) +  # Add colored tip points by category
  scale_color_viridis_d(name = "Category", option = "D") +  # Viridis discrete scale
  theme(axis.ticks = element_blank(), 
        legend.position = "inside",
        legend.position.inside = c(0.25, 0.9),
        legend.justification = c(0, 0),
        legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) + 
  xlim(0, max(6, scale/4)) # Adjust x-axis 


# Extract the order of species from the tree
species_order <- get_taxa_name(tree_plot)

# Extract species, genus, and phylum from taxonomy data
taxonomy_lookup <- bind_rows(lapply(taxonomy_data, as_tibble), .id = "id") %>%
  filter(rank %in% c("species", "genus", "phylum")) %>%
  pivot_wider(names_from = rank, values_from = name) %>%
  mutate(species = factor(species, levels = species[order(species_order)])) %>%
  mutate(species = as.character(species))  # Convert back to character

# Extract genus order based on tree order
genus_order <- taxonomy_lookup %>%
  pull(genus) %>%
  unique()

# Extract phylum order based on tree order
phylum_order <- taxonomy_lookup %>%
  pull(phylum) %>%
  unique()

# Aggregate expected percentages at genus level
genus_percentages <- taxonomy_lookup %>%
  inner_join(ground_truth, by = "species") %>%
  group_by(genus) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE))

# Aggregate expected percentages at phylum level
phylum_percentages <- taxonomy_lookup %>%
  inner_join(ground_truth, by = "species") %>%
  group_by(phylum) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE))

# Compute species-level differences between observed and expected abundance
species_differences <- read_counts %>%
  mutate(across(-species, ~ . / total_reads[[cur_column()]] * 100)) %>% # Convert counts to percentages
  inner_join(ground_truth, by = c("species" = "species")) %>%  # Match species
  rename(species = species) %>%
  mutate(across(-c(species, abundance), ~ (.-abundance)))  %>%  # Compute deviation
  select(-abundance) %>%
  mutate(AVERAGE = rowMeans(across(-species), na.rm = TRUE)) %>%  # Compute row-wise average
  pivot_longer(-species, names_to = "variable", values_to = "value") %>%
  mutate(species = factor(species, levels = rev(species_order), ordered = TRUE))  # Order species for heatmap

# Compute genus-level differences
genus_differences <- taxonomy_lookup %>%
  inner_join(species_differences, by = "species") %>%
  group_by(genus, variable) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")%>%
  mutate(genus = factor(genus, levels = rev(genus_order), ordered = TRUE))  # Order genus for heatmap

# Compute phylum-level differences
phylum_differences <- taxonomy_lookup %>%
  inner_join(species_differences, by = "species") %>%
  group_by(phylum, variable) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(phylum = factor(phylum, levels = rev(phylum_order), ordered = TRUE))  # Order phylum for heatmap

# Function to plot heatmaps 
plot_heatmap <- function(data, y_var, title, lower_bound, upper_bound, threshold) {
  ggplot(data, aes(x = variable, y = .data[[y_var]], fill = value)) +
    geom_tile(color = "grey80") +  # Add grid borders
    geom_text(aes(label = sprintf("%+.2f", value)), color = "white", size = 2) +
    scale_fill_viridis_c(option = "D", name = "% Difference", limits = c(lower_bound, upper_bound), oob = squish) + 
    scale_color_identity() +  # Maintain original colors
    labs(title = title, fill = "% Difference") + 
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.title.x = element_blank(), axis.title.y = element_blank(), 
      axis.text = element_text(size = 8),  # ✅ Missing comma added here
      axis.text.x = element_text(angle = 45, hjust = 1, family = "mono", face = "bold"),  
      axis.text.y = element_text(family = "mono", face = "bold"), # if (y_var == "species") element_blank() else 
      legend.position = "none")
}

# Compute global bounds and threshold
lower_bound <- min(species_differences$value, na.rm = TRUE)
upper_bound <- max(species_differences$value, na.rm = TRUE)
threshold <- abs(lower_bound+0.25*(upper_bound-lower_bound))

# Generate heatmaps
species_heatmap <- plot_heatmap(species_differences, "species", "Species-Level Differences", lower_bound, upper_bound, threshold)
genus_heatmap <- plot_heatmap(genus_differences, "genus", "Genus-Level Differences", lower_bound, upper_bound, threshold)
phylum_heatmap <- plot_heatmap(phylum_differences, "phylum", "Phylum-Level Differences", lower_bound, upper_bound, threshold)

# Combine phylogenetic tree and heatmaps into a single visualization
main_patch <- (tree_plot + species_heatmap)/(genus_heatmap + phylum_heatmap)
legend <- get_legend(
  plot_heatmap(genus_differences, "genus", "", lower_bound, upper_bound, threshold) +
    theme_void() +
    theme(legend.background = element_rect(fill = "white", color = NA),
          legend.title = element_text(size = 9, face="bold"),legend.text = element_text(size = 8)))

# Overlay legend at top-left of full canvas
patchwork <- wrap_elements(full = main_patch) +
  inset_element(legend, left = 0.01, bottom = 0.90, right = 0.18, top = 0.99, on_top = TRUE ) +
  plot_annotation( title = "Heatmaps of Taxa Abundance Changes Across Samples", 
                   theme = theme( plot.title = element_text(size = 15, face = "bold")))

# Save combined plot as an image
output_path <- file.path(dirname(opt$reports), "phylo_classification_comparison.png")
ggsave(output_path, plot = patchwork, width = 15, height = 15, dpi = 300)