#!/usr/bin/env Rscript

# This script generates a phylogenetic tree and heatmaps to visualize taxonomic differences at species, genus, and phylum levels.
# The input requires a combined Bracken report and a ground truth species file, both in CSV format.

# Usage: Rscript script.R -r/--reports <path/to/combined_reports.csv> [-t/--ground-truth <path/to/ground_truth.csv>]

# Load necessary libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ape)
  library(taxize)
  library(ggtree)
  library(treeio)
  library(optparse)
  library(pheatmap)
  library(ggtext)
  library(scales)
  library(cowplot)
  library(viridis)
})

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list=list(
  make_option(c("-r", "--reports"), type="character"),
  make_option(c("-t", "--ground-truth"), type="character", default = NULL))))  # make ground-truth optional

# Check if the reports file is provided
if (is.null(opt$reports)) {
  stop("Usage: Rscript script.R -r/--reports <path/to/combined_reports.csv> [-t/--ground-truth <path/to/ground_truth.csv>]")
} else {
  # Load and clean combined Bracken reports
  if (file.exists(opt$reports)) {
    read_counts <- read_csv(opt$reports, show_col_types = FALSE) %>%
      mutate(
        species = recode(species,  # Rename species to reflect current NCBI taxonomy.
                         "Bacillus intestinalis" = "Bacillus spizizenii",
                         "Cryptococcus gattii VGI" = "Cryptococcus gattii",
                         "Cryptococcus gattii VGII" = "Cryptococcus deuterogattii" )) %>%
      replace(is.na(.), 0) %>%  # Replace any NA values with 0s.
      filter(rowMeans(select(., -species)) > 10)  # Filter species with average read count ≤ 10.
    total_reads <- colSums(select(read_counts, -species))  # Total reads per sample (excluding species column).
  } else {
    stop("❌ Combined reports file not found.")
  }
}

# Load ground truth species or create empty table if not provided
if (!is.null(opt$`ground-truth`)) {
  if (file.exists(opt$`ground-truth`)) {
    ground_truth <- read_csv(opt$`ground-truth`, show_col_types = FALSE)
    gt_flag <- TRUE  # Ground truth is present
  } else {
    stop("❌ Ground truth file path provided does not exist.")
  }
} else {
  message("⚠️ No ground truth provided, creating an empty table.")
  ground_truth <- tibble(species = character(), abundance = numeric())   # Create empty tibble with expected columns
  gt_flag <- FALSE  # Ground truth not provided
}

# Get unique ground truth species
gt_species <- unique(ground_truth$species)

# Identify species in read_counts not listed in ground_truth
missing_species <- setdiff(read_counts$species, gt_species)

# Add missing species with zero abundance
if (length(missing_species) > 0) {
  ground_truth <- bind_rows(
    ground_truth,
    tibble(species = missing_species, abundance = 0))
}

# ---------------------------------------
# TREE CONSTRUCTION & ANNOTATION
# ---------------------------------------

# ---------- HELPER FUNCTIONS ----------

# Get annotation positions (x, y) for taxon labels
get_taxa_annot_positions <- function(labels, rank) {
  taxonomy_lookup %>%
    filter(.data[[rank]] %in% labels) %>%                           # Keep only requested taxa
    left_join(tree_data, by = c("species" = "label")) %>%           # Get tree x/y positions
    group_by(tax_label = .data[[rank]]) %>%                         # Group by taxon name
    summarise(
      x = mean(x, na.rm = TRUE) - ifelse(rank == "genus", 3, 8),    # Offset left for label placement
      y = if (rank == "phylum") max(y, na.rm = TRUE) else mean(y, na.rm = TRUE),  # Use top or center y
      .groups = "drop")
}

# Get MRCA node for each taxon group with ≥2 species in the tree
get_taxon_nodes <- function(rank) {
  taxonomy_lookup %>%
    filter(!is.na(.data[[rank]])) %>%                 # Remove missing rank values
    group_by(taxon = .data[[rank]]) %>%               # Group by rank (e.g., genus or phylum)
    summarise(species_list = list(species), .groups = "drop") %>%
    mutate(
      species_in_tree = map(species_list, ~ intersect(.x, tree$tip.label))  # Keep only species present in tree
    ) %>%
    filter(lengths(species_in_tree) >= 2) %>%         # Keep groups with ≥2 valid species
    mutate(node = map_int(species_in_tree, ~ getMRCA(tree, .x))) %>%  # Get the MRCA node
    select(!!rank := taxon, node)                     # Return rank name and node
}

# ---------- MAIN CONSTRUCTION ----------

# Get taxonomy from NCBI (based on species names)
options(ENTREZ_KEY = "5a8133264ac32a3f11c0f1e666a90d96c908")
taxonomy_data <- classification(read_counts$species, db = "ncbi", http_version = 2L)  

# Build tree from taxonomy
tree <- class2tree(taxonomy_data)$phylo  
scale = length(tree$tip.label)

# Initialize base tree plot
tree_plot <- ggtree(tree, layout = "rectangular", size = 0.3, branch.length = "none") +
  theme(axis.ticks = element_blank())

# Match species order to tree tips
species_order <- get_taxa_name(tree_plot)

# Clean and format taxonomy table
taxonomy_lookup <- bind_rows(lapply(taxonomy_data, as_tibble), .id = "id") %>%
  filter(rank %in% c("species", "genus", "phylum")) %>%
  pivot_wider(names_from = rank, values_from = name) %>%
  mutate(species = factor(species, levels = species[order(species_order)]),
         species = as.character(species))

# Extract unique genus and phylum names
genus_order  <- unique(taxonomy_lookup$genus)
phylum_order <- unique(taxonomy_lookup$phylum)
tree_data    <- tree_plot$data

# Compute label coordinates for genus and phylum
genus_annots  <- get_taxa_annot_positions(genus_order, "genus")
phylum_annots <- get_taxa_annot_positions(phylum_order, "phylum")

# Compute MRCA nodes for genus and phylum
highlight_data <- bind_rows(
  get_taxon_nodes("phylum") %>% mutate(group = phylum, rank = "phylum"),
  get_taxon_nodes("genus")  %>% mutate(group = genus,  rank = "genus"))

# Compose tree with labels and highlighted clades
tree_plot <- tree_plot +
  geom_label2(data = genus_annots,  aes(x = x, y = y, label = tax_label),
              fill = "white", size = 3, label.size = 0, family = "mono") +
  geom_label2(data = phylum_annots, aes(x = x, y = y, label = tax_label),
              fill = "white", size = 3, label.size = 0, family = "mono") +
  geom_hilight(data = highlight_data, aes(node = node, fill = group, alpha = rank)) +
  scale_fill_viridis_d(option = "D", name = "Taxonomic Group") +
  scale_alpha_manual(values = c(phylum = 0.15, genus = 0.3)) 

# Optionally color tree tips by ground truth membership
tree_plot <- tree_plot + if (gt_flag) {
  list( geom_tippoint(aes(color = factor(if_else(label %in% gt_species, "Target species", "Other species"))), size = 3),
        scale_color_viridis_d(name = "Category", option = "D"),
        theme(legend.position = "inside", legend.position.inside = c(0.25, 0.9), legend.justification = c(0, 0),
              legend.title = element_text(size = 9, face = "bold"), legend.text  = element_text(size = 8)))
} else {
  theme(legend.position = "none")
}

# ---------------------------------------
# HEATMAP CALCULATIONS & PLOTTING
# ---------------------------------------

# ---------- HELPER FUNCTIONS ----------

# Aggregate expected abundances at a given taxonomic rank
aggregate_abundance <- function(taxonomy_lookup, ground_truth, rank) {
  taxonomy_lookup %>%
    inner_join(ground_truth, by = "species") %>%
    group_by(.data[[rank]]) %>%
    summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")
}

# Compute abundance differences between observed and expected values at a given rank
compute_taxon_differences <- function(taxonomy_lookup, species_differences, rank, order_levels) {
  taxonomy_lookup %>%
    inner_join(species_differences, by = "species") %>%
    group_by(.data[[rank]], variable) %>%
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(!!rank := factor(.data[[rank]], levels = rev(order_levels), ordered = TRUE))
}

# Generate a heatmap for the provided abundance difference data
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
      axis.text.y = element_text(family = "mono", face = "bold"), legend.position = "none")
}

# ---------- MAIN CALCULATIONS ----------

# Aggregate ground truth abundances per genus and phylum
genus_percentages  <- aggregate_abundance(taxonomy_lookup, ground_truth, "genus")
phylum_percentages <- aggregate_abundance(taxonomy_lookup, ground_truth, "phylum")

# Compute species-level differences: observed minus expected
species_differences <- read_counts %>%
  mutate(across(-species, ~ . / total_reads[[cur_column()]] * 100)) %>% # Convert counts to percentages
  inner_join(ground_truth, by = c("species" = "species")) %>%  # Match species
  rename(species = species) %>%
  mutate(across(-c(species, abundance), ~ (.-abundance)))  %>%  # Compute deviation
  select(-abundance) %>%
  mutate(AVERAGE = rowMeans(across(-species), na.rm = TRUE)) %>%  # Compute row-wise average
  pivot_longer(-species, names_to = "variable", values_to = "value") %>%
  mutate(species = factor(species, levels = rev(species_order), ordered = TRUE))  # Order species for heatmap

# Compute genus and phylum-level differences by aggregating species-level values
genus_differences  <- compute_taxon_differences(taxonomy_lookup, species_differences, "genus", genus_order)
phylum_differences <- compute_taxon_differences(taxonomy_lookup, species_differences, "phylum", phylum_order)

# Determine plot bounds and thresholds for color scale
lower_bound <- min(species_differences$value, na.rm = TRUE)
upper_bound <- max(species_differences$value, na.rm = TRUE)
threshold <- abs(lower_bound+0.25*(upper_bound-lower_bound))

# Generate heatmaps for species, genus, and phylum levels
species_heatmap <- plot_heatmap(species_differences, "species", "Species-Level Differences", lower_bound, upper_bound, threshold)
genus_heatmap <- plot_heatmap(genus_differences, "genus", "Genus-Level Differences", lower_bound, upper_bound, threshold)
phylum_heatmap <- plot_heatmap(phylum_differences, "phylum", "Phylum-Level Differences", lower_bound, upper_bound, threshold)

# ---------------------------------------
# FINAL COMPOSITE PLOT
# ---------------------------------------

# Combine phylogenetic tree and heatmaps into a single visualization
main_patch <- (tree_plot + species_heatmap)/(genus_heatmap + phylum_heatmap)
legend <- get_legend(
  plot_heatmap(genus_differences, "genus", "", lower_bound, upper_bound, threshold) +
    theme_void() +
    theme(legend.title = element_text(size = 9, face="bold"),legend.text = element_text(size = 8)))

# Overlay legend at top-left of full canvas
patchwork <- wrap_elements(full = main_patch) +
  inset_element(legend, left = 0.01, bottom = 0.90, right = 0.07, top = 0.99, on_top = TRUE ) +
  plot_annotation( title = "Heatmaps of Taxa Abundance Changes Across Samples", 
                   theme = theme( plot.title = element_text(size = 15, face = "bold")))

# Save combined plot as an image
output_path <- file.path(dirname(opt$reports), "phylo_classification_comparison.png")
ggsave(output_path, plot = patchwork, width = 15, height = 15, dpi = 300)

# FIGURE OUT THE LEGEND 