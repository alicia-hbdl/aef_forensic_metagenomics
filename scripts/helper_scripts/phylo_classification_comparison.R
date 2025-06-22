#!/usr/bin/env Rscript

# This script generates a phylogenetic tree and heatmaps to visualize taxonomic differences at species, genus, and phylum levels.
# The input requires a combined Bracken report and a ground truth species file, both in CSV format.

# Usage: Rscript script.R -r/--reports <path/to/combined_breports.csv> [-t/--ground-truth <path/to/ground_truth.csv>]
# Note: This script requires internet access to query NCBI via taxize::classification(). Run locally to avoid HPC timeouts.

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
                         "Bacillus subtilis" = "Bacillus spizizenii",
                         "Bacillus intestinalis" = "Bacillus spizizenii",
                         "Cryptococcus gattii VGI" = "Cryptococcus gattii",
                         "Lactobacillus fermentum" = "Limosilactobacillus fermentum",
                         "Cryptococcus gattii VGII" = "Cryptococcus deuterogattii")) %>%
      replace(is.na(.), 0) %>%
      mutate(total_abundance = rowSums(across(-species))) %>%
      slice_max(total_abundance, n = min(55,  nrow(.))) %>%  # Handle fewer than 80 species gracefully
      select(-total_abundance)
    
    total_reads <- colSums(select(read_counts, -species))  # Total reads per sample (excluding species column).
  } else {
    stop("❌ Combined reports file not found.")
  }
}

# Load ground truth species or create empty table if not provided
if (!is.null(opt$`ground-truth`)) {
  if (file.exists(opt$`ground-truth`)) {
    ground_truth <- read_csv(opt$`ground-truth`, show_col_types = FALSE) %>%
      mutate(
        species = recode(species,  # Rename species to reflect current NCBI taxonomy.
                         "Bacillus subtilis" = "Bacillus spizizenii",
                         "Lactobacillus fermentum" = "Limosilactobacillus fermentum" ))
    gt_flag <- TRUE  # Ground truth is present
  } else {
    stop("❌ Ground truth file path provided does not exist.")
  }
} else {
  message("⚠️ No ground truth provided, creating an empty table.")
  ground_truth <- tibble(species = character(), abundance = numeric())   # Create empty tibble with expected columns
  gt_flag <- FALSE  # Ground truth not provided
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
      x = if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE) - ifelse(rank == "genus", 3, 8),
      y = if (all(is.na(y))) NA_real_ else if (rank == "phylum") max(y, na.rm = TRUE) else mean(y, na.rm = TRUE),
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

# Retrieve taxonomy data from NCBI for species in read_counts
taxonomy_full <- classification(read_counts$species, db = "ncbi")

# Filter to keep only entries that are data frames and contain valid species, genus, and phylum information
taxonomy_data <- keep(taxonomy_full, ~ 
                        is.data.frame(.) &&
                        all(c("species", "genus", "phylum") %in% .$rank) &&
                        all(!is.na(.$name[.$rank %in% c("species", "genus", "phylum")]))
)

# Initialize dropped_text with a default message
dropped_text <- "No species were dropped due to incomplete taxonomy."

# Report any species dropped due to incomplete taxonomy
dropped_species <- setdiff(names(taxonomy_full), names(taxonomy_data))
if (length(dropped_species) > 0) {
  dropped_text <- paste0("Dropped species due to incomplete taxonomy: ", 
                         paste(dropped_species, collapse = ", "))
} 
cat(dropped_text)

#  Extract updated species names from taxonomy and map back to read_counts
taxonomy_updates <- map2_dfr(taxonomy_data, names(taxonomy_data), ~ {
  new_name <- .x$name[.x$rank == "species"]
  tibble(original = .y, updated = new_name)
})

# Update species names in read_counts using the NCBI-corrected taxonomy
read_counts <- read_counts %>%
  left_join(taxonomy_updates, by = c("species" = "original")) %>%
  mutate(species = coalesce(updated, species)) %>%
  select(-updated) %>%
  filter(!species %in% dropped_species) %>%  # Drop species with incomplete taxonomy
  group_by(species) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%  # Merge counts for duplicate species
  ungroup()

# Print which species names were updated
name_changes <- taxonomy_updates %>%
  filter(original != updated) %>%
  arrange(original) 

if (nrow(name_changes) > 0) {
  cat("🔄 Updated species names:\n")
  name_changes %>%
    mutate(change = paste(original, "→", updated)) %>%
    pull(change) %>%
    cat(sep = "\n")
} else {
  cat("✅ No species names were updated.\n")
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

# Build tree from taxonomy
tree <- class2tree(taxonomy_data)$phylo  
scale = length(tree$tip.label)

# Initialize base tree plot
base_tree_plot <- ggtree(tree, layout = "rectangular", size = 0.3, branch.length = "none") +
  theme(axis.ticks = element_blank())

# Match species order to tree tips
species_order <- get_taxa_name(base_tree_plot)

# Clean and format taxonomy table
taxonomy_lookup <- bind_rows(lapply(taxonomy_data, as_tibble), .id = "id") %>%
  filter(rank %in% c("species", "genus", "phylum")) %>%
  pivot_wider(names_from = rank, values_from = name) %>%
  mutate(species = factor(species, levels = species[order(species_order)]),
         species = as.character(species))

# Extract unique genus and phylum names
genus_order <- unique(taxonomy_lookup$genus)
phylum_order <- unique(taxonomy_lookup$phylum)
tree_data    <- base_tree_plot$data

# Compute label coordinates for genus and phylum
genus_annots  <- get_taxa_annot_positions(genus_order, "genus")
phylum_annots <- get_taxa_annot_positions(phylum_order, "phylum")

# Compute MRCA nodes for genus and phylum
highlight_data <- bind_rows(
  get_taxon_nodes("phylum") %>% mutate(group = phylum, rank = "phylum"),
  get_taxon_nodes("genus")  %>% mutate(group = genus,  rank = "genus"))

# Compose tree with labels and highlighted clades
tree_plot <- base_tree_plot +
  geom_label2(data = genus_annots,  aes(x = x, y = y, label = tax_label),
              fill = "white", size = 3, label.size = 0, family = "mono") +
  geom_label2(data = phylum_annots, aes(x = x, y = y, label = tax_label),
              fill = "white", size = 3, label.size = 0, family = "mono") +
  geom_hilight(data = highlight_data, aes(node = node, fill = group, alpha = rank),
               show.legend = FALSE) +
  scale_fill_viridis_d(option = "D") +
  scale_alpha_manual(values = c(phylum = 0.15, genus = 0.3)) +
  geom_tippoint(aes(color = factor(if_else(label %in% gt_species, "Target species", "Other species"))),
                size = 3, show.legend = FALSE) +
  scale_color_viridis_d(option = "D")

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

# Aggregate species-level differences to a given taxonomic rank (e.g., genus or phylum)
compute_taxon_differences <- function(taxonomy_lookup, species_differences, rank, order_levels) {
  taxonomy_lookup %>%
    inner_join(species_differences, by = "species") %>%     # Join with species-level difference values
    group_by(.data[[rank]], sample) %>%    # Group by rank and sample
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>% # Sum differences to get total deviation at that rank
    mutate(!!rank := factor(.data[[rank]], levels = rev(order_levels), ordered = TRUE))     # Order rank levels for plotting
}

# Generate a heatmap for the provided abundance difference data
plot_heatmap <- function(data, y_var, title, lower_bound, upper_bound, threshold) {
  ggplot(data, aes(x = sample, y = .data[[y_var]], fill = value)) +
    geom_tile(color = "grey80") +  # Add grid borders
    geom_text(aes(label = sprintf("%+.2f", value)), color = "white", size = 2) +
    scale_fill_viridis_c(option = "D", limits = c(lower_bound, upper_bound), oob = squish) + 
    scale_color_identity() +  # Maintain original colors
    labs(title = title, fill="Legend") + 
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

# Compute and reshape species-level abundance differences for plotting
species_differences <- read_counts %>%
  mutate(across(-species, ~ . / total_reads[[cur_column()]] * 100)) %>%   # Convert read counts to % abundance per sample
  inner_join(ground_truth, by = "species") %>%   # Join with ground truth abundances
  mutate(across(-c(species, abundance), ~ . - abundance)) %>%   # Compute % difference: observed - expected
  pivot_longer(-c(species, abundance), names_to = "sample", values_to = "value") %>%  # Reshape to long format: species, sample, value
  select(-abundance) %>%  # Remove expected abundance column
  mutate(species = factor(species, levels = rev(species_order), ordered = TRUE)) # Order species by tree tip order (reversed)
  #filter(!is.na(species))  %>%  # Drop rows with missing species names
  
# Compute genus and phylum-level differences by aggregating species-level values
genus_differences  <- compute_taxon_differences(taxonomy_lookup, species_differences, "genus", genus_order)
phylum_differences <- compute_taxon_differences(taxonomy_lookup, species_differences, "phylum", phylum_order)

# Determine plot bounds and thresholds for color scale
lower_bound <- min(species_differences$value, na.rm = TRUE)
upper_bound <- max(species_differences$value, na.rm = TRUE)
threshold <- abs(lower_bound+0.25*(upper_bound-lower_bound))

# Generate heatmaps for species, genus, and phylum levels
print(unique(species_differences$species))

species_heatmap <- plot_heatmap(species_differences, "species", "Species", lower_bound, upper_bound, threshold)
genus_heatmap <- plot_heatmap(genus_differences, "genus", "Genus", lower_bound, upper_bound, threshold)
phylum_heatmap <- plot_heatmap(phylum_differences, "phylum", "Phylum", lower_bound, upper_bound, threshold)

# ---------------------------------------
# FINAL COMPOSITE PLOT
# ---------------------------------------
# Combine phylogenetic tree and heatmaps into a single visualization
main_patch <- (tree_plot + species_heatmap) / (genus_heatmap + phylum_heatmap)

legend_hm <- get_legend(
  plot_heatmap(genus_differences, "genus", "", lower_bound, upper_bound, threshold) +
    theme_void() +
    theme(legend.title = element_text(size = 9, face = "bold"),
          legend.text = element_text(size = 8))
)

if (gt_flag) {
  base_tree_plot <- base_tree_plot +
    geom_tippoint(aes(color = factor(if_else(label %in% gt_species, "Target species", "Other species"))), size = 1) +
    scale_color_viridis_d(option = "D", name = "Category")
} else {
  base_tree_plot <- base_tree_plot + theme(legend.position = "none")
}

# Extract legend
legend_tree <- get_legend(base_tree_plot  +
                            theme(legend.title = element_text(size = 9, face = "bold"),
                                  legend.text = element_text(size = 8)))

# Overlay legends at top-left of full canvas
patchwork <- wrap_elements(full = main_patch) +
  inset_element(legend_hm, left = 0.01, bottom = 0.90, right = 0.05, top = 0.99, on_top = TRUE) +
  inset_element(legend_tree, left = 0.01, bottom = 0.80, right = 0.07, top = 0.89, on_top = TRUE) +
  plot_annotation(
    title = "Heatmaps of Taxa Abundance Across Samples",
    caption = dropped_text,
    theme = theme(
      plot.title = element_text(size = 15, face = "bold"),
      plot.caption = element_text(size = 8, face = "italic", hjust = 0)
    )
  )
# Save combined plot as an image
output_path <- file.path(dirname(opt$reports), "phylo_classification_comparison.png")
ggsave(output_path, plot = patchwork, width = 15, height = 15, dpi = 300)

