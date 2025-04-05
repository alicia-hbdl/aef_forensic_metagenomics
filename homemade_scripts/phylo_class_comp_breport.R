# Load required libraries
library(scales)
library(readr)       # Read tabular data
library(dplyr)       # Data manipulation
library(tidyr)       # Data tidying
library(stringr)     # String operations
library(ggplot2)     # Data visualization
library(reshape2)    # Data reshaping
library(tibble)      # Tidy data structures
library(taxize)      # Retrieve taxonomic IDs
library(ape)         # Phylogenetic analysis
library(treeio)      # Manipulate tree data
library(ggtree)      # Visualize phylogenetic trees
library(pheatmap)    # Generate heatmaps
library(patchwork)   # Combine multiple plots
library(ggtext)      # Required for element_markdown()

# Set NCBI API key for taxonomic queries
options(ENTREZ_KEY = "5a8133264ac32a3f11c0f1e666a90d96c908")

# Parse command-line argument (input file path)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript phylo_classification_comparison.R <path/to/pavian/comparison_table.tsv>")
file_path <- args[1]
if (!file.exists(file_path)) stop("Error: Comparison table not found!")

# Define the species that should be bold
highlight_species <- c("Bacillus subtilis", "Enterococcus faecalis", "Escherichia coli", 
                       "Limosilactobacillus fermentum", "Listeria monocytogenes", 
                       "Pseudomonas aeruginosa", "Salmonella enterica", 
                       "Staphylococcus aureus", "Saccharomyces cerevisiae", 
                       "Cryptococcus neoformans")

# Load and clean read count table
read_counts <- read_tsv(file_path, show_col_types = FALSE) %>%
  mutate(name = gsub("root><wbr>...><wbr>", "", name)) %>%  # Remove unwanted text from 'name'
  mutate(
    name = recode(name,
                  "Bacillus intestinalis" = "Bacillus spizizenii",  # Standardize taxonomic names (ensures consistency between the phylogenetic tree and species heatmap)
                  "Cryptococcus gattii VGI" = "Cryptococcus gattii",
                  "Cryptococcus gattii VGII" = "Cryptococcus deuterogattii")) %>%
  rename_with(~ gsub(".cladeReads", "", .x)) %>%  # Clean column names
  select(-c(taxRank, taxID, Max, lineage)) %>%  # Remove unnecessary columns
  { 
    existing <- .
    missing_species <- setdiff(highlight_species, existing$name)
    missing_rows <- tibble(name = missing_species) %>%
      mutate(across(-name, ~ 0))  # Fill other columns with 0
    bind_rows(existing, missing_rows)
  }

read_counts[is.na(read_counts)] <- 0  # Replace all NAs with 0

# Read CSV and transform using pipes
total_reads_path <- file.path(dirname(file_path), "../../../../processed_data/metagenomic/total_reads.csv")
total_reads <- read_csv(total_reads_path) %>%
  pivot_wider(names_from = sample, values_from = total_reads)  # Ensure correct column names

# Retrieve taxonomic hierarchy
tax_ids <- get_uid(read_counts$name)  # Get taxonomic IDs
taxonomy_data <- classification(tax_ids)  # Fetch taxonomy data

# Construct phylogenetic tree
tree <- class2tree(taxonomy_data)$phylo  # Convert classification to tree

# Plot phylogenetic tree
tree_plot <- ggtree(tree, layout = "rectangular", size = 0.3) +
  geom_tiplab(aes(label = label_pad(sub("^(\\w)\\w* (.+)$", "\\1. \\2", ifelse(is.na(label), "Unknown", as.character(label)))),
                  color = "black", 
                  fontface = ifelse(label %in% highlight_species, "bold", "plain")),
              align = TRUE, size = 4, family = 'mono') +  # Add tip labels with conditional formatting
  scale_color_identity() +  # Use specified text colors
  geom_label2(aes(label = label, subset = !isTip),  # Show all taxonomic levels
              hjust = 1.05, fill = "white", size = 3.5, label.size = 0, family = "mono") +
  theme_tree() +  # Apply tree theme
  theme(axis.ticks = element_blank()) +  # Remove axis ticks
  xlim(-mean(tree$edge.length), mean(tree$edge.length)^2/2)  # Adjust x-axis dynamically


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

# Define expected species composition
species_percentages <- tibble(species = read_counts$name) %>%
  mutate(
    expected_percentage = case_when(
      species %in% c("Bacillus subtilis", "Enterococcus faecalis", "Escherichia coli",
                     "Limosilactobacillus fermentum", "Listeria monocytogenes",
                     "Pseudomonas aeruginosa", "Salmonella enterica", "Staphylococcus aureus") ~ 12,  # Bacterial composition
      species %in% c("Saccharomyces cerevisiae", "Cryptococcus neoformans") ~ 2,  # Yeast composition
      TRUE ~ 0))  # Default for others

# Aggregate expected percentages at genus level
genus_percentages <- taxonomy_lookup %>%
  inner_join(species_percentages, by = "species") %>%
  group_by(genus) %>%
  summarise(expected_percentage = sum(expected_percentage, na.rm = TRUE))

# Aggregate expected percentages at phylum level
phylum_percentages <- taxonomy_lookup %>%
  inner_join(species_percentages, by = "species") %>%
  group_by(phylum) %>%
  summarise(expected_percentage = sum(expected_percentage, na.rm = TRUE))

# Compute species-level differences between observed and expected abundance
species_differences <- read_counts %>%
  mutate(across(-name, ~ . / total_reads[[cur_column()]] * 100)) %>% # Convert counts to percentages
  inner_join(species_percentages, by = c("name" = "species")) %>%  # Match species
  rename(species = name) %>%
  mutate(across(-c(species, expected_percentage), ~ (.-expected_percentage)))  %>%  # Compute deviation
  select(-expected_percentage) %>%
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
    geom_text(aes(label = sprintf("%+.2f", value)),
              color = ifelse(is.na(data$value), "transparent", 
                             ifelse(abs(data$value) > threshold, "white", "black")),
              size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                         limits = c(lower_bound, upper_bound), oob = squish) + # Set fixed bounds
    scale_color_identity() +  # Maintain original colors
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, family = "mono", face = "bold"),  
      axis.title.x = element_blank(),
      axis.text.y = if (y_var == "species") element_blank() else element_text(family = "mono", face = "bold"),
      axis.title.y = element_blank(),
      legend.position = if (y_var != "genus") "none" else "left",  
      legend.title = element_text(size = 10)) +
    labs(title = title, fill = "% Difference")
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
layout <- "
AAAAABB
#CCCDDD
"
patchwork <- tree_plot + species_heatmap + genus_heatmap + phylum_heatmap + plot_layout(design = layout) +
  plot_annotation(
    title = "Phylogenetic Heatmap of Taxa Abundance Changes Across Samples",
    caption = "**ZC samples**: Derived from the Zymo Microbial Community Standard, a mixture of live microbial cells used for benchmarking DNA extraction methods.  \n**ZP samples**: Obtained from the Zymo Microbial Community DNA Standard, containing pre-extracted purified DNA from the same microbial community to eliminate DNA extraction variability.",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 19, face = "bold"),
      plot.caption = element_markdown(hjust = 0, size = 12, lineheight=1.2)))

# Save combined plot as an image
output_path <- file.path(dirname(file_path), "phylo_classification_comparison.png")
ggsave(output_path, plot = patchwork, width = 20, height = 16, dpi = 300)
