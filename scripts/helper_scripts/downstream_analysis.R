#!/usr/bin/env Rscript

# Combine Bracken reports and associated metadata into a SummarizedExperiment for downstream analysis of metagenomic pipeline configurations.

# Usage: Rscript downstream_analysis.R -t <path/to/ground_truth.csv> -s <path/to/runs_summary.csv> breport1.csv breport2.csv ...

suppressPackageStartupMessages({
  library(optparse)             # For command-line argument parsing
  library(tidyverse)            # For data manipulation
  library(SummarizedExperiment) # To create the SE object
  library(S4Vectors)            # For SE metadata handling
  library(patchwork)
  library(matrixStats)
  library(Rtsne)
  library(umap)
  library(reshape2)  # for melt()
  library(ggplotify)
  library(pheatmap)
  library(viridis)
  library(ggrepel)
})

# ----------------------------- CLI Parsing ------------------------------

# Define command-line options for ground truth and run metadata
option_list <- list(
  make_option(c("-t", "--ground-truth"), type = "character", help = "Path to ground truth file"),
  make_option(c("-s", "--runs-summary"), type = "character", help = "Path to runs summary CSV")
)

# Parse options and positional arguments (Bracken report paths)
opt <- parse_args(
  OptionParser(option_list = option_list,
               usage = "Usage: %prog -t ground_truth.csv -s runs_summary.csv breport1.csv breport2.csv ..."),
  positional_arguments = TRUE
)

# Extract named and positional arguments
opts <- opt$options
breport_paths <- opt$args

# Validate required arguments
if (length(breport_paths) == 0 || is.null(opts$`ground-truth`) || is.null(opts$`runs-summary`)) {
  print_help(OptionParser(option_list = option_list))
  stop("❌ Must provide at least one Bracken report and both --ground-truth and --runs-summary.", call. = FALSE)
}

# Check all specified files exist
all_paths <- c(breport_paths, opts$`ground-truth`, opts$`runs-summary`)
invisible(lapply(all_paths, function(f) {
  if (!file.exists(f)) stop(paste("❌ File not found:", f))
}))

# ----------------------------- Read Reports -----------------------------

# Initialize an empty table with only species column
combined_table <- data.frame(species = character(), stringsAsFactors = FALSE)

# Iterate over each Bracken report to extract and merge species-level read counts
for (file_path in breport_paths) {
  run_id <- basename(dirname(file_path)) # Use parent directory name as run ID

  # Read the Bracken report and calculate mean read counts across numeric columns
  run_species <- read_csv(file_path, show_col_types = FALSE) %>%
    mutate(!!run_id := rowMeans(across(where(is.numeric), ~ replace_na(.x, 0)))) %>%
    select(species, !!sym(run_id)) # Keep only species and computed column
  
  # Merge with previous samples by species
  combined_table <- if (nrow(combined_table) == 0) {
    run_species
  } else {
    full_join(combined_table, run_species, by = "species") %>%
      mutate(across(where(is.numeric), ~ replace_na(.x, 0)))
  }
}

# Add ground truth to combined table
#truth_data <- read_csv(opts$`ground-truth`, show_col_types = FALSE) #%>%
  #mutate(abundance = abundance)  

#truth_data$ground_truth <- truth_data$abundance
#truth_data$abundance <- NULL

# Add ground truth and replace missing values with 0 for matrix compatibility
#combined_table <- left_join(combined_table, truth_data, by = "species") %>%
  #mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

# Convert combined table to a numeric matrix with species as row names
count_matrix <- combined_table %>%
  column_to_rownames("species") %>%
  as.matrix()

colnames(count_matrix)
# -------------------------- Metadata Assembly ---------------------------

# Create basic row metadata with species names
row_data <- DataFrame(species = rownames(count_matrix))

# Read and aggregate run-level metadata (e.g. runtime, trimming parameters)
col_data <- read_csv(opts$`runs-summary`, show_col_types = FALSE) %>%
  group_by(run_id) %>%
  summarise(
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)) , # Average numeric values
    across(where(~ !is.numeric(.x)), dplyr::first), # Keep first instance for non-numeric
    .groups = "drop"
  ) %>%
  column_to_rownames("run_id") %>%
  DataFrame()

# Create an empty row with the same column structure as col_data
#empty_row <- as.data.frame(as.list(rep(NA, ncol(col_data))))
#colnames(empty_row) <- colnames(col_data)
#rownames(empty_row) <- "ground_truth"

# Add it
#col_data <- rbind(col_data, empty_row)

# Align sample order
ordered_samples <- sort(colnames(count_matrix))
count_matrix <- count_matrix[, ordered_samples]
col_data <- col_data[ordered_samples, , drop = FALSE]

# Create SummarizedExperiment object
se <- SummarizedExperiment(
  assays = list(counts = count_matrix),
  rowData = row_data,
  colData = col_data
)

# -------------------------- Filter Count Data ---------------------------

# Use raw counts assay for filtering
mat <- assay(se, "counts")

# Filter species (rows) with ≥20 reads in ≥50% of samples
min_reads <- 20
min_samples <- ncol(se) / 2
rows_to_keep <- rowSums(mat >= min_reads) >= min_samples
se_filtered <- se[rows_to_keep, ]

# Filter samples (columns) with at least 100 total reads
filtered_mat <- assay(se_filtered, "counts")
cols_to_keep <- colSums(filtered_mat) >= 100
se_filtered <- se_filtered[, cols_to_keep]

# -------------------------- Visualize Library Sizes ---------------------------

# Extract filtered count matrix
filtered_mat <- assay(se_filtered, "counts")

# Compute total classified reads per pipeline run (sample)
lib_sizes <- colSums(filtered_mat)

# Create data frame for plotting
df_lib <- data.frame(Run = names(lib_sizes), Reads = lib_sizes)

# Plot: classified read counts per pipeline run
p1 <- ggplot(df_lib, aes(x = Run, y = Reads, fill = Run)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  labs(x = "Pipeline Run", y = "Classified Reads", title = "Classified Reads per Pipeline Run") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# -------------------------- Normalize Count Data ---------------------------

# Correlation: species abundance vs. library size (raw counts)
cor_vals_raw <- apply(filtered_mat, 1, \(x) cor(x, lib_sizes))
df_cor_raw <- data.frame(Correlation = cor_vals_raw)
p2 <- ggplot(df_cor_raw, aes(x = Correlation)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(cor_vals_raw), sd = sd(cor_vals_raw)), color = "red", linewidth = 1) +
  labs(x = "Correlation Coefficient", y = "Density", title = "Raw Abundance vs. Library Size") +
  theme_minimal()

# Normalize to CPM
filtered_cpm <- t(t(filtered_mat) / lib_sizes * 1e6)
print(filtered_cpm)

# Correlation after CPM normalization
cor_vals_cpm <- apply(filtered_cpm, 1, \(x) cor(x, lib_sizes))
df_cor_cpm <- data.frame(Correlation = cor_vals_cpm)
p3 <- ggplot(df_cor_cpm, aes(x = Correlation)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(cor_vals_cpm), sd = sd(cor_vals_cpm)), color = "red", linewidth = 1) +
  labs(x = "Correlation (CPM)", y = "Number of Species", title = "CPM-Normalized Correlation") +
  theme_minimal()

# -------------------------- Stabilized Normalized Data ---------------------------

# Log mean vs. log standard deviation (CPM scale)
df_var_cpm <- data.frame(Mean = log1p(rowMeans(filtered_cpm)), SD = log1p(rowSds(filtered_cpm)))
cor_cpm <- cor(df_var_cpm$Mean, df_var_cpm$SD)
p4 <- ggplot(df_var_cpm, aes(x = Mean, y = SD)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
  annotate("text", x = min(df_var_cpm$Mean), y = max(df_var_cpm$SD),
           label = paste0("r = ", round(cor_cpm, 3)), hjust = 0, size = 4, fontface = "italic") +
  labs(x = "log(Mean CPM + 1)", y = "log(SD CPM + 1)", title = "Species Variability (CPM)") +
  theme_minimal()

# Log2 transform for downstream analysis
filtered_logcpm <- log2(filtered_cpm+1)

# Log mean vs. log SD (logCPM scale)
df_var_logcpm <- data.frame(Mean = log1p(rowMeans(filtered_logcpm)), SD = log1p(rowSds(filtered_logcpm)))
cor_logcpm <- cor(df_var_logcpm$Mean, df_var_logcpm$SD)
p5 <- ggplot(df_var_logcpm, aes(x = Mean, y = SD)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
  annotate("text", x = max(df_var_logcpm$Mean)-0.2, y = max(df_var_logcpm$SD),
           label = paste0("r = ", round(cor_logcpm, 3)), hjust = 0, size = 4, fontface = "italic") +
  labs(x = "log(Mean logCPM + 1)", y = "log(SD logCPM + 1)", title = "Species Variability (logCPM)") +
  theme_minimal()

combined_transformations <- p1 / (p2 + p3) / (p4 + p5) 
ggsave("combined_transformations.png", combined_transformations, width = 10, height = 10) 

# -------------------------- Clustering ---------------------------

# Use PCA, t-SNE, and UMAP to reduce species-space vectors into 2D.
#  Plot the reduced 2D coordinates of each pipeline setting.
#  Annotate points with their setting name (including “truth”).

# PCA: visualize configurations in reduced space
pca <- prcomp(t(filtered_logcpm), scale. = TRUE)
pca_df <- as.data.frame(pca$x) %>% rownames_to_column("run_id")
p6 <- ggplot(pca_df, aes(x = PC1, y = PC2, label = sample, color = sample)) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = Inf) +
  labs(title = "PCA", x = "PC1", y = "PC2") +
  theme_minimal()

# t-SNE: visualize nonlinear separation
tsne <- Rtsne(t(filtered_logcpm), perplexity = 1)
tsne_df <- as.data.frame(tsne$Y)
colnames(tsne_df) <- c("tsne_1", "tsne_2")
tsne_df$run_id <- colnames(filtered_logcpm)
p7 <- ggplot(tsne_df, aes(x = tsne_1, y = tsne_2, label = sample, color = run_id)) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = Inf) +
  labs(title = "t-SNE", x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()

# UMAP: another non-linear projection
umap_cfg <- umap.defaults; umap_cfg$n_neighbors <- 3
umap_res <- umap(t(filtered_logcpm), config = umap_cfg)
umap_df <- as.data.frame(umap_res$layout)
colnames(umap_df) <- c("umap_1", "umap_2")
umap_df$run_id <- colnames(filtered_logcpm)
p8 <- ggplot(umap_df, aes(x = umap_1, y = umap_2, label = sample, color = run_id)) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = Inf) +
  labs(title = "UMAP", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()
# Helper: wrap pheatmap as a ggplot-compatible grob
pheatmap_grob <- function(mat) {
  p <- pheatmap(
    mat,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    color = viridis(100, option = "D", direction = -1),
    angle_col = 45,  # <<<<< Rotate bottom labels
    silent = TRUE
  )
  as.ggplot(p[[4]])  # Extract heatmap grob and wrap in ggplot-compatible form
}

# Distances in original logCPM space (species abundance profiles)
raw_dists   <- as.matrix(dist(t(filtered_logcpm)))
p9  <- pheatmap_grob(raw_dists)

# Distances in PCA space (first 2 principal components)
pca_coords <- pca_df[, c("PC1", "PC2")]
rownames(pca_coords) <- pca_df$sample
pca_dists <- as.matrix(dist(pca_coords))
pca_dists <- pca_dists / max(pca_dists)
p10 <- pheatmap_grob(pca_dists)

# Distances in UMAP space
rownames(umap_df) <- umap_df$sample
umap_dists  <- as.matrix(dist(umap_df[, c("umap_1", "umap_2")]))
umap_dists <- umap_dists / max(umap_dists)
p11 <- pheatmap_grob(umap_dists)

# Distances in t-SNE space
rownames(tsne_df) <- tsne_df$sample
tsne_dists  <- as.matrix(dist(tsne_df[, c("tsne_1", "tsne_2")]))
tsne_dists <- tsne_dists / max(tsne_dists)
p12 <- pheatmap_grob(tsne_dists)

combined_distances <- (p6 + p7 + p8) / (p10 + p11 + p12)
ggsave("combined_distances.png", combined_distances, width = 12, height = 8) 

# TODO: Use distances from "truth" column to rank best-matching configurations

