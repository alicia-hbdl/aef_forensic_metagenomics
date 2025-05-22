#!/usr/bin/env Rscript

# This script aggregates Bracken output and associated metadata into a SummarizedExperiment object,
# and performs downstream analysis of metagenomic pipeline performance. It generates:
# - Read retention and read length progression plots across pipeline stages,
# - Clustering visualizations of species profiles (PCA, t-SNE, UMAP),
# - Comparative analyses of classification accuracy against ground truth (L2, AUPR, precision/recall).

# Usage: Rscript downstream_analysis.R -t <path/to/ground_truth.csv> -s <path/to/runs_summary.csv> breport1.csv breport2.csv ...

#=========================================================
# Environment Setup
#=========================================================

## ------------------ Load Required Libraries ------------------

suppressPackageStartupMessages({
  library(optparse)             # CLI argument parsing
  library(tidyverse)            # Tidy data manipulation and plotting
  library(SummarizedExperiment) # SummarizedExperiment object 
  library(S4Vectors)            # Metadata handling for SE objects
  library(patchwork)            # Combine multiple ggplots
  library(matrixStats)          # Row/column statistics for matrices
  library(Rtsne)                # t-SNE dimensionality reduction
  library(umap)                 # UMAP dimensionality reduction
  library(reshape2)             # Data reshaping (e.g., melt)
  library(ggplotify)            # Convert plots/tables to ggplot objects
  library(pheatmap)             # Clustered heatmaps
  library(viridis)              # Colorblind-friendly color scales
  library(ggrepel)              # Non-overlapping text labels
  library(ggpubr)               # Statistical tests and plot annotation
  library(scales)
})

## ------------------ Parse Command-Line Arguments ------------------

# Define CLI options for ground truth and run summary files
option_list <- list(
  make_option(c("-t", "--ground-truth"), type = "character", help = "Path to ground truth file"),
  make_option(c("-s", "--runs-summary"), type = "character", help = "Path to runs summary CSV")
)

# Parse CLI options and positional arguments (Bracken report paths)
opt <- parse_args(
  OptionParser(option_list = option_list,
               usage = "Usage: %prog -t ground_truth.csv -s runs_summary.csv breport1.csv breport2.csv ..."),
  positional_arguments = TRUE
)

# Extract parsed arguments
opts <- opt$options
breport_paths <- opt$args

# Validate required inputs
if (length(breport_paths) == 0 || is.null(opts$`ground-truth`) || is.null(opts$`runs-summary`)) {
  print_help(OptionParser(option_list = option_list))
  stop("❌ Must provide at least one Bracken report and both --ground-truth and --runs-summary.", call. = FALSE)
}

# Check that all input files exist
all_paths <- c(breport_paths, opts$`ground-truth`, opts$`runs-summary`)
invisible(lapply(all_paths, function(f) {
  if (!file.exists(f)) stop(paste("❌ File not found:", f))
}))

#=========================================================
# Creating Summarized Experiment 
#=========================================================

## ------------- Load Bracken reports ------------- 

# Initialize empty table to collect species counts
combined_table <- data.frame(species = character(), stringsAsFactors = FALSE)

# Loop through each Bracken report and merge species-level counts
for (file_path in breport_paths) {
  run_id <- basename(dirname(file_path)) # Extract run ID from folder name

  # Read report and compute mean across numeric columns
  run_species <- read_csv(file_path, show_col_types = FALSE) %>%
    mutate(!!run_id := rowMeans(across(where(is.numeric), ~ replace_na(.x, 0)))) %>% # Replace NA with 0 
    select(species, !!sym(run_id)) # Keep only species and new column

  # Merge into the combined table
  combined_table <- if (nrow(combined_table) == 0) {
    run_species
  } else {
    full_join(combined_table, run_species, by = "species") %>%
      mutate(across(where(is.numeric), ~ replace_na(.x, 0))) # Fill missing values with 0
  }
}

# Convert combined table to matrix format with species as rownames
count_matrix <- combined_table %>%
  column_to_rownames("species") %>%
  as.matrix()

## ------------- Assemble metadata ------------- 

# Create row metadata for species
row_data <- DataFrame(species = rownames(count_matrix))

# Read and summarize run metadata (aggregates across all samples in a run)
col_data <- read_csv(opts$`runs-summary`, show_col_types = FALSE) %>%
  group_by(run_id) %>% 
  summarise(
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), # Average numeric values
    across(where(~ !is.numeric(.x)), dplyr::first), # Keep first instance for non-numeric
    .groups = "drop"
  ) %>%
  column_to_rownames("run_id") %>%
  DataFrame()

# Align sample order between count matrix and metadata
ordered_samples <- sort(colnames(count_matrix))
count_matrix <- count_matrix[, ordered_samples]
col_data <- col_data[ordered_samples, , drop = FALSE]

# Create SummarizedExperiment object
se <- SummarizedExperiment(
  assays = list(counts = count_matrix),
  rowData = row_data,
  colData = col_data
)

#=========================================================
# Read Length and Size Progression Analysis 
#=========================================================

## ------------------ Read Retention ------------------

# Define pipeline stages (read counts after major steps)
stage_names <- c("trim_paired", "bt_paired", "kraken2_total", "kraken2_classified", "bracken_total")

# Fallback values if QC/trimming skipped (causes missing values)
mean_trim <- mean(col_data$trim_paired[!is.nan(col_data$trim_paired)], na.rm = TRUE)
mean_bt   <- mean(col_data$bt_paired[!is.nan(col_data$bt_paired)], na.rm = TRUE)

# Aggregate read counts by database (species presence/absence in DB affects classification)
aggregated_summary <- as_tibble(col_data, rownames = "run_id") %>%
  select(db_name, all_of(stage_names)) %>%
  group_by(db_name) %>%                                  
  mutate( 
    trim_paired = if_else(is.nan(trim_paired), mean_trim, trim_paired), # Impute missing
    bt_paired   = if_else(is.nan(bt_paired), mean_bt, bt_paired)
  ) %>%
  summarise(across(all_of(stage_names), ~ mean(.x, na.rm = TRUE)), .groups = "drop") # Equal weight per DB
                
# Reshape & normalize to % of raw reads
long_summary <- aggregated_summary %>%
  pivot_longer(cols = all_of(stage_names), names_to = "Original", values_to = "Reads") %>% 
  mutate(Stage = factor(Original, levels = stage_names)) %>% # Preserve stage order
  group_by(db_name = db_name) %>%
  mutate(Fraction = Reads / Reads[Original == "trim_paired"]) %>% # Normalize to % of raw reads
  mutate(Stage = recode(Stage, # Rename for clarity
                        "trim_paired"        = "Raw",
                        "bt_paired"          = "Trimmed",
                        "kraken2_total"      = "Host Filtered",
                        "kraken2_classified" = "Classified (Kraken)",
                        "bracken_total"      = "Classified (Bracken)"
  )) %>%
  ungroup()

# Mean % retained at each stage across all DBs
mean_summary <- long_summary %>%
  group_by(Stage) %>%
  summarise(Fraction = mean(Fraction, na.rm = TRUE), .groups = "drop")

# Plot read retention progression across pipeline stages
read_retention <- ggplot(long_summary, aes(x = Stage, y = Fraction, group = db_name, color = db_name)) +
  geom_line(linewidth = 0.5, alpha = 0.6) + # Per-database trend
  geom_line(data = mean_summary, aes(x = Stage, y = Fraction, group = 1), 
            color = "grey", linewidth = 1, inherit.aes = FALSE) + # Mean trend 
  geom_point(data = mean_summary, aes(x = Stage, y = Fraction), 
             color = "grey", size = 2, inherit.aes = FALSE) + # Mean values 
  geom_text_repel(data = mean_summary, aes(x = Stage, y = Fraction, label = percent(Fraction, accuracy = 0.01)), 
                  color = "black", size = 3, inherit.aes = FALSE) + # Add % labels to mean
  scale_color_viridis_d(option = "D") +
  labs(title = "Progression of Read Count", y = "Proportion of Raw Reads", x = NULL, color = "Database") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

## ------------------ Read Length Progression ------------------

# Define pipeline stages (read counts after major steps)
bp_cols <- c("preqc_avg_len_r1", "preqc_avg_len_r2", "preqc_med_len_r1", "preqc_med_len_r2",
  "postqc_avg_len_r1", "postqc_avg_len_r2", "postqc_med_len_r1", "postqc_med_len_r2",
  "bt_avg_r1", "bt_avg_r2", "bt_med_r1", "bt_med_r2",
  "clas_avg_r1", "clas_avg_r2", "clas_med_r1", "clas_med_r2",
  "unclas_avg_r1", "unclas_avg_r2", "unclas_med_r1", "unclas_med_r2")

# Calculate per-run min read length per stage (across R1/R2) to account for asymmetry
min_read_data <- as_tibble(col_data, rownames = "run_id") %>%
  select(db_name, all_of(bp_cols)) %>%  
  mutate(
    preqc_avg  = pmin(preqc_avg_len_r1, preqc_avg_len_r2, na.rm = TRUE),
    preqc_med  = pmin(preqc_med_len_r1, preqc_med_len_r2, na.rm = TRUE),
    postqc_avg = pmin(postqc_avg_len_r1, postqc_avg_len_r2, na.rm = TRUE),
    postqc_med = pmin(postqc_med_len_r1, postqc_med_len_r2, na.rm = TRUE),
    bt_avg     = pmin(bt_avg_r1, bt_avg_r2, na.rm = TRUE),
    bt_med     = pmin(bt_med_r1, bt_med_r2, na.rm = TRUE),
    clas_avg   = pmin(clas_avg_r1, clas_avg_r2, na.rm = TRUE),
    clas_med   = pmin(clas_med_r1, clas_med_r2, na.rm = TRUE), 
    unclas_avg   = pmin(unclas_avg_r1, unclas_avg_r2, na.rm = TRUE),
    unclas_med   = pmin(unclas_med_r1, unclas_med_r2, na.rm = TRUE)
  ) 

# Fallback values if QC/trimming skipped (causes missing values)
mean_preqc_avg <- mean(min_read_data$preqc_avg[!is.nan(min_read_data$preqc_avg)], na.rm = TRUE)
mean_preqc_med  <- mean(min_read_data$preqc_med[!is.nan(min_read_data$preqc_med)], na.rm = TRUE)
mean_postqc_avg <- mean(min_read_data$postqc_avg[!is.nan(min_read_data$postqc_avg)], na.rm = TRUE)
mean_postqc_med  <- mean(min_read_data$postqc_med[!is.nan(min_read_data$postqc_med)], na.rm = TRUE)

# Impute NAs and compute average read lengths per stage grouped by database
length_summary <- min_read_data %>%
  group_by(db_name) %>%
  mutate(
    preqc_avg = if_else(is.nan(preqc_avg), mean_preqc_avg, preqc_avg),
    preqc_med = if_else(is.nan(preqc_med), mean_preqc_med, preqc_med),
    postqc_avg = if_else(is.nan(postqc_avg), mean_postqc_avg, postqc_avg),
    postqc_med = if_else(is.nan(postqc_med), mean_postqc_med, postqc_med)
  ) %>%
  summarise(across(c(preqc_avg, preqc_med, postqc_avg, postqc_med, bt_avg, bt_med, clas_avg, clas_med, unclas_avg, unclas_med),
                   ~ mean(.x, na.rm = TRUE)),.groups = "drop")

# Reshape summary to long format for plotting 
length_long <- length_summary %>%
  pivot_longer(cols = -c(db_name, unclas_avg, unclas_med), names_to = "Stage_Metric", values_to = "Length") %>%
  separate(Stage_Metric, into = c("Stage", "Metric"), sep = "_") %>%
  mutate(Stage  = factor(Stage, levels = c("preqc", "postqc", "bt", "clas"), labels = c("Raw", "Trimmed", "Host Filtered", "Classified")),
         Metric = if_else(Metric == "avg", "Mean", "Median"))

# Compute mean read length per stage and metric 
mean_length_summary <- length_long %>%
  group_by(Stage, Metric) %>%
  summarise(Length = mean(Length, na.rm = TRUE), .groups = "drop")

# Plot progression of read length across pipeline stages
read_length <- ggplot(length_long, aes(x = Stage, y = Length, group = interaction(db_name, Metric), color = db_name, linetype = Metric)) +
  geom_line(linewidth = 0.6, alpha = 0.6, na.rm = TRUE) +
  geom_line(data = mean_length_summary, aes(x = Stage, y = Length, group = Metric), color = "grey", linewidth = 1) +
  geom_point(data = mean_length_summary, aes(x = Stage, y = Length), color = "grey", size = 2, inherit.aes = FALSE) +
  geom_text_repel(data = mean_length_summary, aes(x = Stage, y = Length, label = round(Length, 1)), color = "black", size = 3, inherit.aes = FALSE) +
  labs(title = "Progression of Read Length", y = "Read Length (bp)", x = NULL, linetype = "Metric", color = "Database") +
  scale_linetype_manual(values = c("Median" = "dashed", "Mean" = "solid")) +
  scale_color_viridis_d(option = "D") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8))

## ------------------ Boxplot: Classified vs Unclassified ------------------

# Reshape median classified/unclassified read lengths into long format for plotting
clas_unclas <- length_summary %>%
  select(db_name, clas_med, unclas_med) %>%  # Keep relevant columns
  pivot_longer(cols = c(clas_med, unclas_med), names_to = "Type", values_to = "Reads") %>%
  mutate(Type = recode(Type, "clas_med" = "Classified", "unclas_med" = "Unclassified"))  # Rename for clarity

# Boxplot comparing median read lengths of classified vs unclassified reads
# Unpaired Wilcoxon test: non-parametric (handles skewed data), compares medians of independent groups (classified vs unclassified reads)
boxplot <- ggplot(clas_unclas, aes(x = Type, y = Reads, fill = Type)) +
  geom_boxplot(alpha = 0.8, color = "gray", size = 0.5) + 
  geom_point(size = 1, color = "gray") + # Add raw data points
  stat_compare_means(comparisons = list(c("Classified", "Unclassified")), method = "wilcox.test", size = 3) + # Unpaired Wilcoxon test
  stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)), size = 3, color = "black") + # Show median values
  scale_fill_viridis_d(option = "D") + 
  labs(title = "Median Read Length Comparison", x = NULL, y = "Read Length (bp)") +
  theme_bw() +                                                
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Combine all plots (boxplot, read retention, and read length progression) and save
read_progression <- boxplot + read_retention + read_length + 
  plot_annotation(caption = "Note: The boxplot shows the median, and the line plot the mean, of median read lengths per database, hence the slight difference.",
                  theme = theme(plot.caption = element_text(size = 8, hjust = 0, face = "italic")))
            
ggsave("read_progression.png", read_progression, width = 12, height = 4, dpi = 300)


#=========================================================
# Clustering with Ground Truth
#=========================================================

##------------- Filter Count Matrix -------------

se_filtered <- se %>%
  subset(rowSums(assay(., "counts") >= 20) >= (ncol(.) / 2)) %>% # Keep species with ≥20 reads in ≥50% of samples
  subset(, colSums(assay(., "counts")) >= 100) # Keep samples with ≥100 total reads

# Extract raw counts assay
filtered_mat <- assay(se_filtered, "counts")

# Compute total reads per sample (library size)
lib_sizes <- colSums(filtered_mat)

##------------- Visualize Library Sizes -------------

# Wrap into a dataframe for plotting
df_lib <- data.frame(Run = names(lib_sizes), Reads = lib_sizes)

# Barplot: total classified reads per run
p1 <- ggplot(
  data.frame(Run = names(lib_sizes), Reads = lib_sizes), 
  aes(x = Run, y = Reads, fill = Run)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  scale_fill_viridis_d(option = "D") +  
  labs(title = "Classified Reads per Pipeline Run", x = "Pipeline Run", y = "Classified Reads") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1))

##------------- Normalize Count Matrix -------------

# Correlation: how strongly does each species' abundance correlate with library size?

# Raw correlations: species abundance vs. library size (shows library size bias)
cor_vals_raw <- apply(filtered_mat, 1, \(x) cor(x, lib_sizes))
p2 <- tibble(Correlation = cor_vals_raw) %>%
  ggplot(aes(x = Correlation)) +
  geom_histogram(binwidth = 0.05, fill = viridis(1, option = "D"), color = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(cor_vals_raw), sd = sd(cor_vals_raw)), color = "grey", linewidth = 1) + # Normal distribution curve 
  labs(title = "Species Abundance Correlation with Library Size (Raw)", x = "Correlation", y = "Density") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1))

# CPM normalization (adjusts for library size; Zymo is balanced → no genome correction)
assay(se_filtered, "cpm") <- t(t(filtered_mat) / lib_sizes * 1e6)  # Add as new assay

# Correlation after CPM normalization (bias should be reduced)
cor_vals_cpm <- apply(assay(se_filtered, "cpm"), 1, \(x) cor(x, lib_sizes))
p3 <- tibble(Correlation = cor_vals_cpm) %>%
  ggplot(aes(x = Correlation)) +
  geom_histogram(binwidth = 0.05, fill = viridis(1, option = "D"), color = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(cor_vals_cpm), sd = sd(cor_vals_cpm)), color = "grey", linewidth = 1) + # Normal distribution curve 
  labs(title = "Species Abundance Correlation with Library Size (CPM)", x = "Correlation", y = "Density") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1))

## ------------- Log Transformation -------------

# Plot log(Mean) vs log(SD) to examine heteroscedasticity; log(+1) avoids issues with 0s
df_var_cpm <- data.frame(Mean = log1p(rowMeans(assay(se_filtered, "cpm"))), SD = log1p(rowSds(assay(se_filtered, "cpm"))))
p4 <- ggplot(df_var_cpm, aes(x = Mean, y = SD, color = viridis(1, option = "D"))) +
  geom_point(size = 2, color = viridis(1, option = "D")) +  
  geom_smooth(formula = 'y ~ x', method = "lm", se = FALSE, color = "grey", linewidth = 1) +
  annotate("text", x = min(df_var_cpm$Mean), y = max(df_var_cpm$SD),
           label = paste0("r = ", round(cor(df_var_cpm$Mean, df_var_cpm$SD), 3)), hjust = 0, size = 3, fontface = "italic") +
  labs(x = "log(Mean CPM + 1)", y = "log(SD CPM + 1)", title = "Species Abundance Variability (CPM)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Variance stabilization: log2 transform CPM 
assay(se_filtered, "logcpm") <- log2(assay(se_filtered, "cpm")+1)

# Analyze mean–SD relationship again
df_var_logcpm <- data.frame(Mean = rowMeans(assay(se_filtered, "logcpm")), SD = rowSds(assay(se_filtered, "logcpm")))
p5 <- ggplot(df_var_logcpm, aes(x = Mean, y = SD)) +
  geom_point(size = 2, color = viridis(1, option = "D")) +  
  geom_smooth(formula = 'y ~ x', method = "lm", se = FALSE, color = "grey", linewidth = 1) +
  annotate("text", x = max(df_var_logcpm$Mean) - 0.2, y = max(df_var_logcpm$SD),
           label = paste0("r = ", round(cor(df_var_logcpm$Mean, df_var_logcpm$SD), 3)), hjust = +1, size = 3, fontface = "italic") +
  labs(x = "Mean logCPM", y = "SD logCPM", title = "Species Abundance Variability (logCPM)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Combine and export figure
combined_transformations <- p1 / (p2 + p3) / (p4 + p5) 
ggsave("combined_transformations.png", combined_transformations, width = 10, height = 10) 

#=========================================================
# Clustering with Ground Truth
#=========================================================

##------------- Dimensionality Reduction -------------

# PCA
pca <- prcomp(t(assay(se_filtered, "logcpm") ), scale. = TRUE)
pca_df <- as.data.frame(pca$x) %>% rownames_to_column("run_id")
p6 <- ggplot(pca_df, aes(x = PC1, y = PC2, label = sample, color = sample)) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = Inf) +
  labs(title = "PCA", x = "PC1", y = "PC2") +
  theme_minimal()

# t-SNE
tsne <- Rtsne(t(assay(se_filtered, "logcpm") ), perplexity = 3)
tsne_df <- as.data.frame(tsne$Y)
colnames(tsne_df) <- c("tsne_1", "tsne_2")
tsne_df$run_id <- colnames(assay(se_filtered, "logcpm") )
p7 <- ggplot(tsne_df, aes(x = tsne_1, y = tsne_2, label = sample, color = run_id)) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = Inf) +
  labs(title = "t-SNE", x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()

# UMAP
umap_cfg <- umap.defaults; umap_cfg$n_neighbors <- 3
umap_res <- umap(t(assay(se_filtered, "logcpm") ), config = umap_cfg)
umap_df <- as.data.frame(umap_res$layout)
colnames(umap_df) <- c("umap_1", "umap_2")
umap_df$run_id <- colnames(assay(se_filtered, "logcpm") )
p8 <- ggplot(umap_df, aes(x = umap_1, y = umap_2, label = sample, color = run_id)) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = Inf) +
  labs(title = "UMAP", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()

##------------- Distance Computation -------------

# Helper: wrap pheatmap as a ggplot-compatible grob
pheatmap_grob <- function(mat) {
  p <- pheatmap(
    mat,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    color = viridis(100, option = "D", direction = -1),
    angle_col = 45,  # Rotate bottom labels
    silent = TRUE
  )
  as.ggplot(p[[4]])  # Extract heatmap grob and wrap in ggplot-compatible form
}

# Distances in original logCPM space (species abundance profiles)
raw_dists   <- as.matrix(dist(t(assay(se_filtered, "logcpm") )))
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

#=========================================================
# Precision and Recall Graphs 
#=========================================================

# Warning message:
#   One or more parsing issues, call `problems()` on your data frame for details, e.g.:
#   dat <- vroom(...)
# problems(dat) 
# Warning message:
#   There were 5 warnings in `summarise()`.
# The first warning was:
#   ℹ In argument: `across(...)`.
# ℹ In group 1: `db_name = "k2_housepets_250510"`.
# Caused by warning in `mean.default()`:
#   ! argument is not numeric or logical: returning NA
# ℹ Run `dplyr::last_dplyr_warnings()` to see the 4 remaining warnings. 
# Warning messages:
#   1: Removed 5 rows containing non-finite outside the scale range (`stat_boxplot()`). 
# 2: Removed 5 rows containing non-finite outside the scale range (`stat_signif()`). 
# 3: Computation failed in `stat_signif()`.
# Caused by error in `wilcox.test.default()`:
#   ! not enough 'y' observations 
# 4: Removed 5 rows containing missing values or values outside the scale range (`geom_point()`). 
# `geom_smooth()` using formula = 'y ~ x'
# `geom_smooth()` using formula = 'y ~ x'
# Error in Rtsne.default(t(assay(se_filtered, "logcpm") ), perplexity = 1) : 
#   Remove duplicates before running TSNE.
# Calls: Rtsne -> Rtsne.default
# Execution halted