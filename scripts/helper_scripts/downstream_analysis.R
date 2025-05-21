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

  # Read report and compute mean across numeric columns (handle NAs)
  run_species <- read_csv(file_path, show_col_types = FALSE) %>%
    mutate(!!run_id := rowMeans(across(where(is.numeric), ~ replace_na(.x, 0)))) %>%
    select(species, !!sym(run_id)) # Keep only species and new column

  # Merge into the combined table
  combined_table <- if (nrow(combined_table) == 0) {
    run_species
  } else {
    full_join(combined_table, run_species, by = "species") %>%
      mutate(across(where(is.numeric), ~ replace_na(.x, 0)))  # Fill missing values with 0
  }
}

# Convert combined table to matrix format with species as rownames
count_matrix <- combined_table %>%
  column_to_rownames("species") %>%
  as.matrix()

## ------------- Assemble metadata ------------- 

# Create row metadata for species
row_data <- DataFrame(species = rownames(count_matrix))

# Read and summarize run metadata (mean for numeric, first for categorical)
col_data <- read_csv(opts$`runs-summary`, show_col_types = FALSE) %>%
  group_by(run_id) %>%
  summarise(
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)) , # Average numeric values
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

# Define columns representing total read counts at key pipeline stages
stage_names <- c("trim_paired", "bt_paired", "kraken2_total", "kraken2_classified", "bracken_total")

# Compute mean values to impute missing data
mean_trim <- mean(col_data$trim_paired[!is.nan(col_data$trim_paired)], na.rm = TRUE)
mean_bt   <- mean(col_data$bt_paired[!is.nan(col_data$bt_paired)], na.rm = TRUE)

# Aggregate read counts by database
aggregated_summary <- as_tibble(col_data, rownames = "run_id") %>%
  select(db_name, all_of(stage_names)) %>%                           
  group_by(db_name) %>%
  mutate(
    trim_paired = if_else(is.nan(trim_paired), mean_trim, trim_paired),
    bt_paired   = if_else(is.nan(bt_paired), mean_bt, bt_paired)
  ) %>%
  summarise(across(all_of(stage_names), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# Transform data for plotting and normalize to initial stage
long_summary <- aggregated_summary %>%
  pivot_longer(cols = all_of(stage_names), names_to = "Original", values_to = "Reads") %>%
  mutate(Stage = factor(Original, levels = stage_names)) %>%  # Correctly assign Stage
  group_by(db_name = db_name) %>%
  mutate(Fraction = Reads / Reads[Original == "trim_paired"]) %>%
  mutate(Stage = recode(Stage,
                        "trim_paired"        = "Raw Reads",
                        "bt_paired"          = "Trimmed Reads",
                        "kraken2_total"      = "Host Filtered",
                        "kraken2_classified" = "Classified (Kraken)",
                        "bracken_total"      = "Classified (Bracken)"
  )) %>%
  ungroup()

# Compute average proportion per stage
mean_summary <- long_summary %>%
  group_by(Stage) %>%
  summarise(Fraction = mean(Fraction, na.rm = TRUE), .groups = "drop")

# Plot read retention over pipeline stages
read_retention <- ggplot(long_summary, aes(x = Stage, y = Fraction, group = db_name, color = db_name)) +
  geom_line(linewidth = 0.5, linetype = "dashed") +                     # sample lines
  geom_line(data = mean_summary, aes(x = Stage, y = Fraction, group = 1), 
            color = "black", linewidth = 1, inherit.aes = FALSE) +     
  geom_point(data = mean_summary, aes(x = Stage, y = Fraction), 
             color = "black", size = 2, inherit.aes = FALSE) +     
  geom_text_repel(data = mean_summary, 
            aes(x = Stage, y = Fraction, label = scales::percent(Fraction, accuracy = 0.01)),
            color = "black", size = 3.5, vjust = -0.6, hjust = -0.3, inherit.aes = FALSE) +
  theme_bw() +     
  scale_color_brewer(palette = "Dark2") +      
  labs(title = "Progression of Read Count", y = "Proportion of Raw Reads", x = NULL, color = "Database") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")      

## ------------------ Read Length Progression ------------------

# Define raw read length columns
bp_cols <- c("preqc_avg_len_r1", "preqc_avg_len_r2", "preqc_med_len_r1", "preqc_med_len_r2",
  "postqc_avg_len_r1", "postqc_avg_len_r2", "postqc_med_len_r1", "postqc_med_len_r2",
  "bt_avg_r1", "bt_avg_r2", "bt_med_r1", "bt_med_r2",
  "clas_avg_r1", "clas_avg_r2", "clas_med_r1", "clas_med_r2",
  "unclas_avg_r1", "unclas_avg_r2", "unclas_med_r1", "unclas_med_r2")

# Derive min average/median per stage per run
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

# Impute missing values with global means
mean_preqc_avg <- mean(min_read_data$preqc_avg[!is.nan(min_read_data$preqc_avg)], na.rm = TRUE)
mean_preqc_med  <- mean(min_read_data$preqc_med[!is.nan(min_read_data$preqc_med)], na.rm = TRUE)
mean_postqc_avg <- mean(min_read_data$postqc_avg[!is.nan(min_read_data$postqc_avg)], na.rm = TRUE)
mean_postqc_med  <- mean(min_read_data$postqc_med[!is.nan(min_read_data$postqc_med)], na.rm = TRUE)

# Aggregate mean/median lengths per stage by database
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

# Reshape to long format for plotting
length_long <- length_summary %>%
  pivot_longer(cols =-c(db_name, unclas_avg, unclas_med), names_to = "Stage_Metric", values_to = "Length") %>%
  separate(Stage_Metric, into = c("Stage", "Metric"), sep = "_") %>%
  mutate(Stage = factor(Stage,levels = c("preqc", "postqc", "bt", "clas"),
                        labels = c("Raw Reads", "Trimmed Reads", "Host Filtered", "Classified Reads")),
         Metric = if_else(Metric == "avg", "Mean", "Median"))

# Compute average read length per stage and metric
mean_length_summary <- length_long %>%
  group_by(Stage, Metric) %>%
  summarise(Length = mean(Length, na.rm = TRUE), .groups = "drop")

# Plot progression of read length
read_length <- ggplot(length_long, aes(x = Stage, y = Length, group = interaction(db_name, Metric), 
                                       color = db_name, linetype = Metric)) +
  geom_line(linewidth = 0.6, alpha = 0.6, na.rm = TRUE) +
  geom_line(data = mean_length_summary, aes(x = Stage, y = Length, group = Metric), 
            color = "black", linewidth = 1, inherit.aes = FALSE) +
  geom_point(data = mean_length_summary, aes(x = Stage, y = Length), 
             color = "black", size = 2, inherit.aes = FALSE) +
  geom_text_repel(data = mean_length_summary,
            aes(x = Stage, y = Length, label = round(Length, 1)), color = "black", 
            size = 3.2, vjust = -0.6,
            inherit.aes = FALSE) +
  labs(title = "Progression of Read Length", y = "Read Length (bp)", x = NULL, linetype = "Metric") +
  scale_linetype_manual(values = c("Mean" = "dashed", "Median" = "solid")) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +    
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## ------------------ Boxplot: Classified vs Unclassified ------------------

# Prepare long format for median classified/unclassified read length
clas_unclas <- length_summary %>%
  select(db_name, clas_med, unclas_med) %>%  
  pivot_longer(cols = c(clas_med, unclas_med), names_to = "Type", values_to = "Reads") %>%  
  mutate(Type = recode(Type,                
                       "clas_med" = "Classified",
                       "unclas_med" = "Unclassified"))

# Boxplot with Wilcoxon test
boxplot <- ggplot(clas_unclas, aes(x = Type, y = Reads, fill = Type)) +
  geom_boxplot(alpha = 0.8) +                                 
  geom_jitter(width = 0.1, size = 2, color = "black") +       
  stat_compare_means(comparisons = list(c("Classified", "Unclassified")), 
                     method = "wilcox.test", label = "p.format") +  
  scale_fill_brewer(palette = "Set2") +                       
  labs(title = "Median Read Length Comparison", x = NULL, y = "Read Length (bp)") +
  theme_bw() +                                                
  theme(legend.position = "none", axis.text.x = element_text(angle = 15, hjust = 1))

# Save combined read progression figure
read_progression <- boxplot + read_retention + read_length
ggsave("read_progression.png", read_progression, width = 12, height = 6, dpi = 300)

#=========================================================
# Clustering with Ground Truth
#=========================================================

##------------- Filter Count Matrix -------------

# Extract count matrix
mat <- assay(se, "counts")

# Keep species with ≥20 reads in ≥50% of samples
min_reads <- 20
min_samples <- ncol(se) / 2
rows_to_keep <- rowSums(mat >= min_reads) >= min_samples
se_filtered <- se[rows_to_keep, ]

# Keep samples with ≥100 total reads
filtered_mat <- assay(se_filtered, "counts")
cols_to_keep <- colSums(filtered_mat) >= 100
se_filtered <- se_filtered[, cols_to_keep]

##------------- Visualize Library Sizes -------------

# Compute total reads per sample
filtered_mat <- assay(se_filtered, "counts")
lib_sizes <- colSums(filtered_mat)
df_lib <- data.frame(Run = names(lib_sizes), Reads = lib_sizes)

# Barplot: total classified reads per sample
p1 <- ggplot(df_lib, aes(x = Run, y = Reads, fill = Run)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  labs(x = "Pipeline Run", y = "Classified Reads", title = "Classified Reads per Pipeline Run") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##------------- Normalize Count Matrix -------------

# Correlation between raw species counts and library size
cor_vals_raw <- apply(filtered_mat, 1, \(x) cor(x, lib_sizes))
df_cor_raw <- data.frame(Correlation = cor_vals_raw)
p2 <- ggplot(df_cor_raw, aes(x = Correlation)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  stat_function(fun = dnorm, 
                args = list(mean = mean(cor_vals_raw), sd = sd(cor_vals_raw)), 
                color = "red", linewidth = 1) +
  labs(x = "Correlation Coefficient", y = "Density", title = "Raw Abundance vs. Library Size") +
  theme_minimal()

# Convert to CPM
filtered_cpm <- t(t(filtered_mat) / lib_sizes * 1e6)

# Correlation after CPM normalization
cor_vals_cpm <- apply(filtered_cpm, 1, \(x) cor(x, lib_sizes))
df_cor_cpm <- data.frame(Correlation = cor_vals_cpm)
p3 <- ggplot(df_cor_cpm, aes(x = Correlation)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  stat_function(fun = dnorm, 
                args = list(mean = mean(cor_vals_cpm), sd = sd(cor_vals_cpm)), 
                color = "red", linewidth = 1) +
  labs(x = "Correlation (CPM)", y = "Number of Species", title = "CPM-Normalized Correlation") +
  theme_minimal()

## ------------- Log Transformation -------------

# Log mean vs SD (CPM)
df_var_cpm <- data.frame(Mean = log1p(rowMeans(filtered_cpm)), SD = log1p(rowSds(filtered_cpm)))
cor_cpm <- cor(df_var_cpm$Mean, df_var_cpm$SD)
p4 <- ggplot(df_var_cpm, aes(x = Mean, y = SD)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
  annotate("text", x = min(df_var_cpm$Mean), y = max(df_var_cpm$SD), 
           label = paste0("r = ", round(cor_cpm, 3)), hjust = 0, size = 4, fontface = "italic") +
  labs(x = "log(Mean CPM + 1)", y = "log(SD CPM + 1)", title = "Species Variability (CPM)") +
  theme_minimal()

# Log2 CPM
filtered_logcpm <- log2(filtered_cpm+1)

# Log mean vs SD (logCPM)
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

##------------- Dimensionality Reduction -------------

# PCA
pca <- prcomp(t(filtered_logcpm), scale. = TRUE)
pca_df <- as.data.frame(pca$x) %>% rownames_to_column("run_id")
p6 <- ggplot(pca_df, aes(x = PC1, y = PC2, label = sample, color = sample)) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = Inf) +
  labs(title = "PCA", x = "PC1", y = "PC2") +
  theme_minimal()

# t-SNE
tsne <- Rtsne(t(filtered_logcpm), perplexity = 1)
tsne_df <- as.data.frame(tsne$Y)
colnames(tsne_df) <- c("tsne_1", "tsne_2")
tsne_df$run_id <- colnames(filtered_logcpm)
p7 <- ggplot(tsne_df, aes(x = tsne_1, y = tsne_2, label = sample, color = run_id)) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = Inf) +
  labs(title = "t-SNE", x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()

# UMAP
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
# Error in Rtsne.default(t(filtered_logcpm), perplexity = 1) : 
#   Remove duplicates before running TSNE.
# Calls: Rtsne -> Rtsne.default
# Execution halted