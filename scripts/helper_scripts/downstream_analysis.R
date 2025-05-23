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

# Add ground truth to combined table
truth_data <- read_csv(opts$`ground-truth`, show_col_types = FALSE) 
truth_data$ground_truth <- truth_data$abundance
truth_data$abundance <- NULL

# Join ground truth to combined table, preserving all ground truth species
combined_table <- full_join(combined_table, truth_data, by = "species") %>%
  mutate(ground_truth = ground_truth * 2000) %>% # Scale abundance to avoid filtering out
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))  # Fill missing values

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
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),   # Average numeric values
    across(where(~ !is.numeric(.x)), dplyr::first),        # First non-numeric entry
    .groups = "drop"
  ) %>%
  select(-sample) %>%
  column_to_rownames("run_id")

# Create and add an empty "ground_truth" row
empty_row <- as.data.frame(as.list(rep(NA, ncol(col_data))))
colnames(empty_row) <- colnames(col_data)
rownames(empty_row) <- "ground_truth"
col_data <- rbind(col_data, empty_row)
col_data["ground_truth", "db_name"] <- "ground_truth"

# Fill in missing values and finalize format
col_data <- col_data %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(.x, na.rm = TRUE), .))) %>%
  mutate(
    trim_sliding_win = as.character(trim_sliding_win),  # Force character
    runtime = as.character(runtime),  # Force character
    across(where(is.character), ~ ifelse(is.na(.), dplyr::first(na.omit(.)), .))
  ) %>%
  DataFrame()

# Align run order between count matrix and metadata
ordered_runs <- sort(colnames(count_matrix))
count_matrix <- count_matrix[, ordered_runs]
col_data <- col_data[ordered_runs, , drop = FALSE]
col_data["ground_truth", "db_name"] <- "ground_truth"

# Create SummarizedExperiment object
se <- SummarizedExperiment(
  assays = list(counts = count_matrix),
  rowData = row_data,
  colData = col_data
)

# Filter low quality data 
se <- se %>%
  subset(rowSums(assay(., "counts")) >= 150) %>%  # Keep species with ≥150 total reads across all runs (not per-sample filtering; preserves missing species)
  subset(, colSums(assay(., "counts")) >= 150) %>%  # Keep runs with ≥100 total reads
  subset(, !duplicated(t(assay(., "counts"))))  # Remove duplicated count profiles 

#=========================================================
# Read Length and Size Progression Analysis 
#=========================================================

## ------------------ Read Retention ------------------

# Define pipeline stages (read counts after major steps)
stage_names <- c("trim_paired", "bt_paired", "kraken2_total", "kraken2_classified", "bracken_total")

# Aggregate read counts by database (species presence/absence in DB affects classification)
aggregated_summary <- as_tibble(colData(se), rownames = "run_id") %>%
  select(db_name, all_of(stage_names)) %>%
  group_by(db_name) %>%                                  
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
min_read_data <- as_tibble(colData(se), rownames = "run_id") %>%
  select(db_name, all_of(bp_cols)) %>%  
  filter(!is.na(db_name)) %>%
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

# Impute NAs and compute average read lengths per stage grouped by database
length_summary <- min_read_data %>%
  group_by(db_name) %>%
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
read_length <- ggplot(length_long, aes(x = Stage, y = Length, group = interaction(db_name, Metric), 
                                       color = db_name, linetype = Metric)) +
  geom_line(linewidth = 0.6, alpha = 0.6, na.rm = TRUE) +
  geom_line(data = mean_length_summary, aes(x = Stage, y = Length, group = Metric), 
            color = "grey", linewidth = 1) +
  geom_point(data = mean_length_summary, aes(x = Stage, y = Length), 
             color = "grey", size = 2, inherit.aes = FALSE) +
  geom_text_repel(data = mean_length_summary, aes(x = Stage, y = Length, label = round(Length, 1)), 
                  color = "black", size = 3, inherit.aes = FALSE) +
  labs(title = "Progression of Read Length", y = "Read Length (bp)", x = NULL, linetype = "Metric", 
       color = "Database") +
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
            
ggsave("read_progression.png", read_progression, width = 12, height = 4.5, dpi = 300)

#=========================================================
# Clustering with Ground Truth
#=========================================================

##------------- Filter Count Matrix -------------

# Extract raw counts assay
filtered_mat <- assay(se, "counts")

lib_sizes <- setNames(colData(se)$bracken_total, colnames(se)) # Number of reads classified and redistributed by Bracken
# Effective number of reads that contribute to species-level abundance

#lib_sizes      <- setNames(colData(se)$bt_paired, colnames(se)) # Number of reads remaining after Bowtie2 host filtering
# Reflects the raw metagenomic input available before classification

##------------- Visualize Library Sizes -------------

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

if (length(unique(lib_sizes)) != 1) {
  cor_vals_raw <- apply(filtered_mat, 1, \(x) cor(x, lib_sizes))
  p2 <- tibble(Correlation = cor_vals_raw) %>%
    ggplot(aes(x = Correlation)) +
    geom_histogram(binwidth = 0.05, fill = viridis(1, option = "D"), color = "black") +
    stat_function(fun = dnorm, args = list(mean = mean(cor_vals_raw), sd = sd(cor_vals_raw)), 
                  color = "grey", linewidth = 1) + # Normal distribution curve 
    labs(title = "Species Abundance Correlation with Library Size (Raw)", x = "Correlation", y = "Density") +
    theme_minimal() + 
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 8), axis.title = element_text(size = 9),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  # CPM normalization (adjusts for library size; Zymo is balanced → no genome correction)
  assay(se, "cpm") <- t(t(filtered_mat) / lib_sizes * 1e6)  # Add as new assay
  
  # Correlation after CPM normalization (bias should be reduced)
  cor_vals_cpm <- apply(assay(se, "cpm"), 1, \(x) cor(x, lib_sizes))
  p3 <- tibble(Correlation = cor_vals_cpm) %>%
    ggplot(aes(x = Correlation)) +
    geom_histogram(binwidth = 0.05, fill = viridis(1, option = "D"), color = "black") +
    stat_function(fun = dnorm, args = list(mean = mean(cor_vals_cpm), sd = sd(cor_vals_cpm)),
                  color = "grey", linewidth = 1) + # Normal distribution curve 
    labs(title = "Species Abundance Correlation with Library Size (CPM)", x = "Correlation", y = "Density") +
    theme_minimal() + 
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 8), axis.title = element_text(size = 9),
          axis.text.x = element_text(angle = 45, hjust = 1))
} else {
  # CPM normalization (adjusts for library size; Zymo is balanced → no genome correction)
  assay(se, "cpm") <- t(t(filtered_mat) / lib_sizes * 1e6)  # Add as new assay
} 
  
## ------------- Log Transformation -------------

# Plot log(Mean) vs log(SD) to examine heteroscedasticity; log(+1) avoids issues with 0s
df_var_cpm <- data.frame(Mean = log1p(rowMeans(assay(se, "cpm"))), 
                         SD = log1p(rowSds(assay(se, "cpm"))))
p4 <- ggplot(df_var_cpm, aes(x = Mean, y = SD, color = viridis(1, option = "D"))) +
  geom_point(size = 2, color = viridis(1, option = "D")) +  
  geom_smooth(formula = 'y ~ x', method = "lm", se = FALSE, color = "grey", linewidth = 1) +
  annotate("text", x = min(df_var_cpm$Mean), y = max(df_var_cpm$SD),
           label = paste0("r = ", round(cor(df_var_cpm$Mean, df_var_cpm$SD), 3)), 
           hjust = 0, size = 3, fontface = "italic") +
  labs(x = "log(Mean CPM + 1)", y = "log(SD CPM + 1)", title = "Species Abundance Variability (CPM)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Variance stabilization: log2 transform CPM 
assay(se, "logcpm") <- log2(assay(se, "cpm")+1)

# Analyze mean–SD relationship again
df_var_logcpm <- data.frame(Mean = rowMeans(assay(se, "logcpm")), 
                            SD = rowSds(assay(se, "logcpm")))
p5 <- ggplot(df_var_logcpm, aes(x = Mean, y = SD)) +
  geom_point(size = 2, color = viridis(1, option = "D")) +  
  geom_smooth(formula = 'y ~ x', method = "lm", se = FALSE, color = "grey", linewidth = 1) +
  annotate("text", x = max(df_var_logcpm$Mean) - 0.2, y = max(df_var_logcpm$SD),
           label = paste0("r = ", round(cor(df_var_logcpm$Mean, df_var_logcpm$SD), 3)), 
           hjust = +1, size = 3, fontface = "italic") +
  labs(x = "Mean logCPM", y = "SD logCPM", title = "Species Abundance Variability (logCPM)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Combine and export figure
if (length(unique(lib_sizes)) != 1) {
  combined_transformations <- p1 / (p2 + p3) / (p4 + p5) 
  ggsave("combined_transformations_varied.png", combined_transformations, width = 10, height = 10) 
} else {
    combined_transformations <- p1 / (p4 + p5) 
    ggsave("combined_transformations_constant.png", combined_transformations, width = 10, height = 10) 
}

#=========================================================
# Clustering with Ground Truth
#=========================================================

##------------- Dimensionality Reduction -------------

# TODO: set the database for ground_truth to "ground_truth" and color the clusters per database 

# PCA plot: projects runs based on species composition in logCPM space
pca_df <- prcomp(t(assay(se, "logcpm")), scale. = TRUE)$x %>%  # Run PCA and extract principal components
  as.data.frame() %>% 
  rownames_to_column("run_id") # # Add run labels

p6 <- ggplot(pca_df, aes(x = PC1, y = PC2, label = run_id, color = run_id)) +
  geom_point(aes(color = run_id == "ground_truth"), size = 3, show.legend = FALSE) +
  geom_text_repel(aes(color = run_id == "ground_truth"), size = 3, show.legend = FALSE, max.overlaps = Inf) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
  labs(title = "PCA", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# t-SNE (optimized for local structure)
perplexity <- min(30, floor(( ncol(assay(se, "logcpm")) - 1) / 3)) # Set perplexity (30 or lower if needed)
tsne_df <- Rtsne(t(assay(se, "logcpm")), perplexity = perplexity)$Y %>% # Run t-SNE on transposed logCPM matrix
  as.data.frame() %>% 
  setNames(c("tsne_1", "tsne_2")) %>% 
  mutate(run_id = colnames(se)) # Add run labels

p7 <- ggplot(tsne_df, aes(x = tsne_1, y = tsne_2, label = run_id, color = run_id)) +
  geom_point(aes(color = run_id == "ground_truth"), size = 3, show.legend = FALSE) +
  geom_text_repel(aes(color = run_id == "ground_truth"), size = 3, show.legend = FALSE, max.overlaps = Inf) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
  labs(title = "t-SNE", subtitle = paste("Perplexity:", perplexity), x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1))

# UMAP (preserves global and local structure)
n_neighbors <- min(15, max(2, floor(ncol(assay(se, "logcpm")) / 2)))  # Compute safe neighbor count
umap_cfg <- umap.defaults; umap_cfg$n_neighbors <- n_neighbors                 # Update config
umap_df <- umap(t(assay(se, "logcpm")), config = umap_cfg)$layout %>% # Run UMAP and extract results
  as.data.frame() %>% setNames(c("umap_1", "umap_2")) %>% 
  mutate(run_id = colnames(assay(se, "logcpm"))) # Add run labels

p8 <- ggplot(umap_df, aes(x = umap_1, y = umap_2, label = run_id, color = run_id)) + 
  geom_point(aes(color = run_id == "ground_truth"), size = 3, show.legend = FALSE) +
  geom_text_repel(aes(color = run_id == "ground_truth"), size = 3, show.legend = FALSE, max.overlaps = Inf) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
  labs(title = "UMAP", subtitle = paste("N° Neighbors:", n_neighbors), x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1))

##------------- Distance Computation -------------

# TODO: eventually highlight ground truth label in red ? 

# Helper: wrap pheatmap as a ggplot-compatible grob
pheatmap_grob <- function(mat, show_legend=TRUE) {
  p <- pheatmap(mat,clustering_distance_rows = "euclidean",clustering_distance_cols = "euclidean",
    color = viridis(100, option = "D"),angle_col = 45,silent = TRUE,show_colnames = FALSE, 
    legend = show_legend)
  as.ggplot(p[[4]]) 
}

# Distances in original logCPM space (species abundance profiles)
raw_dists   <- as.matrix(dist(t(assay(se, "logcpm"))))
p9  <- pheatmap_grob(raw_dists) # Plot as heatmap 
ggsave("log2cpm_distances.png", p9, width = 12, height = 8) 

# Distances in PCA space (first 2 principal components)
rownames(pca_df) <- pca_df$run_id
pca_dists <- pca_df[, c("PC1", "PC2")] %>%
  dist() %>% # Euclidean distance between runs
  as.matrix()
pca_dists <- pca_dists / max(pca_dists) # Normalize to [0, 1] range for scale consistency
p10 <- pheatmap_grob(pca_dists, show_legend=FALSE) # Plot as heatmap 

# Distances in UMAP space
rownames(umap_df) <- umap_df$run_id
umap_dists  <- as.matrix(dist(umap_df[, c("umap_1", "umap_2")]))
umap_dists <- umap_dists / max(umap_dists) # Normalize to [0, 1] range for scale consistency
p11 <- pheatmap_grob(umap_dists, show_legend=FALSE) # Plot as heatmap 

# Distances in t-SNE space
rownames(tsne_df) <- tsne_df$run_id
tsne_dists  <- as.matrix(dist(tsne_df[, c("tsne_1", "tsne_2")]))
tsne_dists <- tsne_dists / max(tsne_dists)
p12 <- pheatmap_grob(tsne_dists)

if (length(unique(lib_sizes)) != 1) {
  combined_plot <- (p6 + p7 + p8) / (p10 + p11 + p12)
  ggsave("combined_distances_constant.png", combined_plot, width = 12, height = 8) 
} else {
  combined_plot <- (p6 + p7 + p8) / (p10 + p11 + p12)
  ggsave("combined_distances_varied.png", combined_plot, width = 12, height = 8) 
}


# TODO: Use distances from "truth" column to rank best-matching configurations

#=========================================================
# Precision and Recall Graphs 
#=========================================================

#gt_species <- unique(truth_data$species)
#gt_species
