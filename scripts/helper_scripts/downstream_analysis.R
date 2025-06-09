#!/usr/bin/env Rscript

# This script aggregates Bracken output and associated metadata into a SummarizedExperiment object,
# and performs downstream analysis of metagenomic pipeline performance. It generates:
# - Read retention and read length progression plots across pipeline stages,
# - Clustering visualizations of species profiles (PCA, t-SNE, UMAP),
# - Comparative analyses of classification accuracy against ground truth (L2, AUPR, precision/recall).

# Usage: Rscript downstream_analysis.R -t <path/to/ground_truth.csv> -s <path/to/runs_summary.csv> breport1.csv breport2.csv ...

# Rscript downstream_analysis.R \
# -t /Users/aliciahobdell/Desktop/final_project/zymobiomics_folder/raw_data/ground_truth.csv \
# -s /Users/aliciahobdell/Desktop/final_project/zymobiomics_folder/results/runs/runs_summary.csv \
# $(find /Users/aliciahobdell/Desktop/final_project/zymobiomics_folder/results/runs/ -type f -name "combined_breports.csv" | sort)

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
  library(scales)               # Axis scaling and label formatting
  library(precrec)              # Precision-Recall and ROC curve evaluation
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

# Define output directory for plots 
results_dir <- file.path(dirname(opts$`runs-summary`))

#=========================================================
# Creating Summarized Experiment 
#=========================================================

## ------------- Load Bracken reports ------------- 
# Load and merge combined Bracken reports (one per run) into a single species-by-run table
# Note: These reports are already aggregated per run — no sample-level detail is retained

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

## ------------- Add Ground Truth ------------- 
# Merge ground truth profile into the species table and prepare a count matrix for SummarizedExperiment

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

## ------------- Assemble Metadata ------------- 
# Extract species (row) and run-level (column) metadata

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

cols_to_fill <- c(
  "trim_clip", "trim_head", "trim_lead", "trim_crop", "trim_sliding_win", "trim_trail",
  "trim_min_len", "trim_paired", "trim_both", "trim_fw_only", "trim_rev_only", "trim_dropped",
  "bt_prefix", "bt_mode", "bt_sensitivity", "bt_mixed", "bt_discordant",
  "bt_paired", "bt_conc_0", "bt_conc_1", "bt_conc_more"
)

# Fill in missing values and finalize format
col_data <- col_data %>%
  mutate(
    across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)),
    trim_sliding_win = as.character(trim_sliding_win),  # Force character
    runtime = as.character(runtime)                       # Force character
  ) %>%
  fill(all_of(cols_to_fill), .direction = "down") %>%    # Forward fill NAs in specified columns
  DataFrame()

# Align run order between count matrix and metadata
ordered_runs <- sort(colnames(count_matrix))
count_matrix <- count_matrix[, ordered_runs]
col_data <- col_data[ordered_runs, , drop = FALSE]
col_data["ground_truth", "db_name"] <- "ground_truth"

## ------------- Build & Filter SummarizedExperiment ------------- 

# Create SummarizedExperiment object
se <- SummarizedExperiment(
  assays = list(counts = count_matrix),
  rowData = row_data,
  colData = col_data
)

# Filter low quality data 
se <- se %>%
  subset(rowSums(assay(., "counts")) >= 200) %>%  # Keep species with ≥150 total reads across all runs (not per-sample filtering; preserves missing species)
  subset(, colSums(assay(., "counts")) >= 300) %>%  # Keep runs with ≥100 total reads
  subset(, !duplicated(t(assay(., "counts"))))  # Remove duplicated count profiles 

# Create global legend
legend_df <- unique(colData(se)[, "db_name", drop = FALSE])  # Drop = FALSE to keep as dataframe
db_colors <- setNames(viridis(nrow(legend_df), option = "D"), legend_df$db_name)
db_colors["ground_truth"] <- "red"

#=========================================================
# Read Length and Size Progression Analysis 
#=========================================================

## ------------------ Read Retention ------------------
# Plot how the number of reads changes across pipeline stages (trimming, host removal, classification)

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
  filter(db_name != "ground_truth") %>%
  ungroup() 

# Mean % retained at each stage across all DBs
mean_summary <- long_summary %>%
  group_by(Stage) %>%
  summarise(Fraction = mean(Fraction, na.rm = TRUE), .groups = "drop")

# Plot read retention progression across pipeline stages
read_retention <- ggplot(long_summary, aes(x = Stage, y = Fraction, group = db_name, color = db_name)) +
  geom_line(linewidth = 0.5) + # Per-database trend
  geom_line(data = mean_summary, aes(x = Stage, y = Fraction, group = 1), 
            color = "grey", linewidth = 1, inherit.aes = FALSE) + # Mean trend 
  geom_point(data = mean_summary, aes(x = Stage, y = Fraction), 
             color = "grey", size = 2, inherit.aes = FALSE) + # Mean values 
  geom_text_repel(data = mean_summary, aes(x = Stage, y = Fraction, label = percent(Fraction, accuracy = 0.01)), 
                  color = "black", size = 3, inherit.aes = FALSE) + # Add % labels to mean
  scale_color_manual(values = db_colors) +  # Use global color mapping
  labs(title = "Progression of Read Count", y = "Proportion of Raw Reads", x = NULL, color = "Database") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

## ------------------ Read Length Progression ------------------
# Plot how the read length changes across pipeline stages (trimming, host removal, classification)

# Define pipeline stages (read counts after major steps)
bp_cols <- c("preqc_avg_len_r1", "preqc_avg_len_r2", "preqc_med_len_r1", "preqc_med_len_r2",
  "postqc_avg_len_r1", "postqc_avg_len_r2", "postqc_med_len_r1", "postqc_med_len_r2",
  "bt_avg_r1", "bt_avg_r2", "bt_med_r1", "bt_med_r2",
  "clas_avg_r1", "clas_avg_r2", "clas_med_r1", "clas_med_r2",
  "unclas_avg_r1", "unclas_avg_r2", "unclas_med_r1", "unclas_med_r2")

# Calculate per-run min read length per stage (across R1/R2) to account for asymmetry
min_read_data <- as_tibble(colData(se), rownames = "run_id") %>%
  select(db_name, all_of(bp_cols)) %>%  
  filter(db_name != "ground_truth") %>%
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
  geom_line(linewidth = 0.6, na.rm = TRUE) +
  geom_line(data = mean_length_summary, aes(x = Stage, y = Length, group = Metric), 
            color = "grey", linewidth = 1) +
  geom_point(data = mean_length_summary, aes(x = Stage, y = Length), 
             color = "grey", size = 2, inherit.aes = FALSE) +
  geom_text_repel(data = mean_length_summary, aes(x = Stage, y = Length, label = round(Length, 1)), 
                  color = "black", size = 3, inherit.aes = FALSE) +
  labs(title = "Progression of Read Length", y = "Read Length (bp)", x = NULL, linetype = "Metric", 
       color = "Database") +
  scale_linetype_manual(values = c("Median" = "dashed", "Mean" = "solid")) +
  scale_color_manual(values = db_colors) +  # Use global color mapping
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8))

## ------------------ Boxplot: Classified vs Unclassified ------------------
# Compare read lengths between classified and unclassified reads using a boxplot

# Reshape median classified/unclassified read lengths into long format for plotting
clas_unclas <- length_summary %>%
  select(db_name, clas_med, unclas_med) %>%  # Keep relevant columns
  pivot_longer(cols = c(clas_med, unclas_med), names_to = "Type", values_to = "Reads") %>%
  mutate(Type = recode(Type, "clas_med" = "Classified", "unclas_med" = "Unclassified"))  # Rename for clarity

# Boxplot comparing median read lengths of classified vs unclassified reads
# Unpaired Wilcoxon test: non-parametric (handles skewed data), compares medians of independent groups (classified vs unclassified reads)
boxplot <- ggplot(clas_unclas, aes(x = Type, y = Reads, fill = Type)) +
  geom_boxplot(alpha = 0.7, color = "gray", size = 0.5) + 
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
            
ggsave(file.path(results_dir, "read_progression.png"), read_progression, width = 12, height = 4.5, dpi = 300)

#=========================================================
# Data Preprocessing
#=========================================================

#------------- Compare Normalizations -------------
# This section compares L2 distance to ground truth using two normalization methods:
# - Total reads reflect true detection rates but penalize under-classified samples.
# - Classified reads emphasize relative proportions among detected taxa, but can inflate abundances if classification rates vary.
# If classification bias is uniform across taxa, this inflation may cancel out, preserving relative accuracy.

# Get species count matrix from SummarizedExperiment
counts_mat <- assay(se, "counts")

# Get list of run IDs
runs <- colnames(counts_mat)

# Normalize by total reads (bt_paired)
lib_sizes <- setNames(colData(se)$bt_paired, runs)
norm_total <- sweep(counts_mat, 2, lib_sizes, FUN = "/")

# Normalize by classified reads (bracken_total)
class_sizes <- setNames(colData(se)$bracken_total, runs)
norm_class <- sweep(counts_mat, 2, class_sizes, FUN = "/")

# Extract and prepare ground truth profile
gt_profile <- norm_total[, "ground_truth"]
gt_profile <- gt_profile / sum(gt_profile)  # Normalize to relative abundance
norm_total <- norm_total[, runs != "ground_truth"]
norm_class <- norm_class[, runs != "ground_truth"]

# Define L2 distance to ground truth
l2_dist <- function(run_col) sqrt(sum((run_col - gt_profile)^2))
df_l2 <- tibble(
  Run = rep(colnames(norm_total), 2),
  L2 = c(
    apply(norm_class, 2, l2_dist), # Distances using classified read normalization
    apply(norm_total, 2, l2_dist) #  Distances using total read normalization
  ),
  Normalization = rep(c("Classified", "Total"), each = ncol(norm_total))
) %>%
  mutate(Database = colData(se)[Run, "db_name"])

# Plot L2 distances by database
p1 <- ggplot(df_l2, aes(x = Database, y = L2, fill = Normalization)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(aes(group = Normalization), position = position_dodge(width = 0.8), size = 1, color = "gray") +
  labs(title = "Impact of Normalization on Profile Accuracy", y = "L2 Distance to Ground Truth", x = "Database") +
  scale_color_manual(values = db_colors) +  # Use global color mapping
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1))

#------------- Visualize Library Sizes (Classified) -------------

# Barplot: total classified reads per run
p2 <- data.frame(Run = names(class_sizes), Reads = class_sizes) %>%
  filter(Run != "ground_truth") %>%
  ggplot(aes(x = Run, y = Reads, fill = Run)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  scale_color_manual(values = db_colors) +  # Use global color mapping
  labs(title = "Classified Reads per Pipeline Run", x = "Pipeline Run", y = "Classified Reads") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1))

#------------- Data Normalization -------------
# This section examines the correlation of species abundance with classification depth before and after CPM normalization (ground truth excluded).
# CPM adjusts for number of classified reads; TPM is unnecessary as ZymoBIOMICS is genome-balanced.

# Remove ground truth from class sizes
class_sizes_nogt <- class_sizes[names(class_sizes) != "ground_truth"] 

# Prepare count matrix 
counts_mat_nogt <- counts_mat %>%
  .[, colnames(.) != "ground_truth"] %>%      # Exclude ground truth column
  .[rowSums(.) > 0, ]                         # Exclude species with zero counts across all samples

# Compute correlation between species abundance and classified read count (pre-CPM)
cor_vals_raw <- apply(counts_mat_nogt, 1, \(x) cor(x, class_sizes_nogt))

# Plot raw correlation distribution
p3 <- tibble(Correlation = cor_vals_raw) %>%
  ggplot(aes(x = Correlation)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = viridis(1, option = "D"), color = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(cor_vals_raw), sd = sd(cor_vals_raw)), 
                color = "grey", linewidth = 1) + # Normal distribution curve 
  scale_x_continuous(limits = c(-1, 1)) +
  labs(title = "Raw Correlation Distribution", x = "Correlation", y = "Density") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Apply CPM normalization using classified read depth
assay(se, "cpm") <- t(t(counts_mat) / class_sizes * 1e6)  # Add as new assay
cpm_mat <- assay(se, "cpm")

# Prepare normalized matrix
cpm_mat_nogt <- cpm_mat %>%
  .[, colnames(.) != "ground_truth"] %>%      # Exclude ground truth column
  .[rowSums(.) > 0, ]                         # Exclude species with zero CPM across all samples

# Compute correlation after CPM normalization
cor_vals_cpm <- apply(cpm_mat_nogt, 1, \(x) cor(x, class_sizes_nogt))

# Plot CPM-normalized correlation distribution
p4 <- tibble(Correlation = cor_vals_cpm) %>%
  ggplot(aes(x = Correlation)) +
  geom_histogram(binwidth = 0.05, boundary = 0,fill = viridis(1, option = "D"), color = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(cor_vals_cpm), sd = sd(cor_vals_cpm)),
                color = "grey", linewidth = 1) + # Normal distribution curve 
  scale_x_continuous(limits = c(-1, 1)) +
  labs(title = "CPM-Normalized Correlation Distribution", x = "Correlation", y = "Density") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1))

## ------------- Variance Stabilisation -------------
# Plot log(Mean) vs log(SD) to examine heteroscedasticity (when the variance of a variable depends on its mean)

# Compute log-transformed mean and standard deviation per species
df_var_cpm <- data.frame(Mean = log1p(rowMeans(cpm_mat_nogt)), 
                         SD = log1p(rowSds(cpm_mat_nogt)))

# Scatter plot before
p5 <- ggplot(df_var_cpm, aes(x = Mean, y = SD, color = viridis(1, option = "D"))) +
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
assay(se, "logcpm") <- log2(cpm_mat+1) # Log(+1) avoids issues with 0s
logcpm_mat <- assay(se, "logcpm")

# Prepare normalized matrix
logcpm_mat_nogt <- logcpm_mat %>%
  .[, colnames(.) != "ground_truth"] %>%      # Exclude ground truth column
  .[rowSums(.) > 0, ]                         # Exclude species with zero CPM across all samples

# Analyze mean–SD relationship again
df_var_logcpm <- data.frame(Mean = rowMeans(logcpm_mat_nogt), 
                            SD = rowSds(logcpm_mat_nogt))
# Scatter plot after
p6 <- ggplot(df_var_logcpm, aes(x = Mean, y = SD)) +
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

combined_transformations <- (p1 + p2) / (p3 + p4) / (p5 + p6) 
ggsave(file.path(results_dir, "data_preprocessing.png"), combined_transformations, width = 10, height = 10)

#=========================================================
# Database Precision and Recall 
#=========================================================

#------------- Define Labels & Scores -------------

# Define ground truth species
ground_truth_species <- truth_data$species

# Prepare scores per run (exclude ground_truth column)
score_list <- logcpm_mat[, colnames(logcpm_mat) != "ground_truth"] %>%
  as.data.frame() %>%
  as.list()

# Create binary label vector: 1 if species in ground truth, else 0
label_list <- rep(list(as.integer(rownames(logcpm_mat) %in% ground_truth_species)), length(score_list))

# Run identifiers (used as model names)
modnames <- setdiff(colnames(logcpm_mat), "ground_truth")

# Format data for precrec
msmdat <- mmdata(score_list, label_list, modnames = modnames)

#------------- Evaluate Performance -------------

# Calculate ROC and PRC curves
sscurves <- evalmod(msmdat)

# Compute AUCs for all curves
aucs_df <- auc(sscurves)

# Build mapping: run_id ↔︎ database name
db_df <- tibble(modname = colnames(se), Database = as.character(colData(se)$db_name))

# Extract AUPR and match to databases
aupr_df <- aucs_df %>%
  filter(curvetypes == "PRC") %>% 
  left_join(db_df, by = c("modnames" = "modname"))

#------------- Plot AUPR per Database -------------

# Compute median AUPR per database
median_labels <- aupr_df %>%
  group_by(Database) %>%
  summarise(Median = median(aucs, na.rm = TRUE), .groups = "drop")

# Boxplot of AUPR per database
p_aupr <- ggplot(aupr_df, aes(x = Database, y = aucs, fill = Database)) +
  geom_boxplot(size = 0.2) + 
  geom_point(size = 1, color = "gray") +
  scale_fill_manual(values = db_colors) +  # Use global color mapping
  labs(title = "AUPR per Database", x = "Database", y = "AUPR") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_blank(),
        legend.position = "none")

#------------- Plot PRC and ROC Curves -------------

# Convert curves to long-format dataframe
curve_df <- as.data.frame(sscurves) 

# Add database and median AUPR to each run
curve_df <- curve_df %>% 
  left_join(db_df, by = "modname") %>%
  left_join(median_labels, by = "Database") %>%
  mutate(LegendLabel = paste0(Database, " (AUPR=", round(Median, 3), ")"))

# Precision-Recall Curves
p_prc <- ggplot(curve_df %>% filter(type == "PRC"), aes(x = x, y = y, group = modname, color = Database)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = db_colors) +  # Use global color mapping
  labs(title = "Precision-Recall Curves", x = "Recall", y = "Precision", color = "Database (Median AUPR)") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Receiver Operating Characteristic Curves
p_roc <- ggplot(curve_df %>% filter(type == "ROC"), aes(x = x, y = y, group = modname, color = Database)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = db_colors) +  # Use global color mapping
  labs(title = "Receiver Operating Characteristic (ROC)", x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)", color = "Database (Median AUPR)") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8))

# Combine and save plots
curves <- p_aupr + p_prc +  p_roc
ggsave(file.path(results_dir, "precision_recall.png"), curves, width = 12, height = 3.5)

#=========================================================
# Clustering with Ground Truth
#=========================================================

# Create global mappings of run IDs to parameters
db_name_map <- setNames(col_data$db_name, rownames(col_data))
kraken_min_hit_map <- setNames(col_data$kraken2_min_hits, rownames(col_data))
bracken_threshold_map <- setNames(col_data$bracken_thresh, rownames(col_data))

#------------- Dimensionality Reduction -------------

# Plotting function for PCA, t-SNE, UMAP
plot_embedding <- function(df, xvar, yvar, title, subtitle = NULL) {
  
  # Add db_name to the embedding dataframe if not already present
  if (!"db_name" %in% colnames(df)) {
    df$db_name <- db_name_map[match(df$run_id, names(db_name_map))]
  }
  
  ggplot(df, aes_string(x = xvar, y = yvar, label = 'sub("^run_", "", run_id)')) +
    geom_point(data = subset(df, run_id != "ground_truth"), aes(color = db_name), size = 3) +
    geom_point(data = subset(df, run_id == "ground_truth"), color = "red", size = 4) +
    geom_text_repel(data = subset(df, run_id != "ground_truth"), aes(color = db_name), 
                    size = 3, max.overlaps = Inf) +
    scale_color_manual(values = db_colors) +  # Use global color mapping
    labs(title = title, subtitle = subtitle, x = xvar, y = yvar) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 9),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "none")
}

# PCA plot: projects runs based on species composition in logCPM space
pca_df <- prcomp(t(assay(se, "logcpm")), scale. = TRUE)$x %>%
  as.data.frame() %>%
  rownames_to_column("run_id")  # Add run ID as a column

pca_plot <- plot_embedding(pca_df, "PC1", "PC2", "PCA")

# t-SNE (optimized for local structure)
perplexity <- min(30, floor((ncol(assay(se, "logcpm")) - 1) / 3))  # Auto-adjust perplexity
tsne_df <- Rtsne(t(assay(se, "logcpm")), perplexity = perplexity)$Y %>%
  as.data.frame() %>%
  setNames(c("tSNE1", "tSNE2")) %>%
  mutate(run_id = colnames(se))  # Attach run IDs

tsne_plot <- plot_embedding(tsne_df, "tSNE1", "tSNE2", "t-SNE", subtitle = paste("Perplexity:", perplexity))

# UMAP (preserves global and local structure)
n_neighbors <- min(15, max(2, floor(ncol(assay(se, "logcpm")) / 2)))  # Auto-adjust neighbors
umap_cfg <- umap.defaults
umap_cfg$n_neighbors <- n_neighbors

umap_df <- umap(t(assay(se, "logcpm")), config = umap_cfg)$layout %>%
  as.data.frame() %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  mutate(run_id = colnames(assay(se, "logcpm")))  # Attach run IDs

umap_plot <- plot_embedding(umap_df, "UMAP1", "UMAP2", "UMAP", subtitle = paste("N° Neighbors:", n_neighbors))

#------------- Distance Computation -------------

# Helper: wrap pheatmap as a ggplot-compatible grob with db_name annotations
pheatmap_grob <- function(mat, show_legend = TRUE) {
  run_ids <- rownames(mat) # Extract run IDs from the matrix
  is_gt <- db_name_map[run_ids] == "ground_truth" # Identify ground truth
  annotation_row <- data.frame(Database = db_name_map[run_ids],row.names = run_ids) # Create row annotations (database name)
  annotation_col <- data.frame( # Create column annotations (k2_min_hit and b_threshold)
    KrakenMinHit = factor(ifelse(is_gt, "ground_truth", as.character(kraken_min_hit_map[run_ids]))),
    BrackenThreshold = factor(ifelse(is_gt, "ground_truth", as.character(bracken_threshold_map[run_ids]))),
    row.names = run_ids
  )
  
  # Annotation color mappings with red for GT
  ann_colors <- list(
    Database = {col <- db_colors; col["ground_truth"] <- "red"; col[names(col) %in% annotation_row$Database] },
    KrakenMinHit = {lv <- levels(annotation_col$KrakenMinHit);lv_num <- sort(as.numeric(lv[lv != "ground_truth"]))
      col <- setNames(colorRampPalette(c("white", "orange"))(length(lv_num)), as.character(lv_num))
      col["ground_truth"] <- "red";col },
    BrackenThreshold = {lv <- levels(annotation_col$BrackenThreshold); lv_num <- sort(as.numeric(lv[lv != "ground_truth"]))
      col <- setNames(colorRampPalette(c("white", "blue"))(length(lv_num)), as.character(lv_num))
      col["ground_truth"] <- "red"; col }
  )  
  
  # Generate heatmap
  p <- pheatmap(
    mat,color = colorRampPalette(c("white", "black"))(100),
    annotation_row = annotation_row, annotation_col = annotation_col,
    annotation_colors = ann_colors,
    annotation_names_row = FALSE, annotation_names_col = FALSE,
    show_rownames = FALSE, show_colnames = FALSE,
    legend = show_legend, annotation_legend = show_legend
  )
  
  as.ggplot(p[[4]])
}

# Distances in original logCPM space (species abundance profiles)
raw_dists   <- as.matrix(dist(t(assay(se, "logcpm"))))
l2_dist  <- pheatmap_grob(raw_dists, show_legend=TRUE) # Plot as heatmap 

# Distances in PCA space (first 2 principal components)
rownames(pca_df) <- pca_df$run_id
pca_dists <- pca_df[, c("PC1", "PC2")] %>%
  dist() %>% # Euclidean distance between runs
  as.matrix()
pca_dist <- pheatmap_grob(pca_dists, show_legend=FALSE) # Plot as heatmap 

# Distances in UMAP space
rownames(umap_df) <- umap_df$run_id
umap_dists  <- as.matrix(dist(umap_df[, c("UMAP1", "UMAP2")]))
umap_dist <- pheatmap_grob(umap_dists, show_legend=FALSE) # Plot as heatmap 

# Distances in t-SNE space
rownames(tsne_df) <- tsne_df$run_id
tsne_dists  <- as.matrix(dist(tsne_df[, c("tSNE1", "tSNE2")]))
tsne_dist <- pheatmap_grob(tsne_dists, show_legend=FALSE)

design <- "
AAA.
BBCC
DDEE
FFGG
"
combined_plot <- wrap_plots(list(
  A = l2_dist, B = pca_plot, C = pca_dist,
  D = tsne_plot, E = tsne_dist,
  F = umap_plot, G = umap_dist), design = design)

ggsave(file.path(results_dir, "combined_distances.png"), combined_plot, width = 14, height = 22)
# TODO: Use distances from "truth" column to rank best-matching configurations

