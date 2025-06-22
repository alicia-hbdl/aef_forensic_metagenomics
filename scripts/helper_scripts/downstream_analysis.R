#!/usr/bin/env Rscript

# Aggregates Bracken output and metadata into a SummarizedExperiment,
# then analyzes pipeline performance via read retention, classification accuracy 
# (L2, AUPR, precision/recall) and dimensionality reduction (PCA, t-SNE, UMAP),

# Rscript ./scripts/helper_scripts/downstream_analysis.R   -t ./zymobiomics_folder/raw_data/ground_truth.csv   -s ./zymobiomics_folder/results/runs/runs_summary.csv   $(find ./zymobiomics_folder/results/runs/ -type f -name "combined_breports.csv" | sort)

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
  library(gghalves)
  library(ggpp)
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
results_dir <- file.path(dirname(dirname(opts$`runs-summary`)), "benchmarking")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

#=========================================================
# Creating Summarized Experiment 
#=========================================================

## ------------- Load Bracken reports ------------- 
# Loading and merging multiple Bracken reports (already aggregated per run; no sample-level detail)

# Initialize empty table to collect species counts
combined_table <- data.frame(species = character(), stringsAsFactors = FALSE)

# Loop through each Bracken report and merge species-level counts
for (file_path in breport_paths) {
  run_id <- basename(dirname(file_path)) # Extract run ID from folder name

  # Read report and compute mean across numeric columns
  run_species <- read_csv(file_path, show_col_types = FALSE) %>%
    mutate(!!run_id := rowMeans(across(where(is.numeric), ~ replace_na(.x, 0)))) %>% # Replace NA with 0 
    select(species, !!sym(run_id)) %>% # Keep only species and new column
    mutate(species = case_when( 
      species == "Bacillus spizizenii" ~ "Bacillus subtilis",
      species == "Limosilactobacillus fermentum" ~ "Lactobacillus fermentum",
      TRUE ~ species)) %>%
    group_by(species) %>%
    summarise(!!run_id := sum(.data[[run_id]]), .groups = "drop")  # Aggregate rows for the same species name
  
  # Merge into combined table
  combined_table <- if (nrow(combined_table) == 0) {
    run_species  # First run sets the table
  } else {
    full_join(combined_table, run_species, by = "species") %>%
      mutate(across(where(is.numeric), ~ replace_na(.x, 0)))  # Fill missing values with 0
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
  arrange(run_id) %>%                                       # Sort by run_id (ascending)
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

# Filter low-quality data
se <- se %>%
  subset(, !duplicated(as.data.frame(colData(.)[, c("db_name", "kraken2_min_hits", "bracken_thresh")])))

# TODO: SOP CHECK CONSISTENCY OF DATA 

# Create global legend (auto-assign colors for all detected databases)
legend_df <- unique(colData(se)[, "db_name", drop = FALSE])
db_colors <- setNames(viridis(nrow(legend_df), option = "D"), legend_df$db_name)

# Manually assign distinct, readable colors by database group; lighter shades = smaller versions.
db_colors["ground_truth"] <- "#e6194B"
db_colors["k2_core_nt_20240904"] <- "#f58231"
db_colors["k2_eupathdb48_20230407"] <- "#ffe119"
db_colors["k2_pluspfp_20250402"] <- "#267330"
db_colors["k2_pluspfp_16gb_20250402"] <- "#3cb44b"
db_colors["k2_pluspfp_08gb_20241228"] <- "#b3e6b9"
db_colors["k2_standard_20250402"] <- "#1c3387"
db_colors["k2_standard_16gb_20250402"] <- "#4363d8"
db_colors["k2_standard_08gb_20241228"] <- "#b6c2f0"
db_colors["k2_zymobiomics_250509"] <- "#911eb4"
db_colors["k2_housepets_250510"] <- "#f032e6"
db_colors <- db_colors[order(names(db_colors))]

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
  labs(title = "Read Retention Across Stages", y = "Proportion of Raw Reads", x = NULL, color = "Database") +
  theme_minimal() +
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
  labs(title = "Read Length Progression", y = "Read Length (bp)", x = NULL, linetype = "Metric", color = "Database") +
  scale_linetype_manual(values = c("Median" = "dashed", "Mean" = "solid")) +
  scale_color_manual(values = db_colors, guide = "none") +  # Use global color mapping
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 9, face = "bold"), legend.position = "none")

## ------------------ Boxplot: Classified vs Unclassified ------------------
# Compare read lengths between classified and unclassified reads using a boxplot

# Reshape median classified/unclassified read lengths into long format for plotting
clas_unclas <- min_read_data %>%
  select(db_name, clas_med, unclas_med) %>%                            # Select relevant columns
  pivot_longer(cols = c(clas_med, unclas_med),                         # Reshape to long format
               names_to = "Type", values_to = "Reads") %>%
  mutate(Type = recode(Type,                                           # Rename for clarity
                       "clas_med" = "Classified",
                       "unclas_med" = "Unclassified")) %>%
  slice_sample(n = nrow(.))                                            # Shuffle row order

# Boxplot comparing median read lengths of classified vs unclassified reads
# Unpaired Wilcoxon test: non-parametric (handles skewed data), compares medians of independent groups (classified vs unclassified reads)
classification_length <- ggplot(clas_unclas, aes(x = Type, y = Reads, fill = Type)) +
  geom_half_violin(side = "l", color = "black", size = 0.5, trim = FALSE, alpha = 0.5) +  # Half violin on left
  geom_point(aes(color = db_name), size = 0.8, position = position_jitternudge(width = 0.2, x = 0.23, nudge.from = "jittered"))+
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white") +  # Boxplot filled white
  stat_compare_means(comparisons = list(c("Classified", "Unclassified")), method = "wilcox.test", size = 3) + # Unpaired Wilcoxon test
  stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)), size = 2.5, color = "black", vjust = -0.65) + # Show median values
  scale_fill_manual(values = c("Classified" = "#000075", "Unclassified" = "#800000"), guide = "none") +
  scale_color_manual(values = db_colors, name = "Database") +  # Use your existing db_colors
  labs(title = "Classified vs Unclassified Length", x = NULL, y = "Read Length (bp)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold")
  )

# Combine all plots (boxplot, read retention, and read length progression) and save
read_progression <- (read_retention + read_length + classification_length) + 
  plot_annotation(tag_levels = 'A', theme = theme(plot.caption = element_text(size = 8, hjust = 0, face = "italic")))
            
ggsave(file.path(results_dir, "read_progression.png"), read_progression,  width = 12, height = 4, dpi = 600)


#=========================================================
# Data Preprocessing
#=========================================================

# Prepare variables
counts_mat <- assay(se, "counts")
run_ids <- colnames(counts_mat)
non_gt_ids <- run_ids[run_ids != "ground_truth"]
db_names <- colData(se)[non_gt_ids, "db_name"]

# Normalize using classified reads (Bracken total)
class_sizes <- setNames(colData(se)$bracken_total, run_ids)
assay(se, "clas_cpm") <- t(t(counts_mat) / class_sizes * 1e6)  # Add CPM assay

# Determine pseudocount 
min <- min(assay(se, "clas_cpm")[assay(se, "clas_cpm") > 0])
min
log10(min)
max <- max(assay(se, "clas_cpm"))
max 
log10(max)
log10(0.5)

# Log-transform to stabilize variance
assay(se, "log_clas_cpm") <- log10(assay(se, "clas_cpm") + 0.5)  
metadata(se)$log_clas_cpm_dist <- as.matrix(dist(t(assay(se, "log_clas_cpm"))))
dist_clas <- metadata(se)$log_clas_cpm_dist["ground_truth", non_gt_ids]

# Normalize using total reads (bt_paired)
lib_sizes <- setNames(colData(se)$bt_paired, run_ids)
assay(se, "tot_cpm") <- t(t(counts_mat) / lib_sizes * 1e6)  # Add CPM assay
assay(se, "log_tot_cpm") <- log10(assay(se, "tot_cpm") + 0.5)  # Log-transform to stabilize variance
metadata(se)$log_tot_cpm_dist <- as.matrix(dist(t(assay(se, "log_tot_cpm"))))
dist_lib <- metadata(se)$log_tot_cpm_dist["ground_truth", non_gt_ids]

# Construct tidy dataframe
df <- tibble(
  run_id = rep(non_gt_ids, 2),
  Database = rep(db_names, 2),
  Normalization = rep(c("/Classified", "/Total"), each = length(non_gt_ids)),
  L2_Distance = c(dist_clas, dist_lib)
)

normalization <- ggplot(df, aes(x = Database, y = L2_Distance, fill = Normalization)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.4) +
  geom_jitter(aes(group = Normalization, color = Database), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1) +
  stat_compare_means(method = "wilcox.test",label = "p.format",paired = FALSE, size = 3) +
  labs(x = "Database", y = "L2 Distance (logCPM)") +
  scale_color_manual(values = db_colors) +
  scale_fill_manual(values = c("/Classified" = "#000075", "/Total" = "#800000"), name = "Normalization Method") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_blank(), legend.title = element_text(size = 9, face = "bold"), 
        legend.text = element_text(size = 8))

ggsave(file.path(results_dir, "normalization_boxplot.jpg"), normalization,  width = 10, height = 5, dpi = 600)

#=========================================================
# Database Precision and Recall 
#=========================================================

# Remove duplicated **columns** (i.e., runs) based on the log_tot_cpm assay
logcpm_mat <- assay(se, "log_tot_cpm")

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

# Boxplot of AUPR per database
p_aupr <- ggplot(aupr_df, aes(x = Database, y = aucs, fill = Database)) +
  geom_boxplot(alpha = 0.5, size = 0.2) + 
  geom_point(aes(color = Database), size = 1) +  # Correct color mapping
  scale_color_manual(values = db_colors) +
  scale_fill_manual(values = db_colors, guide = "none") +
  labs(title = "AUPR by Database", x = "Database", y = "AUPR") +
  stat_summary(fun = median, geom = "text",
               aes(label = round(after_stat(y), 2)),
               size = 3, color = "black", vjust = -0.5) +
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_blank(), legend.title = element_text(size = 9, face = "bold"), 
        legend.text = element_text(size = 8))

#------------- Plot PRC and ROC Curves -------------

# Compute median AUPR per database
median_labels <- aupr_df %>%
  group_by(Database) %>%
  summarise(Median = median(aucs, na.rm = TRUE), .groups = "drop")

# Convert curves to long-format dataframe
curve_df <- as.data.frame(sscurves) %>%
  left_join(db_df, by = "modname") 

# Define common x-axis grid
x_grid <- seq(0, 1, length.out = 100)

# Interpolate and compute summary stats
prc_summary <- curve_df %>%
  filter(type == "PRC") %>%
  group_by(Database, modname) %>%
  summarise(y_interp = list(approx(x, y, xout = x_grid, rule = 2)$y), .groups = "drop") %>%
  unnest_wider(y_interp, names_sep = "_") %>%
  pivot_longer(cols = starts_with("y_interp_"), names_prefix = "y_interp_", names_to = "x_idx", values_to = "y") %>%
  mutate(x = x_grid[as.integer(x_idx)]) %>%
  group_by(Database, x) %>%
  summarise(
    ymin = min(y, na.rm = TRUE),
    ymax = max(y, na.rm = TRUE),
    ymed = median(y, na.rm = TRUE),
    .groups = "drop"
  )

p_prc <- ggplot(prc_summary, aes(x = x, color = Database, fill = Database)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.3, color = NA) +
  geom_line(aes(y = ymed), size = 0.7) +
  scale_color_manual(values = db_colors) +
  scale_fill_manual(values = db_colors) +
  labs(title = "PR Curves", x = "Recall", y = "Precision") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Interpolate and compute summary stats for ROC
roc_summary <- curve_df %>%
  filter(type == "ROC") %>%
  group_by(Database, modname) %>%
  summarise(y_interp = list(approx(x, y, xout = x_grid, rule = 2)$y), .groups = "drop") %>%
  unnest_wider(y_interp, names_sep = "_") %>%
  pivot_longer(cols = starts_with("y_interp_"), names_prefix = "y_interp_", names_to = "x_idx", values_to = "y") %>%
  mutate(x = x_grid[as.integer(x_idx)]) %>%
  group_by(Database, x) %>%
  summarise(
    ymin = min(y, na.rm = TRUE),
    ymax = max(y, na.rm = TRUE),
    ymed = median(y, na.rm = TRUE),
    .groups = "drop"
  )

# Receiver Operating Characteristic Curves
p_roc <- ggplot(roc_summary, aes(x = x, color = Database, fill = Database)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.1, color = NA) +
  geom_line(aes(y = ymed), size = 0.7) +
  scale_color_manual(values = db_colors) +
  scale_fill_manual(values = db_colors) +
  labs(title = "ROC Curves", x = "FPR", y = "TPR", color = "Database", fill = "Database") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Combine and save plots
curves <- p_roc + p_prc + p_aupr +  
  plot_annotation(tag_levels = 'A')  
ggsave(file.path(results_dir, "precision_recall.png"), curves,  width = 10, height = 3, dpi = 600)


#=========================================================
# Clustering with Ground Truth
#=========================================================

# Create global mappings of run IDs to parameters
db_name_map <- setNames(col_data$db_name, rownames(col_data))
kraken_min_hit_map <- setNames(col_data$kraken2_min_hits, rownames(col_data))
bracken_threshold_map <- setNames(col_data$bracken_thresh, rownames(col_data))
#------------- Distance Computation -------------

# Helper: wrap pheatmap as a ggplot-compatible grob with db_name annotations
pheatmap_grob <- function(mat, show_legend = TRUE, title = NULL) {
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
    Database = {col <- db_colors; col["ground_truth"] <- "#e6194B"; col[names(col) %in% annotation_row$Database] },
    KrakenMinHit = {lv <- levels(annotation_col$KrakenMinHit);lv_num <- sort(as.numeric(lv[lv != "ground_truth"]))
      col <- setNames(colorRampPalette(c("white", "#000075"))(length(lv_num)), as.character(lv_num))
      col["ground_truth"] <- "#e6194B";col },
    BrackenThreshold = {lv <- levels(annotation_col$BrackenThreshold); lv_num <- sort(as.numeric(lv[lv != "ground_truth"]))
      col <- setNames(colorRampPalette(c("white", "#800000"))(length(lv_num)), as.character(lv_num))
      col["ground_truth"] <- "#e6194B"; col }
  )  
  
  # Generate heatmap
  p <- pheatmap(
    mat,
    color = viridis(100, option = "D"),
    annotation_row = annotation_row, annotation_col = annotation_col,
    annotation_colors = ann_colors,
    annotation_names_row = FALSE, annotation_names_col = FALSE,
    show_rownames = FALSE, show_colnames = FALSE,
    legend = show_legend, annotation_legend = show_legend, border_color = "NA"
  )
  
  as.ggplot(p[[4]]) +     
    labs(title = title, y = NULL, x = NULL) +
    theme(plot.title = element_text(size = 13, face = "bold"))
}

# -------------CPM vs log-transformed profiles -------------
# Compute Euclidean distances between CPM-normalized profiles
cpm_distance_matrix <- as.matrix(dist(t(assay(se, "tot_cpm"))))
cpm_heatmap <- pheatmap_grob(cpm_distance_matrix, show_legend = FALSE, title = "CPM-Normalized")

# Compute Euclidean distances between log-CPM-normalized profiles
logcpm_distance_matrix <- as.matrix(dist(t(logcpm_mat)))  # Uses already filtered matrix
logcpm_heatmap <- pheatmap_grob(logcpm_distance_matrix, show_legend = TRUE, title = "logCPM-Normalized")

# Combine and layout the heatmaps to compare raw vs log-transformed profiles
log_transform_comparison_plot <- wrap_plots(list(A = cpm_heatmap, B = logcpm_heatmap),design = "AABBB") + 
  plot_annotation(tag_levels = 'A')

# Save the figure to the benchmarking results directory
ggsave(filename = file.path(results_dir, "log_transform_comparison_heatmap.png"),
  plot = log_transform_comparison_plot, width = 16, height = 7, dpi = 600)

# ------------- Top 40 Profiles -------------
# Identify the 40 configurations closest to the ground truth based on Euclidean distance
closest_run_names <- logcpm_distance_matrix["ground_truth", ] %>%
  sort() %>%
  head(50) %>%
  names()

# Subset the distance matrix to include only the 35 closest runs + ground truth
top_logcpm_matrix <- logcpm_distance_matrix[closest_run_names, closest_run_names]

# Generate heatmap for the selected subset
#logcpm_heatmap <- pheatmap_grob(logcpm_distance_matrix, show_legend = TRUE, title = "All Configurations")
top_logcpm_heatmap <- pheatmap_grob(top_logcpm_matrix,show_legend = TRUE)

# Combine with full logCPM heatmap for comparison
#logcpm_comparison_plot <- wrap_plots(list(A = top_logcpm_heatmap, B = logcpm_heatmap), design = "AABBB") + 
 # plot_annotation(tag_levels = 'A')

# Save the figure
ggsave(filename = file.path(results_dir, "top_50_logcpm_comparison.jpg"),
  plot = top_logcpm_heatmap,  width = 8, height = 5.2, dpi = 600)

#------------- Dimensionality Reduction -------------
set.seed = 42 

# Plotting function for PCA, t-SNE, UMAP
plot_embedding <- function(df, xvar, yvar, title, subtitle = NULL, legend = "none") {
  
  # Add db_name to the embedding dataframe if not already present
  if (!"db_name" %in% colnames(df)) {
    df$db_name <- db_name_map[match(df$run_id, names(db_name_map))]
  }
  
  ggplot(df, aes_string(x = xvar, y = yvar, label = 'sub("^run_", "", run_id)')) +
    geom_point(data = subset(df, run_id != "ground_truth"), aes(color = db_name), size = 2) +
    geom_point(data = subset(df, run_id == "ground_truth"), color = "#e6194B", size = 3) +
    scale_color_manual(values = db_colors, name = "Database") +  # Use global color mapping
    labs(title = title, subtitle = subtitle, x = xvar, y = yvar) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 9),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = legend)
}

# Principal Component Analysis captures global variance structure
pca_df <- prcomp(t(logcpm_mat), scale. = TRUE)$x %>%
  as.data.frame() %>%
  rownames_to_column("run_id")  # Add run ID as a column

pca_plot <- plot_embedding(pca_df, "PC1", "PC2", "PCA")

# t-SNE (optimized for local structure)
perplexity <- min(30, floor((ncol(logcpm_mat) - 1) / 3))  # Dynamically set perplexity

# Remove duplicated runs
se <- se[, !duplicated(t(assay(se, "log_tot_cpm")))]  
logcpm_mat <- assay(se, "log_tot_cpm")

tsne_df <- Rtsne(t(logcpm_mat), perplexity = perplexity)$Y %>%
  as.data.frame() %>%
  setNames(c("tSNE1", "tSNE2")) %>%
  mutate(run_id = colnames(se))  # Re-attach run IDs

tsne_plot <- plot_embedding(tsne_df, "tSNE1", "tSNE2", title = "t-SNE", 
                            subtitle = paste("Perplexity:", perplexity), legend = "right")

clustering_plot <- pca_plot + tsne_plot + plot_annotation(tag_levels = 'A')

ggsave(file.path(results_dir, "pca_tsne_embedding.png"), clustering_plot,  width = 10, height = 4, dpi = 600)

