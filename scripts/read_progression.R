# Load necessary libraries
suppressPackageStartupMessages({
  library(tidyverse)    # for data analysis, manipulation, and visualization (loads dplyr, tidyr, ggplot2, readr, stringr, etc.)
  library(scales)       # for ggplot2 aesthetics
  library(patchwork)    # for combining multiple plots
})

# -- LOAD AND CLEAN DATA --

# Path to input file
filepath <- "/Users/aliciahobdell/Desktop/final_project/zymobiomics_folder/results/runs/new_runs_summary.csv"

# Load main data
runs_summary <- read.csv(filepath, header = TRUE, stringsAsFactors = FALSE)

# -- PLOT 1: READ NUMBER PROGRESSION --

# Define relevant columns
total_cols <- c("sample", "trim_paired", "bt_paired", "kraken2_total", "kraken2_classified", "bracken_total")
# Define stage names
stage_names <- setdiff(total_cols, "sample")

head(runs_summary)

# Aggregate by sample
aggregated_summary <- runs_summary %>%
  select(all_of(total_cols)) %>%
  group_by(sample) %>%  # Grouping by 'sample'
  summarise(across(all_of(stage_names), function(x) mean(x, na.rm = TRUE)), .groups = "drop")


# Reshape data for plotting
long_summary <- aggregated_summary %>%
  pivot_longer(cols = all_of(stage_names), names_to = "Original", values_to = "Reads") %>%
  mutate(Stage = factor(Original, levels = stage_names)) %>%  # Correctly assign Stage
  group_by(sample = sample) %>%
  mutate(Fraction = Reads / Reads[Original == "trim_paired"]) %>%
  ungroup()

# Compute average fraction per stage
mean_summary <- long_summary %>%
  group_by(Stage) %>%
  summarise(Fraction = mean(Fraction, na.rm = TRUE), .groups = "drop")

# Plot read retention
p1 <- ggplot(long_summary, aes(x = Stage, y = Fraction, group = sample, color = sample)) +
  geom_line(linewidth = 0.5, linetype = "dashed") +                     # sample lines
  geom_line(data = mean_summary, aes(x = Stage, y = Fraction, group = 1),
            color = "black", linewidth = 1, inherit.aes = FALSE) +      # Mean line
  geom_point(data = mean_summary, aes(x = Stage, y = Fraction), 
             color = "black", size = 2, inherit.aes = FALSE) +     # Mean points
  geom_text(data = mean_summary, 
            aes(x = Stage, y = Fraction, label = scales::percent(Fraction, accuracy = 0.01)),
            color = "black", size = 3.5, vjust = -0.6, hjust = -0.3, inherit.aes = FALSE) +
  theme_bw() +                                                     # White background
  scale_y_continuous(labels = percent_format()) +                 # Y-axis in %
  labs(title = "Proportion of Reads Retained at Each Pipeline Stage",
       y = "Proportion of Raw Reads") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))         # Rotated x labels

# Save plot
ggsave("/Users/aliciahobdell/Desktop/final_project/zymobiomics_folder/results/read_number_progression.png", 
       width = 8, height = 6, dpi = 300)

# -- PLOT 2: READ LENGTH PROGRESSION -- 
# WORK ON THIS ONCE WE HAVE ALL THE DATA FROM BRACKEN ETC. 

# Define columns for read lengths
bp_cols <- c("sample", 
             "preqc_avg_len_r1","preqc_avg_len_r2", "preqc_med_len_r1","preqc_med_len_r2",
             "postqc_avg_len_r1","postqc_avg_len_r2","postqc_med_len_r1","postqc_med_len_r2",
             "bt_avg_len","bt_med_len",
             "kraken2_avg_len","kraken2_med_len")


# Compute per-sample read length summaries 
aggregated_summary <- runs_summary %>%
  select(all_of(bp_cols)) %>%
  mutate(
    preqc_avg_len = pmin(preqc_avg_len_r1, preqc_avg_len_r2, na.rm = TRUE),
    preqc_med_len = pmin(preqc_med_len_r1, preqc_med_len_r2, na.rm = TRUE),
    postqc_avg_len = pmin(postqc_avg_len_r1, postqc_avg_len_r2, na.rm = TRUE),
    postqc_med_len = pmin(postqc_med_len_r1, postqc_med_len_r2, na.rm = TRUE)
  ) %>%
  group_by(sample) %>%
  summarise(across(
    c(preqc_avg_len, preqc_med_len,postqc_avg_len, postqc_med_len,
      bt_avg_len, bt_med_len, kraken2_avg_len, kraken2_med_len),
    ~mean(.x, na.rm = TRUE)
  ), .groups = "drop")

head(aggregated_summary)
# Stage mapping for plotting
stage_names <- c(
  "FastQC_1" = "RawReads",
  "FastQC_2" = "TrimmedReads",
  "Bowtie2"  = "HostFilteredReads",
  "Kraken2"  = "ClassifiedReads"
)

# Reshape read length data
length_long <- aggregated_summary %>%
  pivot_longer(cols = -sample, names_to = "Original", values_to = "Length") %>%
  separate(Original, into = c("Tool", "MetricType"), sep = "_(?=[^_]+$)", remove = FALSE) %>%
  mutate(
    Stage  = factor(stage_names[Tool], levels = stage_names),
    Metric = if_else(str_detect(MetricType, "Avg"), "Mean", "Median")
  ) %>%
  select(sample, Stage, Metric, Length)

# Compute mean read length per stage and metric
mean_length_summary <- length_long %>%
  group_by(Stage, Metric) %>%
  summarise(Length = mean(Length, na.rm = TRUE), .groups = "drop")

# Plot read length progression
p2 <- ggplot(length_long, aes(x = Stage, y = Length, group = interaction(sample, Metric),
                        color = sample, linetype = Metric)) +
  geom_line(linewidth = 0.5, na.rm = TRUE) +                               # sample lines
  geom_point(size = 2) +
  geom_line(data = mean_length_summary, aes(x = Stage, y = Length, group = Metric),
            color = "black", linewidth = 1, inherit.aes = FALSE) +         # Mean lines
  geom_point(data = mean_length_summary, aes(x = Stage, y = Length),
             color = "black", size = 2, inherit.aes = FALSE) +             # Mean points
  geom_text(data = mean_length_summary,
            aes(x = Stage, y = Length, label = round(Length, 1)),
            color = "black", size = 3.5, vjust = -0.6, hjust = -0.3, inherit.aes = FALSE) +
  theme_bw() +                                                             # White background
  labs(title = "Read Length Across Pipeline Stages", 
       y = "Read Length (bp)", linetype = "Metric") +
  scale_linetype_manual(values = c("Median" = "dashed", "Mean" = "dotted")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plot
ggsave("/Users/aliciahobdell/Desktop/final_project/zymobiomics_folder/results/read_length_progression.png", 
       width = 8, height = 6, dpi = 300)

# Combine plots side-by-side
combined_plot <- p1 / p2 + plot_layout(ncol = 1)
# Save plot
ggsave("/Users/aliciahobdell/Desktop/final_project/zymobiomics_folder/results/combined_progression.png", 
       width = 8, height = 6, dpi = 300)
