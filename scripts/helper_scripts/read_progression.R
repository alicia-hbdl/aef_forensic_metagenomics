#!/usr/bin/env Rscript

# This script generates two plots to visualize read number and read length progression across pipeline stages from a runs summary CSV.

# Usage: Rscript read_progression.R <path/to/runs_summary.csv># Load necessary libraries

# Load necessary libraries
suppressPackageStartupMessages({
  library(tidyverse)    # for data analysis, manipulation, and visualization (loads dplyr, tidyr, ggplot2, readr, stringr, etc.)
  library(scales)       # for ggplot2 aesthetics
  library(patchwork)    # for combining multiple plots
})

# -- LOAD AND CLEAN DATA --

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript read_progression.R <path/to/runs_summary.csv>")
file_path <- args[1]
if (!file.exists(file_path)) stop("Error: Run summary file not found!")

runs_summary <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE) 

# -- PLOT 1: READ NUMBER PROGRESSION --

# Define relevant columns
total_cols <- c("sample", "trim_paired", "bt_paired", "kraken2_total", "kraken2_classified", "bracken_total")

# Define stage names
stage_names <- setdiff(total_cols, "sample")

# Load and process the data
aggregated_summary <- runs_summary %>%
  select(all_of(total_cols)) %>%
  mutate(across(all_of(stage_names), as.numeric)) %>%  # Convert to numeric
  group_by(sample) %>%
  summarise(across(all_of(stage_names), \(x) mean(x, na.rm = TRUE)), .groups = "drop")  

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
output_path <- file.path(dirname(file_path), "read_number_progression.png")
ggsave(output_path, width = 8, height = 6, dpi = 300)

# -- PLOT 2: READ LENGTH PROGRESSION -- 
# WORK ON THIS ONCE WE HAVE ALL THE DATA FROM BRACKEN ETC. 

# Define relevant columns
bp_cols <- c(
  "sample",
  "preqc_avg_len_r1", "preqc_avg_len_r2", "preqc_med_len_r1", "preqc_med_len_r2",
  "postqc_avg_len_r1", "postqc_avg_len_r2", "postqc_med_len_r1", "postqc_med_len_r2",
  "bt_avg_r1", "bt_avg_r2", "bt_med_r1", "bt_med_r2",
  "clas_avg_r1", "clas_avg_r2", "clas_med_r1", "clas_med_r2"
)

# Summarize mean and median read lengths per stage
length_summary <- runs_summary %>%
  select(all_of(bp_cols)) %>%
  mutate(across(-sample, as.numeric)) %>%  # Ensure all length columns are numeric
  mutate(
    preqc_avg = pmin(preqc_avg_len_r1, preqc_avg_len_r2, na.rm = TRUE),
    preqc_med = pmin(preqc_med_len_r1, preqc_med_len_r2, na.rm = TRUE),
    postqc_avg = pmin(postqc_avg_len_r1, postqc_avg_len_r2, na.rm = TRUE),
    postqc_med = pmin(postqc_med_len_r1, postqc_med_len_r2, na.rm = TRUE),
    bt_avg = pmin(bt_avg_r1, bt_avg_r2, na.rm = TRUE),
    bt_med = pmin(bt_med_r1, bt_med_r2, na.rm = TRUE),
    clas_avg = pmin(clas_avg_r1, clas_avg_r2, na.rm = TRUE),
    clas_med = pmin(clas_med_r1, clas_med_r2, na.rm = TRUE),
  ) %>%
  select(sample, preqc_avg, preqc_med, postqc_avg, postqc_med,
         bt_avg, bt_med, clas_avg, clas_med)

# Reshape to long format
length_long <- length_summary %>%
  pivot_longer(cols = -sample, names_to = "Stage_Metric", values_to = "Length") %>%
  separate(Stage_Metric, into = c("Stage", "Metric"), sep = "_") %>%
  mutate(
    Stage = factor(Stage,
                   levels = c("preqc", "postqc", "bt", "clas"),
                   labels = c("Raw Reads", "Trimmed Reads", "Host Filtered", "Classified Reads")),
    Metric = if_else(Metric == "avg", "Mean", "Median")
  )

# Mean summary
mean_length_summary <- length_long %>%
  group_by(Stage, Metric) %>%
  summarise(Length = mean(Length, na.rm = TRUE), .groups = "drop")

print(mean_length_summary)

# Plot
p2 <- ggplot(length_long, aes(x = Stage, y = Length, group = interaction(sample, Metric),
                              color = sample, linetype = Metric)) +
  geom_line(linewidth = 0.6, alpha = 0.6, na.rm = TRUE) +
  geom_point(size = 1.8, na.rm = TRUE) +
  geom_line(data = mean_length_summary, aes(x = Stage, y = Length, group = Metric),
            color = "black", linewidth = 1, inherit.aes = FALSE) +
  geom_point(data = mean_length_summary, aes(x = Stage, y = Length),
             color = "black", size = 2, inherit.aes = FALSE) +
  geom_text(data = mean_length_summary,
            aes(x = Stage, y = Length, label = round(Length, 1)), color = "black", size = 3.2, vjust = -0.6,
            inherit.aes = FALSE) +
  labs(title = "Read Length Across Pipeline Stages", y = "Read Length (bp)", x = NULL, linetype = "Metric") +
  scale_linetype_manual(values = c("Median" = "dashed", "Mean" = "solid")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# Save
output_path <- file.path(dirname(file_path), "read_length_progression.png")
ggsave(output_path, p2, width = 8, height = 6, dpi = 300)

# Save combined format 
combined_plot <- p1 + p2 

output_path <- file.path(dirname(file_path), "combined_progression.png")
ggsave(output_path, combined_plot,  width = 16, height = 8, dpi = 300)

