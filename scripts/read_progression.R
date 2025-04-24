library(tidyverse)
library(scales)
library(patchwork)

# -- LOAD AND CLEAN DATA --

# Path to input file
filepath <- "/Users/aliciahobdell/Desktop/final_project/zymobiomics_folder/results/runs/runs_summary.csv"

# Read header rows
headers <- read.csv(filepath, nrows = 2, header = FALSE, stringsAsFactors = FALSE)

# Combine multi-level headers
combined_headers <- paste(headers[1, ], headers[2, ], sep = "_")

# Load main data
runs_summary <- read.csv(filepath, skip = 2, header = FALSE, stringsAsFactors = FALSE)

# Apply cleaned column names
colnames(runs_summary) <- make.names(combined_headers, unique = TRUE)

# -- PLOT 1: READ NUMBER PROGRESSION --

# Define relevant columns
total_cols <- c("Metadata_Sample", "Trimmomatic_InputReadPairs", "Bowtie2_TotReads", 
                "Kraken2_TotalReads", "Kraken2_KrakenClassified", "Bracken_BrackenTot")

# Aggregate by sample
aggregated_summary <- runs_summary %>%
  select(all_of(total_cols)) %>%
  group_by(Sample = Metadata_Sample) %>%
  summarise(across(-Metadata_Sample, \(x) mean(x, na.rm = TRUE)), .groups = "drop")

# Define stage names
stage_names <- c(
  "Trimmomatic_InputReadPairs" = "RawReads",
  "Bowtie2_TotReads"           = "TrimmedReads",
  "Kraken2_TotalReads"         = "HostFilteredReads",
  "Kraken2_KrakenClassified"   = "ClassifiedReads", 
  "Bracken_BrackenTot"         = "QuantifiedReads"
)

# Reshape data for plotting
long_summary <- aggregated_summary %>%
  pivot_longer(cols = names(stage_names), names_to = "Original", values_to = "Reads") %>%
  mutate(Stage = factor(stage_names[Original], levels = stage_names)) %>%
  group_by(Sample) %>%
  mutate(Fraction = Reads / Reads[Original == "Trimmomatic_InputReadPairs"]) %>%
  ungroup()

# Compute average fraction per stage
mean_summary <- long_summary %>%
  group_by(Stage) %>%
  summarise(Fraction = mean(Fraction, na.rm = TRUE), .groups = "drop")

# Plot read retention
p1 <- ggplot(long_summary, aes(x = Stage, y = Fraction, group = Sample, color = Sample)) +
  geom_line(linewidth = 0.5, linetype = "dashed") +                     # Sample lines
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

# Define columns for read lengths
bp_cols <- c("Metadata_Sample", "FastQC_1_MedReadLength_1", "FastQC_1_MedReadLength_2", 
             "FastQC_1_AvgReadLength_1", "FastQC_1_AvgReadLength_2", 
             "FastQC_2_MedReadLength_1", "FastQC_2_MedReadLength_2", 
             "FastQC_2_AvgReadLength_1", "FastQC_2_AvgReadLength_2", 
             "Bowtie2_MedReadLen", "Bowtie2_AvgReadLen",
             "Kraken2_MedReadLen", "Kraken2_AvgReadLen")


# Compute per-sample read length summaries 
aggregated_summary <- runs_summary %>%
  select(all_of(bp_cols)) %>%
  rename(Sample = Metadata_Sample) %>%  # 🔁 Rename before any grouping/select
  mutate(
    FastQC_1_MedReadLength = pmin(FastQC_1_MedReadLength_1, FastQC_1_MedReadLength_2, na.rm = TRUE),
    FastQC_1_AvgReadLength = pmin(FastQC_1_AvgReadLength_1, FastQC_1_AvgReadLength_2, na.rm = TRUE),
    FastQC_2_MedReadLength = pmin(FastQC_2_MedReadLength_1, FastQC_2_MedReadLength_2, na.rm = TRUE),
    FastQC_2_AvgReadLength = pmin(FastQC_2_AvgReadLength_1, FastQC_2_AvgReadLength_2, na.rm = TRUE)
  ) %>%
  group_by(Sample) %>%
  summarise(across(
    c(FastQC_1_MedReadLength, FastQC_1_AvgReadLength,FastQC_2_MedReadLength, FastQC_2_AvgReadLength,
      Bowtie2_MedReadLen, Bowtie2_AvgReadLen, Kraken2_MedReadLen, Kraken2_AvgReadLen),
    ~mean(.x, na.rm = TRUE)
  ), .groups = "drop")

# Stage mapping for plotting
stage_names <- c(
  "FastQC_1" = "RawReads",
  "FastQC_2" = "TrimmedReads",
  "Bowtie2"  = "HostFilteredReads",
  "Kraken2"  = "ClassifiedReads"
)

# Reshape read length data
length_long <- aggregated_summary %>%
  pivot_longer(cols = -Sample, names_to = "Original", values_to = "Length") %>%
  separate(Original, into = c("Tool", "MetricType"), sep = "_(?=[^_]+$)", remove = FALSE) %>%
  mutate(
    Stage  = factor(stage_names[Tool], levels = stage_names),
    Metric = if_else(str_detect(MetricType, "Avg"), "Mean", "Median")
  ) %>%
  select(Sample, Stage, Metric, Length)

# Compute mean read length per stage and metric
mean_length_summary <- length_long %>%
  group_by(Stage, Metric) %>%
  summarise(Length = mean(Length, na.rm = TRUE), .groups = "drop")

# Plot read length progression
p2 <- ggplot(length_long, aes(x = Stage, y = Length, group = interaction(Sample, Metric),
                        color = Sample, linetype = Metric)) +
  geom_line(linewidth = 0.5, na.rm = TRUE) +                               # Sample lines
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
