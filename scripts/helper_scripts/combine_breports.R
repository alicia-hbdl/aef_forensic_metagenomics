#!/usr/bin/env Rscript

# Combines multiple Bracken TSV reports into a single CSV, filtering for species-level data.

# Usage: Rscript combine_breports.R <path/to/breports>

# Load necessary libraries
suppressPackageStartupMessages({
  library(tidyverse)    # for data analysis, manipulation, and visualization (loads dplyr, tidyr, ggplot2, readr, stringr, etc.)
})


# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript combine_breports.R <path/to/breports>") # At least one argument (file paths) is provided
file_paths <- args

for (file_path in file_paths) { # Check if all provided file paths exist
  if (!file.exists(file_path)) stop(paste("Error: File not found:", file_path))
}

# Define the expected column names 
headers <- c('ReadFraction', 'TotalReads', 'KrakenReads', 'LevelID', 'TaxID', 'species')

# Initialize an empty data frame to store the combined results
combined_table <- data.frame(species = character(), stringsAsFactors = FALSE)

# Loop over each file path and process the contents
for (file_path in file_paths) {
  breport <- read_tsv(file_path, col_names = headers, show_col_types = FALSE) %>%
    filter(LevelID == 'S') %>% # Filter rows where 'LevelID' is 'S'
    mutate(species = trimws(species)) %>%  # Trim whitespaces from species names
    select(species, TotalReads) %>%
    rename(!!gsub("\\.breport$", "", basename(file_path)) := TotalReads)  # Rename 'TotalReads' dynamically based on the file name
  
  # Merge the current file's data with the existing combined table by 'species'
  if (nrow(combined_table) == 0) {
    combined_table <- breport
  } else {
    combined_table <- left_join(combined_table, breport, by = "species")
  }
}

# Write the combined data table to a CSV file in the results directory
write_csv(combined_table, file.path(dirname(dirname(file_paths[1])), "combined_breports.csv"))