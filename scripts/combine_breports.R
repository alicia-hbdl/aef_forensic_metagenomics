# Load necessary library
library(tidyverse)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript combine_breports.R <path/to/breports>") # # At least one argument (file paths) is provided
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
output_dir <- dirname(dirname(file_paths[1]))
output_file <- file.path(output_dir, "combined_breports.csv")
write_csv(combined_table, output_file)
print(combined_table)