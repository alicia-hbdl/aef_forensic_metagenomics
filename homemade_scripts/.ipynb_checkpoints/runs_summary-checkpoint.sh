#!/bin/bash

# Set the root comparison directory
RUNS_DIR="../zymobiomics_folder/results/runs"

# Output CSV file path
OUTPUT="$RUNS_DIR/runs_summary.csv"

# Write the CSV header
echo "Run,Sample,Database,Runtime,MinHitGroups,TotalReads,KrakenClassified,KrakenUnclassified,BrackenThreshold,NumSpecies,SpeciesAboveThres,SpeciesBelowThres,KeptReads,DiscardedReads,RedistributedReads,NotRedistributedReads" > "$OUTPUT"

# Find all kraken_pipeline.log files under comparisons/
LOG_FILES=$(find "$RUNS_DIR" -type f -name "kraken_pipeline.log")

# Process each log file
for log in $LOG_FILES; do

    # Set run name to the name of the parent directory of the log file
    RUN=$(basename "$(dirname "$(dirname "$log")")")    
    
    # Extract the database name (only basename of the path)
    DB=$(grep "Kraken2/Bracken Database Path:" "$log" | awk -F': ' '{print $2}' | xargs basename)
    
    # Extract runtime string (e.g., 00:01:45)
    RUNTIME=$(grep "Metagenomic classification completed in:" "$log" | awk -F': ' '{print $2}')
    
    # Extract the value passed to --minimum-hit-groups in Kraken2 (e.g., 2)
    MIN_HIT=$(grep -m1 "kraken2 " "$log" | grep -oE -- "--minimum-hit-groups [0-9]+" | awk '{print $2}')

    # Parse the log file using awk, block-separated by "Processing sample: "
    awk -v run="$RUN" -v db="$DB" -v runtime="$RUNTIME" -v minhit="$MIN_HIT" -v out="$OUTPUT" '
    BEGIN { FS = "\n"; RS = "Processing sample: " } # Use sample block separator
    
    NR > 1 {
        sample = $1 # First line in block is the sample name
        
        # Initialize all variables to placeholder "-"
        total = classified = unclassified = "-"
        thresh = species = above = below = kept = discard = redist = not_redist = "-"

        # Loop over each line within the block
        for (i = 1; i <= NF; i++) {
            
            # Total number of sequences processed
            if ($i ~ /sequences \(/) {
              total = $i
              gsub(/[^0-9].*$/, "", total) # Keep only the first number (strip from first non-digit onward)
            } 
            
            # Classified reads: strip percentage in parentheses and keep digits
            if ($i ~ /sequences classified \(/) { 
              gsub(/\([^)]*\)/, "", $i) # Remove parentheses and contents
              gsub(/[^0-9]/, "", $i) # Remove any non-digit characters
              classified = $i
            }
            
            # Unclassified reads
            if ($i ~ /sequences unclassified \(/) {
              gsub(/\([^)]*\)/, "", $i)
              gsub(/[^0-9]/, "", $i)
              unclassified = $i
            }
            
            # Bracken summary values
            if ($i ~ /Threshold:/)                        { thresh = $i; gsub(/[^0-9]/, "", thresh) }
            if ($i ~ /Number of species in sample:/)      { species = $i; gsub(/[^0-9]/, "", species) }
            if ($i ~ /reads > threshold:/)                { above = $i; gsub(/[^0-9]/, "", above) }
            if ($i ~ /reads < threshold:/)                { below = $i; gsub(/[^0-9]/, "", below) }
            if ($i ~ /reads kept at species level/)       { kept = $i; gsub(/[^0-9]/, "", kept) }
            if ($i ~ /reads discarded/)                   { discard = $i; gsub(/[^0-9]/, "", discard) }
            if ($i ~ /Reads distributed:/)                { redist = $i; gsub(/[^0-9]/, "", redist) }
            if ($i ~ /Reads not distributed/)             { not_redist = $i; gsub(/[^0-9]/, "", not_redist) }
        }
        
        # Print all extracted fields as a CSV row
        print run "," sample "," db "," runtime "," minhit "," total "," classified "," unclassified "," thresh "," species "," above "," below "," kept "," discard "," redist "," not_redist "," unclassified_reads >> out
    }
    ' "$log"
done

echo "✅ Summary saved to: $OUTPUT"