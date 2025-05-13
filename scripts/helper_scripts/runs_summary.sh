#!/bin/bash

# This script processes Kraken2 pipeline log files, extracting relevant statistics for each sample, including metadata, quality control, trimming, Bowtie2 alignment, Kraken2 classification, and Bracken abundance estimation, and saves the data into a CSV file for further analysis.

# Usage: ./runs_summary.sh <path/to/runs/directory>

# To Do: 
# - Print QC stats to log file
# - Get length of Kraken classified and unclassified reads
# - Print Bowtie2 command to log file

# --- Argument Parsing ---
# Ensure that the `rundir` argument is provided
if [ -z "$1" ]; then
  echo "❌ Error: Missing 'rundir' argument. Usage: $0 <path/to/runs/directory>"
  exit 1
fi

# Set the run directory to the provided argument
RUNS_DIR=$(realpath "$1")
OUTPUT="$RUNS_DIR/runs_summary.csv"

# Arrays to store metadata and other variables
metadata_vars=(run_id db_name runtime)
preqc_vars=(preqc_pct_dups_r1 preqc_pct_dups_r2 preqc_pct_gc_r1 preqc_pct_gc_r2 preqc_avg_len_r1 preqc_avg_len_r2 preqc_med_len_r1 preqc_med_len_r2 preqc_fail_r1 preqc_fail_r2)
global_trim_vars=(trim_clip trim_head trim_lead trim_crop trim_sliding_win trim_trail trim_min_len)
local_trim_vars=(trim_paired trim_both trim_fw_only trim_rev_only trim_dropped)
postqc_vars=(postqc_pct_dups_r1 postqc_pct_dups_r2 postqc_pct_gc_r1 postqc_pct_gc_r2 postqc_avg_len_r1 postqc_avg_len_r2 postqc_med_len_r1 postqc_med_len_r2 postqc_fail_r1 postqc_fail_r2)
global_bt_vars=(bt_prefix bt_mode bt_sensitivity bt_mixed bt_discordant) 
local_bt_vars=(bt_paired bt_conc_0 bt_conc_1 bt_conc_more bt_avg_len bt_med_len)
kraken2_vars=(kraken2_min_hits kraken2_total kraken2_classified kraken2_unclassified kraken2_avg_len kraken2_med_len)
bracken_vars=(bracken_thresh bracken_species bracken_species_above_thresh bracken_species_below_thresh bracken_kept bracken_discarded bracken_redistributed bracken_not_redistributed bracken_total)

# Create CSV header if file doesn't exist
if [[ ! -f "$OUTPUT" ]]; then
  {
    IFS=,
    echo "${metadata_vars[*]},sample,${preqc_vars[*]},${global_trim_vars[*]},${local_trim_vars[*]},${postqc_vars[*]},${global_bt_vars[*]},${local_bt_vars[*]},${kraken2_vars[*]},${bracken_vars[*]}"
    unset IFS
  } > "$OUTPUT"
fi

# Find all kraken_pipeline.log files
LOG_FILES=$(find "$RUNS_DIR" -type f -name "*.log")

# Process each log file
for log in $LOG_FILES; do
  # -- RUN METADATA -- 
  run_id=$(basename "$(dirname "$(dirname "$log")")") # Get run ID

  # Skip if run already exists in the output
  if grep -q "^$run_id," "$OUTPUT"; then
      echo "Skipping $run_id (already in summary)"
      continue
  fi

  # Get database and runtime info
  db_name=$(grep "Kraken2/Bracken Database:" "$log" | awk -F': ' '{print $2}' | xargs basename)
  runtime=$(grep "Metagenomic classification completed in:" "$log" | awk -F': ' '{print $2}')    
  
  # -- KRAKEN & BRACKEN METADATA --
  # Process Kraken and Bracken metadata first to capture sample IDs, as this step is always executed.
  kraken2_min_hits=$(grep "kraken2" "$log" | grep -oE -- '--minimum-hit-groups[ =][0-9]+' | grep -oE '[0-9]+' | head -n1)
  bracken_thresh=$(grep "Threshold:" "$log" | sed -n 's/.*Threshold: \([0-9]*\).*/\1/p' | head -n1)
  
  # Parse Kraken2 & Bracken summary blocks and store stats for each sample
  kraken_stats=()
  bracken_stats=()
  sample_ids=()
  while IFS=',' read -r sample_id minhit total classified unclassified avg med \
                      thresh species above below kept discard redist notredist total_brack; do
    # Store sample-specific stats in associative arrays (corrected quotes)
    kraken_stats+=("$minhit,$total,$classified,$unclassified,$avg,$med")
    bracken_stats+=("$thresh,$species,$above,$below,$kept,$discard,$redist,$notredist,$total_brack")
    sample_ids+=("$sample_id")
  done < <(
    awk -v b_thresh="$bracken_thresh" -v minhits="$kraken2_min_hits" -v RS="Processing sample: " '    
    BEGIN { FS = "\n" }
    NR > 1 {
      sample = $1 # Sample ID is first line of each block

      # Initialize output fields with placeholders
      k_total = k_class = k_unclass = k_avg = k_med = "-"
      b_species = b_above = b_below = b_kept = b_discard = b_redist = b_noredist = b_total = "-"
      
      # Loop over each line in the block to extract relevant stats
      # Extract Kraken2 classification summary values
      for (i = 1; i <= NF; i++) {
        # Kraken2 values
        if ($i ~ /sequences \(/) {
          k_total = $i
          gsub(/[^0-9].*$/, "", k_total) # Keep digits before first non-digit
        } 
        if ($i ~ /sequences classified \(/) {
          k_class = $i
          gsub(/\([^)]*\)/, "", k_class) # Remove percentage in parentheses
          gsub(/[^0-9]/, "", k_class) # Keep only digits
        }
        if ($i ~ /sequences unclassified \(/) {
          k_unclass = $i
          gsub(/\([^)]*\)/, "", k_unclass)
          gsub(/[^0-9]/, "", k_unclass)
        }

        # Extract Bracken species abundance summary values
        if ($i ~ /Number of species in sample:/)      { b_species = $i; gsub(/[^0-9]/, "", b_species) }
        if ($i ~ /reads > threshold:/)                { b_above = $i; gsub(/[^0-9]/, "", b_above) }
        if ($i ~ /reads < threshold:/)                { b_below = $i; gsub(/[^0-9]/, "", b_below) }
        if ($i ~ /reads kept at species level/)       { b_kept = $i; gsub(/[^0-9]/, "", b_kept) }
        if ($i ~ /reads discarded/)                   { b_discard = $i; gsub(/[^0-9]/, "", b_discard) }
        if ($i ~ /Reads distributed:/)                { b_redist = $i; gsub(/[^0-9]/, "", b_redist) }
        if ($i ~ /Reads not distributed/)             { b_noredist = $i; gsub(/[^0-9]/, "", b_noredist) }
        
        # Compute total reads considered at species level
        b_total = b_kept + b_redist

        k_avg= "TODO"
        k_med= "TODO"
      }

      # Output all values in a single CSV-formatted line
      if (sample != "-" && kraken2_total != "-") {
      printf "%s,%s,%s,%s,%s,%s,%s,", sample, minhits, k_total, k_class, k_unclass, k_avg, k_med
      printf "%s,%s,%s,%s,%s,%s,%s,%s,%s\n", b_thresh, b_species, b_above, b_below, b_kept, b_discard, b_redist, b_noredist, b_total}
    }' "$log"
  )

  # -- FASTQC & TRIMMING METADATA --
  preqc_stats=()
  trimming_stats=()
  postqc_stats=()

  if grep -q "Quality Control & Trimming: Enabled" "$log"; then
    
    for i in "${!sample_ids[@]}"; do
      preqc_stats+=("TODO,TODO,TODO,TODO,TODO,TODO,TODO,TODO,TODO,TODO")
      postqc_stats+=("TODO,TODO,TODO,TODO,TODO,TODO,TODO,TODO,TODO,TODO")
    done

    # Extract Trimmomatic arguments (from first occurrence)
    trimmomatic_args=$(grep -m1 "TrimmomaticPE: Started with arguments:" -A1 "$log" | tail -n1)

    # ILLUMINACLIP: basename + parameters
    full_clip=$(echo "$trimmomatic_args" | sed -n 's/.*ILLUMINACLIP:\([^ ]*\).*/\1/p')
    trim_clip="$(basename "${full_clip%%:*}"):${full_clip#*:}"

    # Other parameters
    trim_head=$(echo "$trimmomatic_args" | sed -n 's/.*HEADCROP:\([0-9]*\).*/\1/p')
    trim_lead=$(echo "$trimmomatic_args" | sed -n 's/.*LEADING:\([0-9]*\).*/\1/p')
    trim_crop=$(echo "$trimmomatic_args" | sed -n 's/.*CROP:\([0-9]*\).*/\1/p')
    trim_sliding_win=$(echo "$trimmomatic_args" | sed -n 's/.*SLIDINGWINDOW:\([0-9]*:[0-9]*\).*/\1/p')
    trim_trail=$(echo "$trimmomatic_args" | sed -n 's/.*TRAILING:\([0-9]*\).*/\1/p')
    trim_min_len=$(echo "$trimmomatic_args" | sed -n 's/.*MINLEN:\([0-9]*\).*/\1/p')    
    
    while IFS= read -r line; do
      if [[ "$line" == Input\ Read\ Pairs:* ]]; then
        combined_line="$line"
      elif [[ "$line" == ✅\ Trimmed* ]]; then
        sample_id=${line#✅ Trimmed }
        sample_id=${sample_id%!}

        if [[ -n "$sample_id" ]]; then
          trim_paired=$(echo "$combined_line" | sed -n 's/.*Input Read Pairs: \([0-9]*\).*/\1/p')
          trim_both=$(echo "$combined_line" | sed -n 's/.*Both Surviving: \([0-9]*\).*/\1/p')
          trim_fw_only=$(echo "$combined_line" | sed -n 's/.*Forward Only Surviving: \([0-9]*\).*/\1/p')
          trim_rev_only=$(echo "$combined_line" | sed -n 's/.*Reverse Only Surviving: \([0-9]*\).*/\1/p')
          trim_dropped=$(echo "$combined_line" | sed -n 's/.*Dropped: \([0-9]*\).*/\1/p')

          trimming_stats+=("$trim_paired,$trim_both,$trim_fw_only,$trim_rev_only,$trim_dropped")
        fi
      fi
    done < "$log"
  else
    
    # Set global trimming variables to "NA" if trimming was not performed
    for var in "${global_trim_vars[@]}"; do
      eval "$var='NA'"
    done

    # Set local qc and trimming variables to "NA" for each sample
    for i in "${!sample_ids[@]}"; do
      preqc_stats+=("NA,NA,NA,NA,NA,NA,NA,NA,NA,NA")
      trimming_stats+=("NA,NA,NA,NA,NA")
      postqc_stats+=("NA,NA,NA,NA,NA,NA,NA,NA,NA,NA")
    done
    
  fi

   # -- HOST-DNA REMOVAL METADATA --
  bowtie_stats=()

  if grep -q "Host DNA Removal: Enabled" "$log"; then

    # Extract Bowtie2 arguments (from first occurrence)
    bt_prefix="TODO"
    bt_mode="TODO"
    bt_sensitivity="TODO"
    bt_mixed="TODO"
    bt_discordant="TODO"

    # Parse Bowtie2 stats line-by-line for each sample
    while read -r line; do
    # grep -A7 "Processing sample:" "$log" | while read -r line; do
      [[ $line =~ ^Processing\ sample:\  ]] && sample_id=${line#*: }  && continue
      [[ $line =~ reads\; ]] && bt_paired=$(echo "$line" | sed -E 's/^([0-9]+) .*/\1/') && continue
      [[ $line =~ concordantly\ 0\ times ]] && bt_conc_0=$(echo "$line" | sed -E 's/^[[:space:]]*([0-9]+).*/\1/') && continue
      [[ $line =~ exactly\ 1\ time ]] && bt_conc_1=$(echo "$line" | sed -E 's/^[[:space:]]*([0-9]+).*/\1/') && continue
      if [[ $line =~ \>1\ times ]]; then
        bt_conc_more=$(echo "$line" | sed -E 's/^[[:space:]]*([0-9]+).*/\1/')
        bt_avg_len="TODO"
        bt_med_len="TODO"
      # Store stats in sample-specific variable once full block is complete       
      bowtie_stats+=("$bt_paired,$bt_conc_0,$bt_conc_1,$bt_conc_more,$bt_avg_len,$bt_med_len")
      fi        
    done < <(grep -A7 "Processing sample:" "$log")
    # done
  else 

    for var in "${global_bt_vars[@]}"; do
      eval "$var='NA'"
    done

    for i in "${!sample_ids[@]}"; do
      bowtie_stats+=("NA,NA,NA,NA,NA,NA")
    done
  fi 

# Assign values in the sheet 
for i in "${!sample_ids[@]}"; do
{
    # Print metadata (global per-run)
    for var in "${metadata_vars[@]}"; do
      printf "%s," "${!var}"
    done

    # Print sample ID
    printf "%s," "${sample_ids[$i]}"

    printf "%s," "${preqc_stats[$i]}"

    # Print other per-sample global and local stats
    for var in "${global_trim_vars[@]}"; do
      printf "%s," "${!var}"
    done
    
    printf "%s," "${trimming_stats[$i]}"
    
    printf "%s," "${postqc_stats[$i]}"

    for var in "${global_bt_vars[@]}"; do
      printf "%s," "${!var}"
    done

    printf "%s," "${bowtie_stats[$i]}"
    printf "%s," "${kraken_stats[$i]}"
    printf "%s\n" "${bracken_stats[$i]}"

  } >> "$OUTPUT"
  done 
done