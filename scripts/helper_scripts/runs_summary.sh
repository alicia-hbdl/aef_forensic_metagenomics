#!/bin/bash

# Script to extract run and sample statistics from Kraken2 pipeline logs.
# Outputs a single CSV file for downstream analysis of QC, trimming, host removal, and taxonomic profiling.

# Usage: ./runs_summary.sh <path/to/runs/directory>

# To Do: 
# - Print QC stats to log file

# --- ARGUMENT PARSING ---
if [ -z "$1" ] || [ ! -d "$1" ]; then  # Check argument and directory validity
  echo "❌ Missing or invalid run directory. Usage: $0 <path/to/runs_directory>"
  exit 1
fi
RUNS_DIR=$(realpath "$1")              # Get absolute path to run directory
OUTPUT="$RUNS_DIR/runs_summary.csv"    # Set output CSV file path

# Define global variables (common to all samples) and local variables (sample-specific).
metadata_vars=(run_id db_name runtime)
global_trim_vars=(trim_clip trim_head trim_lead trim_crop trim_sliding_win trim_trail trim_min_len)
local_trim_vars=(trim_paired trim_both trim_fw_only trim_rev_only trim_dropped)
global_bt_vars=(bt_prefix bt_mode bt_sensitivity bt_mixed bt_discordant) 
local_bt_vars=(bt_paired bt_conc_0 bt_conc_1 bt_conc_more)
kraken2_vars=(kraken2_min_hits kraken2_total kraken2_classified kraken2_unclassified)
bracken_vars=(bracken_thresh bracken_species bracken_species_above_thresh bracken_species_below_thresh bracken_kept \
              bracken_discarded bracken_redistributed bracken_not_redistributed bracken_total)
# Read length statistics: raw, post-trimming, post-host removal, and for classified vs unclassified reads
read_len_vars=(preqc_avg_len_r1 preqc_med_len_r1 preqc_avg_len_r2 preqc_med_len_r2 \
              postqc_avg_len_r1 postqc_med_len_r1 postqc_avg_len_r2 postqc_med_len_r2 \
              bt_avg_r1 bt_med_r1 bt_avg_r1 bt_med_r2 \
              clas_avg_r1 clas_med_r1 clas_avg_r2 clas_med_r2 \
              unclas_avg_r1 unclas_avg_r2 unclas_med_r1 unclas_med_r2)

# Create CSV header if file not existing
if [[ ! -f "$OUTPUT" ]]; then
  {
    IFS=,
    echo "${metadata_vars[*]},sample,${global_trim_vars[*]},${local_trim_vars[*]},${global_bt_vars[*]},${local_bt_vars[*]},${kraken2_vars[*]},${bracken_vars[*]},${read_len_vars[*]}"
    unset IFS
  } > "$OUTPUT"
fi

# Locate and iterate over all pipeline log files (assumes standard format output from pipeline.sh)
LOG_FILES=$(find "$RUNS_DIR" -type f -name "*.log")
for log in $LOG_FILES; do
  # -- RUN METADATA -- 
  run_id=$(basename "$(dirname "$(dirname "$log")")") # Get run ID

  # Skip run if already processed
  if grep -q "^$run_id," "$OUTPUT"; then
      echo "Skipping $run_id (already in summary)"
      continue
  fi

  # Get database name and runtime
  db_name=$(grep "Kraken2/Bracken Database:" "$log" | awk -F': ' '{print $2}' | xargs basename)
  runtime=$(grep "Metagenomic classification completed in:" "$log" | awk -F': ' '{print $2}')    
  
  # -- KRAKEN & BRACKEN METADATA --
  # Process Kraken2/Bracken first to get unique sample IDs, as QC/trimming may not run but classification always does.
  
  # Extract Kraken2 minimum hit groups and Bracken threshold (same for all samples in run)
  kraken2_min_hits=$(grep "kraken2" "$log" | grep -oE -- '--minimum-hit-groups[ =][0-9]+' | grep -oE '[0-9]+' | head -n1)
  bracken_thresh=$(grep "Threshold:" "$log" | sed -n 's/.*Threshold: \([0-9]*\).*/\1/p' | head -n1)
  
  # Initialize arrays to store per-sample Kraken2 and Bracken statistics
  kraken_stats=()
  bracken_stats=()
  sample_ids=()
  
  # Parse each sample block from log file and extract stats
  while IFS=',' read -r sample_id minhit total classified unclassified \
                      thresh species above below kept discard redist notredist total_brack; do
    # Store extracted stats for the sample
    kraken_stats+=("$minhit,$total,$classified,$unclassified")
    bracken_stats+=("$thresh,$species,$above,$below,$kept,$discard,$redist,$notredist,$total_brack")
    sample_ids+=("$sample_id")
  done < <(
    # AWK block to extract Kraken2 and Bracken summary values per sample
    awk -v b_thresh="$bracken_thresh" -v minhits="$kraken2_min_hits" -v RS="Processing sample: " '    
    BEGIN { FS = "\n" }
    NR > 1 {
      sample = $1 # Sample ID is the first line of each block

      # Initialize all fields to placeholder values
      k_total = k_class = k_unclass = "-"
      b_species = b_above = b_below = b_kept = b_discard = b_redist = b_noredist = b_total = "-"
      
      # Loop through lines in block to extract stats
      for (i = 1; i <= NF; i++) {
        # Extract Kraken2 stats
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

        # Extract Bracken species abundance stats
        if ($i ~ /Number of species in sample:/)      { b_species = $i; gsub(/[^0-9]/, "", b_species) }
        if ($i ~ /reads > threshold:/)                { b_above = $i; gsub(/[^0-9]/, "", b_above) }
        if ($i ~ /reads < threshold:/)                { b_below = $i; gsub(/[^0-9]/, "", b_below) }
        if ($i ~ /reads kept at species level/)       { b_kept = $i; gsub(/[^0-9]/, "", b_kept) }
        if ($i ~ /reads discarded/)                   { b_discard = $i; gsub(/[^0-9]/, "", b_discard) }
        if ($i ~ /Reads distributed:/)                { b_redist = $i; gsub(/[^0-9]/, "", b_redist) }
        if ($i ~ /Reads not distributed/)             { b_noredist = $i; gsub(/[^0-9]/, "", b_noredist) }
        
        # Compute total reads considered at species level
        b_total = b_kept + b_redist
      }

      # Print CSV-formatted line if valid sample block
      if (sample != "-" && k_total != "-") {
        printf "%s,%s,%s,%s,%s,", sample, minhits, k_total, k_class, k_unclass
        printf "%s,%s,%s,%s,%s,%s,%s,%s,%s\n", b_thresh, b_species, b_above, b_below, b_kept, b_discard, b_redist, b_noredist, b_total}
      }' "$log"
  )

  read_len_stats=()
  
  # Function to extract mean and median for a given pattern
  get_read_stats() {
      local pattern="$1"
      # -m1 ensures only the first match counts (useful when raw vs trimmed match the same pattern)
      local mean=$(grep "$pattern" "$log" | grep -m1 "mean =" | awk -F'mean = ' '{print $2}' | cut -d',' -f1 | xargs)
      local median=$(grep "$pattern" "$log" | grep -m1 "mean =" | awk -F'median = ' '{print $2}' | xargs)
      echo "$mean,$median"
  }
  
  for i in "${!sample_ids[@]}"; do
    sample="${sample_ids[$i]}"
  
    if grep -q "Quality Control & Trimming: Enabled" "$log"; then  # QC/trimming enabled
      # Match any suffix between sample name and .f* filename
      pre_stats_1=$(get_read_stats "${sample}.*R1.*\.f.*: mean" "$log")
      pre_stats_2=$(get_read_stats "${sample}.*R2.*\.f.*: mean" "$log")
      post_stats_1=$(get_read_stats "${sample}_R1_paired.f*" "$log")
      post_stats_2=$(get_read_stats "${sample}_R2_paired.f*" "$log")
    else
      pre_stats_1="NA,NA"
      pre_stats_2="NA,NA"
      post_stats_1="NA,NA"
      post_stats_2="NA,NA"
    fi
  
    # Extract values per file type (metagenomic, classified, unclassified)
    bt_stats_1=$(get_read_stats "${sample}_metagenomic.1" "$log")
    bt_stats_2=$(get_read_stats "${sample}_metagenomic.2" "$log")
    clas_stats_1=$(get_read_stats "${sample}_classified_1.fastq" "$log")
    clas_stats_2=$(get_read_stats "${sample}_classified_2.fastq" "$log")
    unclas_stats_1=$(get_read_stats "${sample}_unclassified_1.fastq" "$log")
    unclas_stats_2=$(get_read_stats "${sample}_unclassified_2.fastq" "$log")
  
    # Concatenate all QC stats and append to array
    read_len_stats+=("${pre_stats_1},${pre_stats_2},${post_stats_1},${post_stats_2},${bt_stats_1},${bt_stats_2},${clas_stats_1},${clas_stats_2},${unclas_stats_1},${unclas_stats_2}")
  done
  
  # -- FASTQC & TRIMMING METADATA --
  
  # Initialize arrays to store pre-trimming QC, trimming, and post-trimming QC stats
  preqc_stats=()
  trimming_stats=()
  postqc_stats=()

  # Check if trimming was enabled for this run
  if grep -q "Quality Control & Trimming: Enabled" "$log"; then

    # Extract global Trimmomatic parameters
    # Get full argument line for TrimmomaticPE (first occurrence)
    trimmomatic_args=$(grep -m1 "TrimmomaticPE: Started with arguments:" -A1 "$log" | tail -n1)

    # Extract ILLUMINACLIP adapter file and parameters
    full_clip=$(echo "$trimmomatic_args" | sed -n 's/.*ILLUMINACLIP:\([^ ]*\).*/\1/p')
    trim_clip="$(basename "${full_clip%%:*}"):${full_clip#*:}"

    # Extract additional trimming parameters
    trim_head=$(echo "$trimmomatic_args" | sed -n 's/.*HEADCROP:\([0-9]*\).*/\1/p')
    trim_lead=$(echo "$trimmomatic_args" | sed -n 's/.*LEADING:\([0-9]*\).*/\1/p')
    trim_crop=$(echo "$trimmomatic_args" | sed -n 's/.*CROP:\([0-9]*\).*/\1/p')
    trim_sliding_win=$(echo "$trimmomatic_args" | sed -n 's/.*SLIDINGWINDOW:\([0-9]*:[0-9]*\).*/\1/p')
    trim_trail=$(echo "$trimmomatic_args" | sed -n 's/.*TRAILING:\([0-9]*\).*/\1/p')
    trim_min_len=$(echo "$trimmomatic_args" | sed -n 's/.*MINLEN:\([0-9]*\).*/\1/p')    
    
    # Loop through log to extract per-sample trimming results
    while IFS= read -r line; do
      # Capture input read pair line
      if [[ "$line" == Input\ Read\ Pairs:* ]]; then
        combined_line="$line"
      # Capture end of sample block and extract per-sample stats
      elif [[ "$line" == ✅\ Trimmed* ]]; then
        sample_id=${line#✅ Trimmed }
        sample_id=${sample_id%!}

        if [[ -n "$sample_id" ]]; then
          trim_paired=$(echo "$combined_line" | sed -n 's/.*Input Read Pairs: \([0-9]*\).*/\1/p')
          trim_both=$(echo "$combined_line" | sed -n 's/.*Both Surviving: \([0-9]*\).*/\1/p')
          trim_fw_only=$(echo "$combined_line" | sed -n 's/.*Forward Only Surviving: \([0-9]*\).*/\1/p')
          trim_rev_only=$(echo "$combined_line" | sed -n 's/.*Reverse Only Surviving: \([0-9]*\).*/\1/p')
          trim_dropped=$(echo "$combined_line" | sed -n 's/.*Dropped: \([0-9]*\).*/\1/p')

          # Store sample-specific trimming results
          trimming_stats+=("$trim_paired,$trim_both,$trim_fw_only,$trim_rev_only,$trim_dropped")
        fi
      fi
    done < "$log"
  else
    
    # Trimming not performed → set all values to NA
    # Set global trimming parameter variables to NA
    for var in "${global_trim_vars[@]}"; do
      eval "$var='NA'"
    done

    # Set per-sample trimming and QC stats to NA for each sample
    for i in "${!sample_ids[@]}"; do
      preqc_stats+=("NA,NA,NA,NA,NA,NA")
      trimming_stats+=("NA,NA,NA,NA,NA")
      postqc_stats+=("NA,NA,NA,NA,NA,NA")
    done
    
  fi

   # -- HOST-DNA REMOVAL METADATA --
  bowtie_stats=()

  if grep -q "Host DNA Removal: Enabled" "$log"; then

    # Get full Bowtie2 command (first occurrence)
    bt_args=$(grep -m1 "^bowtie2 -x" "$log")

    # Extract Bowtie2 values from command string
    bt_prefix=$(echo "$bt_args" | sed -n 's/.*-x \([^ ]*\).*/\1/p' | xargs basename)
    bt_mode=$(echo "$bt_args" | grep -q -- '--local' && echo "local" || echo "end-to-end")
    bt_sensitivity=$(echo "$bt_args" | sed -En 's/.*--(very-sensitive|sensitive|fast|very-fast).*/\1/p')
    bt_mixed=$(echo "$bt_args" | grep -q -- '--no-mixed' && echo "off" || echo "on")
    bt_discordant=$(echo "$bt_args" | grep -q -- '--no-discordant' && echo "off" || echo "on")
    
    # Loop through log lines to extract per-sample Bowtie2 alignment stats
    while read -r line; do
      [[ $line =~ ^Processing\ sample:\  ]] && sample_id=${line#*: }  && continue
      [[ $line =~ reads\; ]] && bt_paired=$(echo "$line" | sed -E 's/^([0-9]+) .*/\1/') && continue
      [[ $line =~ concordantly\ 0\ times ]] && bt_conc_0=$(echo "$line" | sed -E 's/^[[:space:]]*([0-9]+).*/\1/') && continue
      [[ $line =~ exactly\ 1\ time ]] && bt_conc_1=$(echo "$line" | sed -E 's/^[[:space:]]*([0-9]+).*/\1/') && continue
      if [[ $line =~ \>1\ times ]]; then
        bt_conc_more=$(echo "$line" | sed -E 's/^[[:space:]]*([0-9]+).*/\1/')
        # Store final sample-specific Bowtie2 stats
        bowtie_stats+=("$bt_paired,$bt_conc_0,$bt_conc_1,$bt_conc_more")
      fi        
    done < <(grep -A8 "Processing sample:" "$log")
  else 
    # Set Bowtie2 global and per-sample stats to NA if host removal was not performed
    for var in "${global_bt_vars[@]}"; do
      eval "$var='NA'"
    done

    for i in "${!sample_ids[@]}"; do
      bowtie_stats+=("NA,NA,NA,NA")
    done
  fi 

  # Write sample statistics to CSV
  for i in "${!sample_ids[@]}"; do
  {
        # Global run metadata
        for var in "${metadata_vars[@]}"; do
          printf "%s," "${!var}"
        done
    
        # Sample ID
        printf "%s," "${sample_ids[$i]}"
        
        # Global trimming parameters
        for var in "${global_trim_vars[@]}"; do
          printf "%s," "${!var}"
        done
        
        # Per-sample trimming stats
        printf "%s," "${trimming_stats[$i]}"
        
        # Global Bowtie2 parameters
        for var in "${global_bt_vars[@]}"; do
          printf "%s," "${!var}"
        done
    
        # Per-sample Bowtie2 stats
        printf "%s," "${bowtie_stats[$i]}"
    
        # Kraken2 stats
        printf "%s," "${kraken_stats[$i]}"
        
        # Bracken stats
        printf "%s," "${bracken_stats[$i]}"
        
      # Read length stats
      printf "%s\n" "${read_len_stats[$i]}"
  
  } >> "$OUTPUT"
  done
done
