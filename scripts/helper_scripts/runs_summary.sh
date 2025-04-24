#!/bin/bash

# Set the root comparison directory
RUNS_DIR=$(realpath "../../zymobiomics_folder/results/runs")

# Output CSV file path
OUTPUT="$RUNS_DIR/new_runs_summary.csv"
SAMPLES=()

metadata_vars=(run_id db_name runtime)
preqc_vars=(preqc_pct_dups_r1 preqc_pct_dups_r2 preqc_pct_gc_r1 preqc_pct_gc_r2 preqc_avg_len_r1 preqc_avg_len_r2 preqc_med_len_r1 preqc_med_len_r2 preqc_fail_r1 preqc_fail_r2)
trim_vars=(trim_clip trim_head trim_lead trim_crop trim_sliding_win trim_trail trim_min_len trim_paired trim_both trim_fw_only trim_rev_only trim_dropped)
postqc_vars=(postqc_pct_dups_r1 postqc_pct_dups_r2 postqc_pct_gc_r1 postqc_pct_gc_r2 postqc_avg_len_r1 postqc_avg_len_r2 postqc_med_len_r1 postqc_med_len_r2 postqc_fail_r1 postqc_fail_r2)
bt_vars=(bt_prefix bt_mode bt_sensitivity bt_mixed bt_discordant bt_paired bt_conc_0 bt_conc_1 bt_conc_more bt_avg_len bt_med_len)
kraken2_vars=(kraken2_min_hits kraken2_total kraken2_classified kraken2_unclassified kraken2_avg_len kraken2_med_len)
bracken_vars=(bracken_thresh bracken_species bracken_species_above_thresh bracken_species_below_thresh bracken_kept bracken_discarded bracken_redistributed bracken_not_redistributed bracken_total)
eval_metrics=(precision recall)

  # Write CSV header only if file doesn't exist
if [[ ! -f "$OUTPUT" ]]; then
  {
    IFS=,
    echo "${metadata_vars[*]},${preqc_vars[*]},${trim_vars[*]},${postqc_vars[*]},${bt_vars[*]},${kraken2_vars[*]},${bracken_vars[*]},${eval_metrics[*]}"
    unset IFS
  } > "$OUTPUT"
fi

# Find all kraken_pipeline.log files
LOG_FILES=$(find "$RUNS_DIR" -type f -name "kraken_pipeline.log")

# Process each log file
for log in $LOG_FILES; do
    
  # -- RUN METADATA -- 

  # Get run ID from log file's grandparent directory
  run_id=$(basename "$(dirname "$(dirname "$log")")") 

  # Skip if run already exists in the output
  if grep -q "^$run_id," "$OUTPUT"; then
      echo "Skipping $run_id (already in summary)"
      continue
  fi

  # Get database name (basename of the path)
  db_name=$(grep "Kraken2/Bracken Database Path:" "$log" | awk -F': ' '{print $2}' | xargs basename)

  # Get runtime string (e.g., 00:01:45)
  runtime=$(grep "Metagenomic classification completed in:" "$log" | awk -F': ' '{print $2}')    
    
  # -- FASTQC & TRIMMING METADATA --

  if grep -q "Quality Control & Trimming: Enabled" "$log"; then
    echo "Processing FastQC & Trimming for $run_id"  
    # Initialize pre-trimming FastQC variables (placeholders)
    for var in ${preqc_vars[@]}; do
      declare "$var=_"
    done
    echo "Pre-trimming FastQC variables initialized"

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
    
    # Check output
    echo "$trim_clip,$trim_head,$trim_lead,$trim_crop,$trim_sliding_win,$trim_trail,$trim_min_len"

    # Extract per-sample trimming stats and sample IDs
    grep -A2 "Input Read Pairs:" "$log" | while read -r line; do
      # Skip grep's -- separator line
      [[ "$line" == "--" ]] && continue

      # Save the line with trimming stats
      if [[ "$line" =~ ^Input\ Read\ Pairs: ]]; then
        combined_line="$line"
        continue
      fi

      # Extract sample ID and parse trimming values
      if [[ "$line" =~ ^✅\ Trimmed ]]; then
        sample_id=$(echo "$line" | sed -n 's/^✅ Trimmed \(.*\)!.*/\1/p')

        # Extract values from saved stats line
        trim_paired=$(echo "$combined_line" | sed -n 's/.*Input Read Pairs: \([0-9]*\).*/\1/p')
        trim_both=$(echo "$combined_line" | sed -n 's/.*Both Surviving: \([0-9]*\).*/\1/p')
        trim_fw_only=$(echo "$combined_line" | sed -n 's/.*Forward Only Surviving: \([0-9]*\).*/\1/p')
        trim_rev_only=$(echo "$combined_line" | sed -n 's/.*Reverse Only Surviving: \([0-9]*\).*/\1/p')
        trim_dropped=$(echo "$combined_line" | sed -n 's/.*Dropped: \([0-9]*\).*/\1/p')

        # Store stats for sample as a single variable
        declare "trimming_stats_$sample_id=$trim_paired,$trim_both,$trim_fw_only,$trim_rev_only,$trim_dropped"
        varname="trimming_stats_$sample_id"
        echo "${!varname}"
      fi
    done
    
    # Initialize post-trimming FastQC variables (placeholders)
    for var in ${postqc_vars[@]}; do
      declare "$var=_"
    done
    echo "Post-trimming FastQC variables initialized"
  fi

  # -- HOST-DNA REMOVAL METADATA --
  if grep -q "Host DNA Removal: Enabled" "$log"; then
    echo "Host DNA Removal for $run_id"  

    # Extract Bowtie2 arguments (from first occurrence)
    # CREATE PLACEHOLDER UNTIL THE COMMAND IS ECHOED IN THE ACTUAL SCRIPT 
    for var in $bt_prefix $bt_mode $bt_sensitivity $bt_mixed $bt_discordant $bt_paired; do
      declare "$var=_"
    done  
    echo "Global Bowtie2 variables initialized"
    
    # Parse Bowtie2 stats line-by-line for each sample
    grep -A7 "Processing sample:" "$log" | while read -r line; do
      [[ $line =~ ^Processing\ sample:\  ]] && sample_id=${line#*: } && continue
      [[ $line =~ reads\; ]] && bt_paired=$(echo "$line" | sed -E 's/^([0-9]+) .*/\1/') && continue
      [[ $line =~ concordantly\ 0\ times ]] && bt_conc_0=$(echo "$line" | sed -E 's/^[[:space:]]*([0-9]+).*/\1/') && continue
      [[ $line =~ exactly\ 1\ time ]] && bt_conc_1=$(echo "$line" | sed -E 's/^[[:space:]]*([0-9]+).*/\1/') && continue
      if [[ $line =~ \>1\ times ]]; then
        bt_conc_more=$(echo "$line" | sed -E 's/^[[:space:]]*([0-9]+).*/\1/')
        
        # Store stats in sample-specific variable once full block is complete       
        declare "bowtie_stats_$sample_id=$bt_paired,$bt_conc_0,$bt_conc_1,$bt_conc_more"
        varname="bowtie_stats_$sample_id"
        echo "${!varname}"
      fi        
    done

  fi 

# -- KRAKEN & BRACKEN METADATA -- 
kraken2_min_hits=$(grep "kraken2" "$log" | grep -oE -- '--minimum-hit-groups[ =][0-9]+' | grep -oE '[0-9]+' | head -n1)
bracken_thresh=$(grep "Threshold:" "$log" | sed -n 's/.*Threshold: \([0-9]*\).*/\1/p' | head -n1)

# Parse Kraken2 & Bracken summary blocks from log (split by "Processing sample: ")
# Collect samples from awk
while IFS=',' read -r sample_id minhit total classified unclassified avg med thresh species above below total_brack; do
  declare "kraken_bracken_stats_$sample_id=$minhit,$total,$classified,$unclassified,$avg,$med,$thresh,$species,$above,$below,$total_brack"
  varname="kraken_bracken_stats_$sample_id"
  echo "${!varname}"
done < <(
  awk -v bracken_thresh="$bracken_thresh" -v minhits="$kraken2_min_hits" -v RS="Processing sample: " '    
  BEGIN { FS = "\n" }
  NR > 1 {
    sample = $1 # Sample ID is first line of each block
    print sample  

    # Initialize output fields with placeholders
    kraken2_total = kraken2_classified = kraken2_unclassified = kraken2_avg_len = kraken2_med_len = "-"
    bracken_species = bracken_species_above_thresh = bracken_species_below_thresh = bracken_kept = "-"
    bracken_discarded = bracken_redistributed = bracken_not_redistributed = bracken_total = "-"
    
    # Loop over each line in the block to extract relevant stats
    # Extract Kraken2 classification summary values
    for (i = 1; i <= NF; i++) {
      # Kraken2 values
      if ($i ~ /sequences \(/) {
        kraken2_total = $i
        gsub(/[^0-9].*$/, "", kraken2_total) # Keep digits before first non-digit
      } 
      if ($i ~ /sequences classified \(/) {
        kraken2_classified = $i
        gsub(/\([^)]*\)/, "", kraken2_classified) # Remove percentage in parentheses
        gsub(/[^0-9]/, "", kraken2_classified) # Keep only digits
      }
      if ($i ~ /sequences unclassified \(/) {
        kraken2_unclassified = $i
        gsub(/\([^)]*\)/, "", kraken2_unclassified)
        gsub(/[^0-9]/, "", kraken2_unclassified)
      }

      # Extract Bracken species abundance summary values
      if ($i ~ /Number of species in sample:/)      { bracken_species = $i; gsub(/[^0-9]/, "", bracken_species) }
      if ($i ~ /reads > threshold:/)                { bracken_species_above_thresh = $i; gsub(/[^0-9]/, "", bracken_species_above_thresh) }
      if ($i ~ /reads < threshold:/)                { bracken_species_below_thresh = $i; gsub(/[^0-9]/, "", bracken_species_below_thresh) }
      if ($i ~ /reads kept at species level/)       { bracken_kept = $i; gsub(/[^0-9]/, "", bracken_kept) }
      if ($i ~ /reads discarded/)                   { bracken_discarded = $i; gsub(/[^0-9]/, "", bracken_discarded) }
      if ($i ~ /Reads distributed:/)                { bracken_redistributed = $i; gsub(/[^0-9]/, "", bracken_redistributed) }
      if ($i ~ /Reads not distributed/)             { bracken_not_redistributed = $i; gsub(/[^0-9]/, "", bracken_not_redistributed) }
      
      # Compute total reads considered at species level
      bracken_total = bracken_kept + bracken_redistributed
    }

    # Output all values in a single CSV-formatted line
    if (sample != "-" && kraken2_total != "-") {
    printf "%s,%s,%s,%s,%s,%s,%s,", sample, minhits, kraken2_total, kraken2_classified, kraken2_unclassified, kraken2_avg_len, kraken2_med_len
    printf "%s,%s,%s,%s,%s,%s,%s,%s,%s\n", bracken_thresh, bracken_species, bracken_species_above_thresh, bracken_species_below_thresh, bracken_kept, bracken_discarded, bracken_redistributed, bracken_not_redistributed, bracken_total
}
  }' "$log"
)
echo "${kraken_bracken_stats_ZC1_S4}"
done


