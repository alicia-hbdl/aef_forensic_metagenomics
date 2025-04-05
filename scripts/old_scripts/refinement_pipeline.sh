#!/bin/bash

# Ensure a basename is provided as an argument  
if [ -z "$1" ]; then  
    echo -e "\n❗Usage: $0 <basename> [-r | --remove-human-dna]"  
    exit 1  
fi  
BASENAME="$1"  

REMOVE_HUMAN_DNA=false
# Check if the remove human DNA flag is provided
if [ "$2" == "-r" ] || [ "$2" == "--remove-human-dna" ]; then
    REMOVE_HUMAN_DNA=true
fi

# Define root directory & project folder  
ROOT_DIR="/scratch/users/k24087895/final_project/metagenome_analysis"
LOG_FILE="$ROOT_DIR/homemade_scripts/logs/forensics_pipeline.log"  # Path to the pipeline log file
FOLDER="$ROOT_DIR/${BASENAME}_folder"

# Check & Download Kraken2 database if missing  
KRAKEN_DATABASE="$ROOT_DIR/data/databases/k2protocol_db/"
if [ ! -d "$KRAKEN_DATABASE" ] || [ -z "$(ls -A "$KRAKEN_DATABASE" 2>/dev/null)" ]; then  
    echo -e "\n❗ Kraken2 database not found. Downloading..."  
    mkdir -p "$KRAKEN_DATABASE"  
    wget -O "$KRAKEN_DATABASE/k2_db.tar.gz" "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_eupath_20201202.tar.gz" && \
    tar -xvzf "$KRAKEN_DATABASE/k2_db.tar.gz" -C "$KRAKEN_DATABASE" && \
    rm "$KRAKEN_DATABASE/k2_db.tar.gz" && \
    echo -e "✅ Kraken2 database ready!" || { echo -e "❌ Error: Kraken2 download failed!"; exit 1; }
else  
    echo -e "\n✅ Kraken2 database found."
fi  

# Check & Download GRCh38 index if missing if remove human DNA is enabled
GRCh38_INDEX="$ROOT_DIR/data/bowtie_index/GRCh38_noalt_as/GRCh38_noalt_as"
if [ "$REMOVE_HUMAN_DNA" = true ] && [ ! -f "${GRCh38_INDEX}.1.bt2" ]; then  
    echo "\n❗GRCh38 index not found. Downloading..."  
    mkdir -p "$ROOT_DIR/data/bowtie_index"  
    wget -O "$ROOT_DIR/data/bowtie_index/GRCh38_noalt_as.zip" "https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip" && \
    unzip -d "$ROOT_DIR/data/bowtie_index/" "$ROOT_DIR/data/bowtie_index/GRCh38_noalt_as.zip" && \
    rm "$ROOT_DIR/data/bowtie_index/GRCh38_noalt_as.zip" && \
    echo -e "✅ GRCh38 reference genome index ready!" || { echo -e "❌ Error: GRCh38 download failed!"; exit 1; }
else  
    echo -e "\n✅ GRCh38 reference genome index found."
fi  

# Define directories  
DATA="$FOLDER/data"
HUMAN_READS="$DATA/human_reads"
METAGENOMIC_READS="$DATA/metagenomic_reads"
TRIMMED_FASTQ="$DATA/raw_fastq/trimmed"

# Ensure TRIMMED_FASTQ directory exists
if [ ! -d "$TRIMMED_FASTQ" ]; then
  echo -e "\n❌ Error: Directory '$TRIMMED_FASTQ' does not exist!"
  exit 1
fi

# Detect and validate FASTQ file extension
file_ext=$(find "$TRIMMED_FASTQ" -maxdepth 1 -type f | grep -Eo '\.(fastq\.gz|fq\.gz)$' | head -n 1)
if [[ "$file_ext" =~ \.(fastq\.gz|fq\.gz) ]]; then
  echo -e "\n✅ Valid $file_ext files detected in $TRIMMED_FASTQ."
else
  echo -e "\n❌ Error: No valid files found in '$TRIMMED_FASTQ'!"
  exit 1
fi

# Define subdirectories  
ALIGNED_SAM="$HUMAN_READS/aligned_sam"
FILTERED_FASTQ="$METAGENOMIC_READS/filtered_fastq"
KRAKEN2="$METAGENOMIC_READS/kraken2"
K2REPORT="$METAGENOMIC_READS/k2report"
BRACKEN="$METAGENOMIC_READS/bracken"
BREPORT="$METAGENOMIC_READS/breport"
KRONA_TXT="$METAGENOMIC_READS/krona_txt"
TOTAL_READS="$METAGENOMIC_READS/total_reads"
CLASSIFIED="$METAGENOMIC_READS/classified"
UNCLASSIFIED="$METAGENOMIC_READS/unclassified"
# Define output directories  
OUTPUTS="$FOLDER/outputs"
DIVERSITY="$OUTPUTS/diversity"
KRONA="$OUTPUTS/krona"
PAVIAN="$OUTPUTS/pavian"

# Create required directories  
mkdir -p "$ALIGNED_SAM" "$FILTERED_FASTQ" "$K2REPORT" "$BRACKEN" \
         "$BREPORT" "$KRONA_TXT" "$TOTAL_READS" "$DIVERSITY" "$KRONA" \
         "$PAVIAN" "$UNCLASSIFIED" "$CLASSIFIED" "$KRAKEN2"
echo "✅ Directory structure created under: $FOLDER"

# Define paths to external tools  
TOOLS="$ROOT_DIR/tools"
KREPORT2KRONA="$TOOLS/KrakenTools/kreport2krona.py"
IMPORT_TEXT="$TOOLS/Krona/KronaTools/scripts/ImportText.pl"
A_DIVERSITY="$TOOLS/KrakenTools/DiversityTools/alpha_diversity.py"
B_DIVERSITY="$TOOLS/KrakenTools/DiversityTools/beta_diversity.py"
B_DIVERSITY_HEATMAP="$ROOT_DIR/homemade_scripts/b_diversity_heatmap.R"
HUMAN_DNA_ANALYSIS="$ROOT_DIR/homemade_scripts/human_alignment.sh"

cd "$TRIMMED_FASTQ" || { echo "❌ Error: Cannot access '$TRIMMED_FASTQ'"; exit 1; } # Navigate to TRIMMED_FASTQ directory & exit if it fails
#mapfile -t SAMPLES < <(find . -maxdepth 1 -type f -name "*$file_ext" | sed -E "s/_L001_R[12]_001.*$//" | sort -u) # Extract unique base sample names (ignoring paired/unpaired)
mapfile -t SAMPLES < <(ls *$file_ext | sed -E "s/_L001_R[12]_001.*$//" | sort -u)

echo -e "\nSamples: ${SAMPLES[@]}" # Display the extracted sample names

# Loop through each sample and process with Kraken2, Bracken, and diversity tools
for SAMPLE in "${SAMPLES[@]}"; do
    
    echo -e "\nProcessing sample $SAMPLE:"

    # Remove human reads with Bowtie2 (GRCh38, 8 threads)  
    # --end-to-end ensures full-read alignment | --very-sensitive enhances accuracy but slows alignment  
    # --no-mixed prevents single-end mismapping | --no-discordant retains only properly paired reads
    # --un-conc saves paired-end reads that fail to align concordantly | -S stores all reads (aligned + unaligned) in SAM format
    if [ "$REMOVE_HUMAN_DNA" = true ]; then
	      echo -e "Removing host human DNA with Bowtie2..."
        bowtie2 -x "$GRCh38_INDEX" -p 8 -q \
            --end-to-end --very-sensitive --no-mixed --no-discordant \
            -1 "$TRIMMED_FASTQ/${SAMPLE}_L001_R1_001_paired$file_ext" -2 "$TRIMMED_FASTQ/${SAMPLE}_L001_R2_001_paired$file_ext" \
            --un-conc "$FILTERED_FASTQ/${SAMPLE}_nonhuman" \
            -S "$ALIGNED_SAM/${SAMPLE}_alignments.sam" 2>&1
        echo " ✅ Human reads removed with Bowtie2."

        # Gzip and rename output files
        gzip -c "$FILTERED_FASTQ/${SAMPLE}_nonhuman.1" > "$FILTERED_FASTQ/${SAMPLE}_nonhuman.1$file_ext"
        gzip -c "$FILTERED_FASTQ/${SAMPLE}_nonhuman.2" > "$FILTERED_FASTQ/${SAMPLE}_nonhuman.2$file_ext"
        rm "$FILTERED_FASTQ/${SAMPLE}_nonhuman.1" "$FILTERED_FASTQ/${SAMPLE}_nonhuman.2"
    
        echo " ✅ Compressed and renamed non-human read files."
    else
        cp "$TRIMMED_FASTQ/${SAMPLE}_L001_R1_001_paired$file_ext" "$FILTERED_FASTQ/${SAMPLE}_nonhuman.1$file_ext"
        cp "$TRIMMED_FASTQ/${SAMPLE}_L001_R2_001_paired$file_ext" "$FILTERED_FASTQ/${SAMPLE}_nonhuman.2$file_ext"
        echo "❗Skipped human read removal with Bowtie2."
    fi
    
    # Classify non-human reads with Kraken2 (8 threads)  
    # --minimum-hit-groups 3: Require ≥3 k-mers | --report: Save report | --report-minimizer-data: Add minimizer stats  
    kraken2 --db "$KRAKEN_DATABASE" --threads 8 --report "$K2REPORT/${SAMPLE}.k2report" \
            --report-minimizer-data --paired --minimum-hit-groups 2 --gzip-compressed \
            --classified-out "$CLASSIFIED/${SAMPLE}#.fq" --unclassified-out "$UNCLASSIFIED/${SAMPLE}#.fq" \
            "$FILTERED_FASTQ/${SAMPLE}_nonhuman.1$file_ext" "$FILTERED_FASTQ/${SAMPLE}_nonhuman.2$file_ext" \
            --output "$KRAKEN2/${SAMPLE}.kraken2" --use-names 2>&1 
    rm -f "$FILTERED_FASTQ/${SAMPLE}_nonhuman.1$file_ext" "$FILTERED_FASTQ/${SAMPLE}_nonhuman.2$file_ext"
    echo " ✅ Non-human reads classified with Kraken2."

    # Estimate microbial abundance with Bracken  
    # -r 100: Read length | -l S: Species-level abundance | -t 10: Min reads for estimation | -w: Human-readable report  
    bracken -d "$KRAKEN_DATABASE" -i "$K2REPORT/${SAMPLE}.k2report" -l S \
            -r 100 -t 10 -o "$BRACKEN/${SAMPLE}.bracken" -w "$BREPORT/${SAMPLE}.breport" 2>&1
    echo " ✅ Microbial abundance estimated with Bracken."

    # Generate Krona visualization  
    # --no-intermediate-ranks: Major taxonomic levels only  
    python "$KREPORT2KRONA" -r "$BREPORT/${SAMPLE}.breport" -o "$KRONA_TXT/${SAMPLE}.krona.txt" --no-intermediate-ranks && \
    "$IMPORT_TEXT" "$KRONA_TXT/${SAMPLE}.krona.txt" -o "$KRONA/${SAMPLE}.krona.html" && \
    rm "$KRONA_TXT/${SAMPLE}.krona.txt" 
    echo " ✅ Krona plot generated."

    # Compute Alpha Diversity for each metric  
    for METRIC in "BP" "Sh" "F" "Si" "ISi"; do  
        python "$A_DIVERSITY" -f "$BRACKEN/${SAMPLE}.bracken" -a "$METRIC"  
    done  
    echo " ✅ Alpha diversity analysis completed."

done

# Extract the total number of reads processed for each sample for normalization 
awk '
  BEGIN {print "sample,total_reads"}  # Print CSV header
  /Processing sample/ {sample = $3; sub(/:$/, "", sample)}  # Extract sample name and clean trailing colon
  /Total reads in sample:/ {print sample "," $6}  # Extract total reads and print
  ' "$LOG_FILE" > "$TOTAL_READS/total_reads.csv"  # Read from log file and write to CSV
echo "✅ Total reads extracted in: $TOTAL_READS/total_reads.csv"

# Run Beta Diversity and Alignment Summary if multiple samples exist  
if [ "${#SAMPLES[@]}" -gt 1 ]; then  
  # Extracting the path to all bracken utput files
  INPUT_FILES=""  # Initialize variable
  for SAMPLE in "${SAMPLES[@]}"; do
    INPUT_FILES+="$BRACKEN/${SAMPLE}.bracken "  
  done
  
  # Compute b-diversity with KrakenTools
  python "$B_DIVERSITY" -i ${INPUT_FILES} --type bracken > "$DIVERSITY/beta_diversity_matrix.tsv" && \
  rm -f ${INPUT_FILES}  
  
  # Generate heatmap to illustrate b-divesrity
  Rscript "$B_DIVERSITY_HEATMAP" "$DIVERSITY/beta_diversity_matrix.tsv" 
  echo "✅ Beta Diversity heatmap generated."

  # Analyze human read alignment 
  if [ "$REMOVE_HUMAN_DNA" = true ]; then
      echo "Running alignment summary on human DNA..."
      "$HUMAN_DNA_ANALYSIS" "$FOLDER"
  else  
      echo -e "❗ No human DNA removed.\n"
  fi  

else  
    echo "❗Skipping Beta Diversity and Alignment Summary: Only one sample detected."
fi

echo "Pipeline completed."
