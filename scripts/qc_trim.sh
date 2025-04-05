#!/bin/bash

# Ensure the FASTQ directory is provided and valid
if [[ -z "$1" || "$1" == -* || ! -d "$1" || -z "$(find "$1" -maxdepth 1 -type f \( -name "*.fq" -o -name "*.fastq" -o -name "*.fq.gz" -o -name "*.fastq.gz" \) 2>/dev/null)" ]]; then
    echo "❌ Error: Invalid or missing FASTQ directory. Ensure it contains .fq, .fastq, .fq.gz, or .fastq.gz files."
    exit 1
fi

# Set the RAW_FASTQ_DIR to the directory passed as argument
RAW_FASTQ_DIR="$1"

# Define and create output directories for FastQC and trimmed data
# These are created in the main script but are also initialized here for standalone use
FASTQC_DIR="$RAW_FASTQ_DIR/../results/fastqc"
mkdir -p "$FASTQC_DIR/pre_trimming" "$FASTQC_DIR/post_trimming"
TRIMMED_DIR="$RAW_FASTQ_DIR/../processed_data/trimmed"
mkdir -p "$TRIMMED_DIR/paired" "$TRIMMED_DIR/unpaired"

# Run FastQC on all raw reads
echo -e "\nRunning FastQC on all raw reads..."
fastqc "$RAW_FASTQ_DIR"/*.fastq.gz --outdir "$FASTQC_DIR/pre_trimming" &>/dev/null || {
    echo "❌ FastQC failed!"
    exit 1
}
echo "✅ FastQC completed successfully."

# Run MultiQC to summarize FastQC reports and clean up zip files
multiqc "$FASTQC_DIR/pre_trimming" --no-data-dir -o "$FASTQC_DIR/pre_trimming" --force 2>&1 || {
    echo "❌ MultiQC failed!"
    exit 1
}
echo -e "✅ MultiQC report generated successfully.\n"
rm -f "$FASTQC_DIR/pre_trimming"/*.zip  

# Define Trimmomatic adapter file, trimming parameters, and steps
ADAPT_FA="$ROOT_DIR/data/adapters/zymobiomics_adaptors.fa"
TRIM_PARAM="PE -phred33 -threads 4"
STEPS="ILLUMINACLIP:$ADAPT_FA:2:30:10:8:true HEADCROP:5 LEADING:5 CROP:240 SLIDINGWINDOW:4:15 TRAILING:10 MINLEN:35"

echo -e "Running Trimmomatic..."
for R1 in "$RAW_FASTQ_DIR"/*_R1*.fastq.gz; do  # Handle variations of sample names
    [[ -f "$R1" ]] || {
        echo "❌ Error: No matching files found in $RAW_FASTQ_DIR"
        exit 1
    }
    
    # Extract the base sample name by removing lane, read, and other suffixes
    base=$(basename "$R1" | sed -E 's/_L[0-9]+_R[12]_?[0-9]*\.fastq\.gz$//; s/_R[12]_?[0-9]*\.fastq\.gz$//')
    R2=$(echo "$R1" | sed 's/_R1/_R2/')  # Identify the corresponding reverse read

    # Check that the reverse read file exists
    if [[ ! -f "$R2" ]]; then
        echo "⚠️ Warning: Reverse read missing for $R1. Skipping..."
        continue
    fi

    # Define output file paths for paired and unpaired reads
    R1_paired="$TRIMMED_DIR/paired/${base}_R1_paired.fastq.gz"
    R1_unpaired="$TRIMMED_DIR/unpaired/${base}_R1_unpaired.fastq.gz"
    R2_paired="$TRIMMED_DIR/paired/${base}_R2_paired.fastq.gz"
    R2_unpaired="$TRIMMED_DIR/unpaired/${base}_R2_unpaired.fastq.gz"

    # Run Trimmomatic and log success or failure
    if trimmomatic $TRIM_PARAM "$R1" "$R2" "$R1_paired" "$R1_unpaired" "$R2_paired" "$R2_unpaired" $STEPS 2>&1; then
        echo -e "✅ Trimmed $base!\n"
    else
        echo -e "❌ Trimmomatic failed for $base!"
        exit 1
    fi
done

# Run FastQC on trimmed paired reads
echo -e "\nRunning FastQC on paired trimmed reads..."
fastqc "$TRIMMED_DIR/paired"/*.fastq.gz --outdir "$FASTQC_DIR/post_trimming" &>/dev/null || {
    echo "❌ FastQC failed!"
    exit 1
}
echo "✅ FastQC completed successfully."

# Run MultiQC to summarize post-trimming FastQC results and clean up
multiqc "$FASTQC_DIR/post_trimming" --no-data-dir -o "$FASTQC_DIR/post_trimming" --force 2>&1 || {
    echo "❌ MultiQC failed!"
    exit 1
}
echo -e "✅ MultiQC report generated successfully.\n"
rm -f "$FASTQC_DIR/post_trimming"/*.zip  

echo "✅ All processing steps completed successfully."