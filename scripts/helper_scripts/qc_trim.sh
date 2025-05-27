#!/bin/bash

# This is a pipeline that performs quality control and trimming of FASTQ files using FastQC and Trimmomatic.
# It generates reports before and after trimming, and handles paired-end reads.

# Usage: ./qc_trim.sh -r|--reads <FASTQ_DIR> -a|--adapters <adapter_file.fa>

# TODO: Eventually allow the user to specify TRIM_PARAM and STEPS as arguments

# Basic parsing only (adapter file and reads directory are validated in the main pipeline)
while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--reads) RAW_FASTQ_DIR="$2"; shift 2 ;; # Set the RAW_FASTQ_DIR to the directory passed as argument
    -a|--adapters) ADAPTER_FILE="$2"; shift 2 ;;
    *) echo "❌ Unknown parameter: $1"; exit 1 ;;
  esac
done

# Define and create output directories for FastQC and trimmed data
# These are also initialized here to support standalone use of this script
FASTQC_DIR="$RAW_FASTQ_DIR/../results/fastqc"
mkdir -p "$FASTQC_DIR/pre_trimming" "$FASTQC_DIR/post_trimming"
TRIMMED_DIR="$RAW_FASTQ_DIR/../processed_data/trimmed"
mkdir -p "$TRIMMED_DIR/paired" "$TRIMMED_DIR/unpaired"

# ------------ PRE-TRIMMING QC ------------ 

echo -e "\nRunning FastQC on all raw reads..."
fastqc "$RAW_FASTQ_DIR"/*.fastq.gz --outdir "$FASTQC_DIR/pre_trimming" >&2 || { echo "❌ FastQC failed!"; exit 1; }
echo "✅ FastQC completed successfully."

# Run MultiQC to summarize FastQC reports and clean up zip files
multiqc "$FASTQC_DIR/pre_trimming" --no-data-dir -o "$FASTQC_DIR/pre_trimming" --force 2>&1 || { echo "❌ MultiQC failed!"; exit 1; }
echo -e "✅ MultiQC report generated successfully.\n"
rm -f "$FASTQC_DIR/pre_trimming"/*.zip  

# ------------ TRIMMING ------------ 

# Define trimming parameters and steps
TRIM_PARAM="PE -phred33 -threads 4"
STEPS="ILLUMINACLIP:$ADAPTER_FILE:2:30:10:8:true HEADCROP:5 LEADING:5 CROP:240 SLIDINGWINDOW:4:15 TRAILING:10 MINLEN:35"

echo -e "Running Trimmomatic..."
for R1 in "$RAW_FASTQ_DIR"/*_R1*.fastq.gz; do  # Handle variations of sample names
    
    # Checks if the forward read file exists
    [[ -f "$R1" ]] || { echo "❌ Error: No matching files found in $RAW_FASTQ_DIR"; exit 1; }

    # Extract the base sample name by removing lane, read, and other suffixes
    base=$(basename "$R1" | sed -E 's/_L[0-9]+_R[12]_?[0-9]*\.fastq\.gz$//; s/_R[12]_?[0-9]*\.fastq\.gz$//')
    R2=$(echo "$R1" | sed 's/_R1/_R2/')  # Identify the corresponding reverse read

    # Check that the reverse read file exists
    [[ -f "$R2" ]] || { echo "⚠️ Warning: Reverse read missing for $R1. Skipping..."; continue; }
  
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

# ------------ POST-TRIMMING QC ------------

echo -e "\nRunning FastQC on paired trimmed reads..."
fastqc "$TRIMMED_DIR/paired"/*.fastq.gz --outdir "$FASTQC_DIR/post_trimming" >&2 || { echo "❌ FastQC failed!"; exit 1; }
echo "✅ FastQC completed successfully."

multiqc "$FASTQC_DIR/post_trimming" --no-data-dir -o "$FASTQC_DIR/post_trimming" --force 2>&1 || { echo "❌ MultiQC failed!"; exit 1; }
echo -e "✅ MultiQC report generated successfully.\n"
rm -f "$FASTQC_DIR/post_trimming"/*.zip  
