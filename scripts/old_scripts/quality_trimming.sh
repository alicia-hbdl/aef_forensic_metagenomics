#!/bin/bash
#SBATCH --output=logs/trimming.log  
#SBATCH --error=logs/trimming.err  
#SBATCH --job-name=trimmomatic
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=msc_appbio

# Check if exactly one argument is provided
if [ "$#" -ne 1 ]; then
    echo -e "\n❗ Usage: $0 <path/to/raw_reads/>"
    exit 1
fi

# Validate that the argument is a directory
raw_reads="$1"
if [ ! -d "$raw_reads" ]; then
    echo -e "\n❌ Error: Directory does not exist!"
    exit 1
fi

# Check for uncompressed FASTQ files
uncompressed_files=$(find "$raw_reads" -type f \( -name '*.fastq' -o -name '*.fq' \))
if [ -n "$uncompressed_files" ]; then
    echo "\nCompressing uncompressed FASTQ files..."
    find "$raw_reads" -type f \( -name '*.fastq' -o -name '*.fq' \) -exec gzip {} +
    echo "✅ Compression complete."
fi

# Check for `.fastq.gz` files
file_ext=$(find "$raw_reads" -maxdepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' \) | head -n 1 | grep -Eo '\.(fastq.gz|fq.gz)' || echo "")

if [ -n "$file_ext" ]; then
    echo -e "\n✅ ${file_ext#.} files detected in $raw_reads."

else
    echo -e "\n❌ No FASTQ files found in $raw_reads!"
    exit 1
fi

# Define output directories
fastqc_pre="$raw_reads/fastqc_pre_trim"
fastqc_post="$raw_reads/fastqc_post_trim"
trimmed_reads="$raw_reads/trimmed"

# Create and clean directories efficiently
for dir in "$fastqc_pre" "$fastqc_post" "$trimmed_reads"; do
    rm -rf "$dir" && mkdir -p "$dir"
done

# Run FastQC on raw reads before trimming
# 1 represents stdout (standard output), 2 represents stderr (standard error), and /dev/null discards the redirected output.
echo -e "\nRunning FastQC on all raw reads..."
if fastqc "$raw_reads"/*"$file_ext" --outdir "$fastqc_pre" 2>&1; then
    echo "✅ Pre-trimming FastQC completed successfully."
else
    echo "❌ FastQC failed!" 
    exit 1
fi

if multiqc "$fastqc_pre" --no-data-dir -o "$fastqc_pre" 2>&1; then
    echo "✅ Pre-trimming MultiQC report generated successfully." && rm -f "$fastqc_pre"/*.zip  
else
    echo "❌ MultiQC failed!" 
    exit 1
fi

# Define adapter file for trimming
ADAPT_FA="/scratch/users/k24087895/final_project/metagenome_analysis/data/adapters/zymobiomics_adaptors.fa"

# Define Trimmomatic parameters and steps 
TRIM_PARAM="PE -phred33 -threads 4"
STEPS="ILLUMINACLIP:$ADAPT_FA:2:30:10:8:true HEADCROP:5 LEADING:5 CROP:240 SLIDINGWINDOW:4:15 TRAILING:10 MINLEN:30"
# Look at maxinfo 

# Loop through all forward reads
echo -e "\nRunning Trimmomatic..."
for R1 in $raw_reads/*_L001_R1_001$file_ext; do
    [[ -f "$R1" ]] || continue  # Skip if no matching file is found
    
    base=$(basename "$R1" "_L001_R1_001$file_ext")  
    R2="$raw_reads/${base}_L001_R2_001$file_ext"

    # Check if the matching reverse read exists
    if [ ! -f "$R2" ]; then
        echo "❗ Warning: Reverse read missing for $R1. Skipping..."
        continue
    fi

    # Define output file paths for Trimmomatic
    R1_paired="$trimmed_reads/${base}_L001_R1_001_paired$file_ext"
    R1_unpaired="$trimmed_reads/${base}_L001_R1_001_unpaired$file_ext"
    R2_paired="$trimmed_reads/${base}_L001_R2_001_paired$file_ext"
    R2_unpaired="$trimmed_reads/${base}_L001_R2_001_unpaired$file_ext"

    # Run Trimmomatic **on compressed FASTQ files**
    trimmomatic $TRIM_PARAM "$R1" "$R2" \
        "$R1_paired" "$R1_unpaired" \
        "$R2_paired" "$R2_unpaired" \
        $STEPS 2>&1 && echo -e  "✅ Trimmed $base!" || echo -e "❌ Trimmomatic failed for $base!"
done 

# Run FastQC on trimmed reads to check resulting quality
echo -e "\nRunning FastQC on paired trimmed readss..."
if fastqc "$raw_reads"/trimmed/*_paired"$file_ext" --outdir "$fastqc_post" 2>&1; then
    echo "✅ Post-trimming FastQC completed successfully."
else
    echo "❌ FastQC failed!" 
    exit 1
fi

if multiqc "$fastqc_post" --no-data-dir -o "$fastqc_post" 2>&1; then
    echo "✅ Post-trimming MultiQC report generated successfully." && rm -f "$fastqc_post"/*.zip  
else
    echo "❌ MultiQC failed!" 
    exit 1
fi

