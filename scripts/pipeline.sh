#!/bin/bash

#SBATCH --job-name=pipeline
#SBATCH --output=logs/pipeline.log 
#SBATCH --error=logs/pipeline.err
#SBATCH --time=02:00:00
#SBATCH --mem=44G # Check sacct to determine appropriate values 
#SBATCH --cpus-per-task=8  # Match requested threads 
#SBATCH --ntasks=1  
#SBATCH --partition=cpu 

set -e # Exit on error 

<< 'COMMENT'
=============================================================================================
Usage:
    ./pipeline.sh --raw-fastq/-fq <reads_dir> [--database/-db <database_path>] [-t | --trim] [-r | --remove-host-dna <index_path>]

Arguments:
    -fq, --raw-fastq <path>       Path to the directory containing raw FASTQ reads to process.
    -db, --database <path>        (Optional) Path to the Kraken2/Bracken database used for metagenomic classification.
    -t, --trim                    (Optional) Enable adapter and quality trimming of reads.
    -r, --remove-host-dna <path>  (Optional) Remove host DNA contamination.

Notes:
    - Ensure that all necessary dependencies (Kraken2, Trimmomatic, etc.) are installed before running.
=============================================================================================
COMMENT

echo -e "\n================================================= CONDA ENVIRONMENT ACTIVATION =================================================="

# Check for Conda installation and initialize shell integration
if ! command -v conda &> /dev/null; then
    echo "❌ Conda not found. Please install Conda."
    exit 1
fi

eval "$(conda shell.bash hook)" || { echo "❌ Conda shell integration not initialized. Run 'conda init bash'."; exit 1; }

# Create environment if it doesn't exist
if ! conda env list | grep -q "final_project"; then
    echo "❌ 'final_project' environment not found. Creating from environment.yml..."
    conda env create -f environment.yml || { echo "❌ Failed to create 'final_project' environment."; exit 1; }
    echo "✅ 'final_project' environment created."
fi 

# Activate the environment
conda activate final_project || { echo "❌ Failed to activate 'final_project' environment."; exit 1; }
echo "✅ 'final_project' environment activated successfully."

# Use Conda's Python and libraries instead of the system's Python
export PATH="/users/k24087895/.conda/envs/final_project/bin:$PATH"
export PYTHONPATH="/users/k24087895/.conda/envs/final_project/lib/python3.13/site-packages:$PYTHONPATH"

echo -e "\n================================================= ARGUMENT PARSING & VALIDATION ================================================="

# Default values
TRIM=false
REMOVE_HOST_DNA=false

# Parse command-line arguments
while [[ $# -gt 0 ]]; do 
    case "$1" in
    # Process raw FASTQ directory
        -fq|--raw-fastq) 
            # Validate FASTQ directory and its contents
            if [[ -z "$2" || "$2" == -* || ! -d "$2" || -z "$(find "$2" -maxdepth 1 -type f \( -name "*.fq" -o -name "*.fastq" -o -name "*.fq.gz" -o -name "*.fastq.gz" \) 2>/dev/null)" ]]; then
              echo "❌ Error: Invalid or missing FASTQ directory. Ensure it contains .fq, .fastq, .fq.gz, or .fastq.gz files."
              exit 1
            fi
            
            # Compress uncompressed FASTQ files
            for file in "$2"/*.{fq,fastq}; do [[ -f "$file" ]] && gzip "$file" && echo "Compressing '$file'..."; done
            
            # Standardize file extension to .fastq.gz
            for file in "$2"/*.fq.gz; do [[ -f "$file" ]] && mv "$file" "${file%.*}.fastq.gz" && echo "Renaming '$file'..."; done
                    
            # Set raw FASTQ directory path
            RAW_FASTQ_DIR="$2" && echo "Raw FASTQ directory path: $RAW_FASTQ_DIR" 
            
            # Define project root directory and default database/index paths
            PROJ_DIR=$(dirname "$RAW_FASTQ_DIR")  # Current subproject directory
            ROOT_DIR=$(dirname "$PROJ_DIR")  # Root directory with all subprojects and tools
            DATABASE="$ROOT_DIR/data/databases/k2_standard_16gb_20241228"
            BOWTIE_PREFIX="$ROOT_DIR/data/bowtie_index/GRCh38_noalt_as/GRCh38_noalt_as"           
            
            shift 2
            ;;
    
        # Process Kraken2/Bracken database directory argument
        -db|--database) 
            if [[ -z "$2" || "$2" == -* || ! -d "$2" || ! -f "$2/hash.k2d" || ! -f "$2/taxo.k2d" || ! -f "$2/opts.k2d" || ! -f "$2/*kmer_distrib" ]]; then
                echo "⚠️ Warning: Invalid or missing database. Using default: $DATABASE"
            else
                DATABASE="$2"
                echo "✅ Database set to: $DATABASE"
            fi
            shift 2
            ;;
            
      # Enable trimming
      -t|--trim)
          TRIM=true
          echo "✅ Trimming enabled."
          shift
          ;;            
        
      # Process host DNA Bowtie2 index prefix
      -r|--remove-host-dna)  
          if [[ -z "$2" || "$2" == -* || ! -f "$2.1.bt2" || ! -f "$2.2.bt2" || ! -f "$2.3.bt2" || ! -f "$2.4.bt2" || ! -f "$2.rev.1.bt2" || ! -f "$2.rev.2.bt2" ]]; then
              echo "⚠️ Warning: Incomplete or missing Bowtie2 index. Using default: $BOWTIE_PREFIX"
          else
              BOWTIE_PREFIX="$2"
              echo "✅ Bowtie2 index set to: $BOWTIE_PREFIX"
          fi
          REMOVE_HOST_DNA=true
          shift 2
          ;;        
        
      # Handle unknown arguments
      *) 
          echo "❌ Unknown argument: $1. Usage: $0 --raw-fastq/-fq <reads_dir> [--database/-db <database_path>] [-t | --trim] [-r | --remove-host-dna <index_path>]"
          exit 1
          ;;    
    esac
done

# Ensure --raw-fastq is provided
if [[ -z "$RAW_FASTQ_DIR" ]]; then  
    echo "❌ Error: --raw-fastq/-fq is required. Usage: $0 --raw-fastq/-fq <reads_dir> [--database/-db <database_path>] [-t | --trim] [-r | --remove-host-dna <index_path>]"
    exit 1
fi

echo -e "\n======================================================= PIPELINE SUMMARY ========================================================" 

# Display quality control and trimming status
echo -e "Quality Control & Trimming: $( [[ $TRIM == true ]] && echo 'Enabled ✅' || echo 'Disabled ❌' )"
echo "⚠️ Optional, but must be done at least once on the raw reads for the pipeline to work."

echo ""

# Display host DNA removal status and index
if [[ $REMOVE_HOST_DNA == true ]]; then
    echo -e "Host DNA Removal: Enabled ✅"
    echo "Bowtie2 Index Prefix: $BOWTIE_PREFIX"
else
    echo -e "Host DNA Removal: Disabled ❌"
fi
echo "⚠️ Optional, but must be done at least once on the trimmed reads for the pipeline to work."
echo ""

# Display Kraken2/Bracken database path
echo "Kraken2/Bracken Database: $DATABASE"
echo "================================================================================================================================="

echo -e "\n====================================================== PROJECT STRUCTURE ======================================================"
# Define output directories
RESULTS_DIR="$PROJ_DIR/results"  # Stores general results of the entire project
FASTQC_DIR="$RESULTS_DIR/fastqc"  # Stores FastQC reports for quality control analysis
HOST_DNA_ANALYSIS_DIR="$RESULTS_DIR/host_dna_analysis"  # Stores host DNA-related outputs like BLAST, karyotype, etc.
echo "Results are in $RESULTS_DIR (FastQC: $FASTQC_DIR, Host DNA: $HOST_DNA_ANALYSIS_DIR)."

PROCESSED_DIR="$(dirname "$RAW_FASTQ_DIR")/processed_data"  # Root directory for processed data
TRIMMED_DIR="$PROCESSED_DIR/trimmed"  # Stores trimmed FASTQ files

# Define human DNA alignment outputs
# These directories store files that are created once for the entire sub-project, not for each individual run.
HUMAN_DIR="$PROCESSED_DIR/human"  # Stores human DNA-related outputs
ALIGNED_SAM_DIR="$HUMAN_DIR/aligned_sam"  # Stores SAM files of host-aligned reads
SORTED_BAM_DIR="$HUMAN_DIR/sorted_bam"  # Stores sorted BAM files of host-alignments
BED_FILES_DIR="$HUMAN_DIR/bed_files"  # Stores BED files for host DNA reads
tree -d "$HUMAN_DIR" 
echo "Human DNA-related files (SAM, BAM, and BED) in HUMAN_DIR are created once per sub-project and do not change with each pipeline run."

# Define metagenomic Kraken2/Bracken output directories
# Unlike the files in $HUMAN_DIR, these files are generated and overwritten with each pipeline execution.
METAGENOMIC_DIR="$PROCESSED_DIR/metagenomic"  # Stores host-filtered reads and metagenomic classification results
FILTERED_FASTQ_DIR="$METAGENOMIC_DIR/filtered_fastq"  # Stores host-removed paired FASTQ files
KRAKEN2_DIR="$METAGENOMIC_DIR/kraken2"  # Stores Kraken2 classification results
CLASSIFIED_DIR="$METAGENOMIC_DIR/classified"  # Stores classified reads after Kraken2
UNCLASSIFIED_DIR="$METAGENOMIC_DIR/unclassified"  # Stores unclassified reads from Kraken2
BRACKEN_DIR="$METAGENOMIC_DIR/bracken"  # Stores Bracken abundance estimation results
echo "Intermediate metagenomic files in $METAGENOMIC_DIR are generated and overwritten with each pipeline execution and should be inspected beforehand if desired."

# Define run-specific directories
# These directories store final output for each run (e.g., Kraken2 classification, Bracken results, and others).
RUN_DIR="$RESULTS_DIR/runs/$(date +"run_%d%m_%H%M")"  # Unique directory for each pipeline run, based on timestamp
REPORTS_DIR="$RUN_DIR/reports"  # Stores Kraken2 and Bracken reports (k2report, breport, etc.)
KRONA_DIR="$RUN_DIR/krona"  # Stores Krona plots for visualizing taxonomic classification
DIVERSITY_DIR="$RUN_DIR/diversity"  # Stores diversity analysis results (e.g., alpha/beta diversity)
echo "Run-specific results (e.g., Kraken2, Bracken) are in $RUN_DIR."

# Create all necessary directories
mkdir -p "$RESULTS_DIR" "$HOST_DNA_ANALYSIS_DIR" "$FASTQC_DIR/pre_trimming" "$FASTQC_DIR/post_trimming" \
    "$TRIMMED_DIR/paired" "$TRIMMED_DIR/unpaired" "$ALIGNED_SAM_DIR" "$FILTERED_FASTQ_DIR" \
    "$KRAKEN2_DIR" "$CLASSIFIED_DIR" "$UNCLASSIFIED_DIR" "$BRACKEN_DIR" "$REPORTS_DIR" \
    "$KRONA_DIR" "$DIVERSITY_DIR" "$SORTED_BAM_DIR" "$BED_FILES_DIR"

echo "✅ All output directories created successfully."
if [[ "$TRIM" == true ]]; then
    
    echo -e "\n=================================================== QUALITY CONTROL & TRIMMING ==================================================="

    # Run FastQC on all raw reads
    echo -e "\nRunning FastQC on all raw reads..."
    fastqc "$RAW_FASTQ_DIR"/*.fastq.gz --outdir "$FASTQC_DIR/pre_trimming" &>/dev/null || {
        echo "❌ FastQC failed!"; exit 1;
    }
    echo "✅ FastQC completed successfully."

    # Run MultiQC to summarize FastQC reports and clean up zip files
    multiqc "$FASTQC_DIR/pre_trimming" --no-data-dir -o "$FASTQC_DIR/pre_trimming" --force 2>&1 || {
        echo "❌ MultiQC failed!"; exit 1;
    }
    echo -e "✅ MultiQC report generated successfully.\n"
    rm -f "$FASTQC_DIR/pre_trimming"/*.zip  
    
    # Define Trimmomatic adapter file, trimming parameters, and steps
    ADAPT_FA="$ROOT_DIR/data/adapters/zymobiomics_adaptors.fa"
    TRIM_PARAM="PE -phred33 -threads 4"
    STEPS="ILLUMINACLIP:$ADAPT_FA:2:30:10:8:true HEADCROP:5 LEADING:5 CROP:240 SLIDINGWINDOW:4:15 TRAILING:10 MINLEN:35"
## EVENTUALLY, LOOK AT WHICH MINLEN IS NECESSARY FOR ACCURATE CLASSIFICATION WITH KRAKEN2
    echo -e "Running Trimmomatic..."
    for R1 in "$RAW_FASTQ_DIR"//*_R1*.fastq.gz; do # Handle variations of sample names 
        [[ -f "$R1" ]] || { 
          echo "❌ Error: No matching files found in $RAW_FASTQ_DIR"; exit 1; 
        }
        
        # Extract the base sample name by removing lane, read, and other suffixes
        base=$(basename "$R1" | sed -E 's/_L[0-9]+_R[12]_?[0-9]*\.fastq\.gz$//; s/_R[12]_?[0-9]*\.fastq\.gz$//')
        R2=$(echo "$R1" | sed 's/_R1/_R2/')     # Identify the corresponding reverse read

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
        echo "❌ FastQC failed!"; exit 1;
    }
    echo "✅ FastQC completed successfully."
  
    # Run MultiQC to summarize post-trimming FastQC results and clean up
    multiqc "$FASTQC_DIR/post_trimming" --no-data-dir -o "$FASTQC_DIR/post_trimming" --force 2>&1 || {
        echo "❌ MultiQC failed!"; exit 1;
    }
    echo -e "✅ MultiQC report generated successfully.\n"
    rm -f "$FASTQC_DIR/post_trimming"/*.zip  
fi 

echo -e "\n================================================= METAGENOMIC ABUNDANCE ESTIMATION ================================================="

echo "⚠️ This step assumes reads have already been trimmed with Trimmomatic and host DNA removed with Bowtie2. Processed reads must be in the correct directory."

# Start timer for metagenomic classification
START_TIME=$(date +%s)

for R1 in "$TRIMMED_DIR"/paired/*_R1_paired.fastq.gz; do  
    base=$(basename "$R1" "_R1_paired.fastq.gz") 
    echo -e "\nProcessing sample: $base"

    # Host DNA removal with Bowtie2
    if [[ "$REMOVE_HOST_DNA" == true ]]; then
        echo -e "\nRemoving human reads with Bowtie2..."
        bowtie2 -x "$BOWTIE_PREFIX" -p 8 -q --end-to-end --very-sensitive --no-mixed --no-discordant \
            -1 "$TRIMMED_DIR/paired/${base}_R1_paired.fastq.gz" -2 "$TRIMMED_DIR/paired/${base}_R2_paired.fastq.gz" \
            --un-conc "$FILTERED_FASTQ_DIR/${base}_metagenomic" -S "$ALIGNED_SAM_DIR/${base}_human.sam" 2>&1 
        
        # Compress filtered metagenomic reads
        for i in 1 2; do # Compress and rename metagenomic reads
            gzip -f "$FILTERED_FASTQ_DIR/${base}_metagenomic.$i"
        done &&
        echo "✅ Host reads removed and filtered reads compressed."
    fi
    
    # Taxonomic classification with Kraken2
    echo -e "\nClassifying metagenomic reads with Kraken2..."
    
    KRAKEN_CMD="kraken2 --db \"$DATABASE\" --threads 8 --report \"$REPORTS_DIR/${base}.k2report\" \
               --report-minimizer-data --paired --minimum-hit-groups 2 --gzip-compressed \
               --classified-out \"$CLASSIFIED_DIR/${base}_classified#.fastq\" --unclassified-out \"$UNCLASSIFIED_DIR/${base}_unclassified#.fastq\" \
               \"$FILTERED_FASTQ_DIR/${base}_metagenomic.1.gz\" \"$FILTERED_FASTQ_DIR/${base}_metagenomic.2.gz\" \
               --output \"$KRAKEN2_DIR/${base}.kraken2\" --use-names 2>&1"
    
    echo "$KRAKEN_CMD" 
    eval $KRAKEN_CMD
    
    echo "✅ Classification complete."

    # Abundance estimation with Bracken
    echo -e "\nEstimating species abundance with Bracken..."
    bracken -d "$DATABASE" -i "$REPORTS_DIR/${base}.k2report" -l S \
            -r 100 -t 10 -w "$REPORTS_DIR/${base}.breport" -o "$BRACKEN_DIR/${base}.bracken" 2>&1 
    echo "✅ Abundance estimated."
    # Generate Krona interactive plot
    echo -e "\nGenerating Krona visualization..."
    python "$ROOT_DIR/tools/KrakenTools/kreport2krona.py" -r "$REPORTS_DIR/${base}.breport" -o "$KRONA_DIR/${base}.krona.txt" --no-intermediate-ranks && \
    "$ROOT_DIR/tools/Krona/KronaTools/scripts/ImportText.pl" "$KRONA_DIR/${base}.krona.txt" -o "$KRONA_DIR/${base}.krona.html" && rm "$KRONA_DIR/${base}.krona.txt" && \
    echo "✅ Krona plot generated."
done

# End timer for metagenomic classification
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

# Format the duration into HH:MM:SS
printf -v FORMATTED_DURATION '%02d:%02d:%02d' $((DURATION/3600)) $(( (DURATION%3600)/60 )) $((DURATION%60))

echo -e "\nMetagenomic classification completed in: $FORMATTED_DURATION"

echo -e "\n================================================= TAXONOMY HEATMAP ================================================="

# GENERATE THE TAXONOMY HEATMAP 

echo -e "\n================================================= METAGENOMIC DIVERSITY ANALYSIS ================================================="

# Define alpha diversity metrics
METRICS=("BP" "Sh" "F" "Si" "ISi") # Shannon's, Berger-Parker's, Simpson's, Inverse Simpson's, Fisher's

# Create output file with header
echo -e "Sample\t${METRICS[*]}" | tr ' ' '\t' > "$DIVERSITY_DIR/alpha_diversity.tsv"

echo -e "Calculating alpha and beta diversity..."

INPUT_FILES=()  # Initialize input file array for beta diversity

# Loop over all Bracken files to compute diversity
for file in "$BRACKEN_DIR"/*.bracken; do  
    [[ -f "$file" ]] || { echo "❌ Error: No Bracken files found in $BRACKEN_DIR"; exit 1; }
    
    base=$(basename "$file" ".bracken")
    DIVERSITY_RESULTS=("$base")  # Start with sample name
    
    # Compute each alpha diversity metric
    for METRIC in "${METRICS[@]}"; do  
        VALUE=$(python "$ROOT_DIR/tools/KrakenTools/DiversityTools/alpha_diversity.py" -f "$file" -a "$METRIC" | awk -F': ' '{if (NF>1) print $2}')
        DIVERSITY_RESULTS+=("$VALUE")
    done  

    # Save results to TSV
    echo -e "${DIVERSITY_RESULTS[*]}" | tr ' ' '\t' >> "$DIVERSITY_DIR/alpha_diversity.tsv"

    # Collect files for beta diversity
    INPUT_FILES+=("$file")
done  

echo "✅ Alpha diversity calculated for all samples."

# Run beta diversity analysis (Bray-Curtis dissimilarity beta diversity)
if (( ${#INPUT_FILES[@]} < 2 )); then
  echo "⚠️ Skipping beta diversity – fewer than 2 samples."
else
    python "$ROOT_DIR/tools/KrakenTools/DiversityTools/beta_diversity.py" -i "${INPUT_FILES[@]}" --type bracken > "$DIVERSITY_DIR/beta_diversity_matrix.tsv"

    # Generate beta diversity heatmap
    Rscript "$ROOT_DIR/homemade_scripts/b_diversity_heatmap.R" "$DIVERSITY_DIR/beta_diversity_matrix.tsv"
    echo "✅ Beta diversity heatmap generated."

fi

echo -e "\n==================================================== TOTAL READ EXTRACTION ===================================================="

LOG_FILE="$ROOT_DIR/homemade_scripts/logs/kraken_pipeline.log"
TOTAL_READS="$METAGENOMIC_DIR/total_reads.csv"

# Extract total read counts from the log file for normalization
awk '
  BEGIN { print "sample,total_reads" }
  /Processing sample:/ {sample = $NF}  # Capture the sample name
  /sequences \(/ && sample {
    print sample "," $1  # Print sample name and number of reads
    sample = ""  # Reset sample after use
  }
' "$LOG_FILE" > "$TOTAL_READS"
echo "✅ Total read counts saved to: $TOTAL_READS"

# Proceed if host DNA removal was performed
if [[ "$REMOVE_HOST_DNA" == true ]]; then    
    echo -e "\n====================================================== HUMAN DNA ANALYSIS ======================================================"
    
    # Convert and process each SAM file
    for file in "$ALIGNED_SAM_DIR"/*_human.sam; do  
        [[ -f "$file" ]] || { echo "❌ No matching SAM files in $ALIGNED_SAM_DIR"; exit 1; }
        base=$(basename "$file" "_human.sam")
       
        echo -e "\nProcessing host reads from sample: $base"
        # Convert SAM to sorted BAM (excluding unmapped reads)
        samtools view -b -h -F 4 "$file" | samtools sort -@8 -o "$SORTED_BAM_DIR/${base}_human.sorted.bam"
        echo "✅ BAM sorted."
        
        samtools index "$SORTED_BAM_DIR/${base}_human.sorted.bam"
        echo "✅ BAM indexed."
    
        # Convert BAM to sorted BED
        bedtools bamtobed -i "$SORTED_BAM_DIR/${base}_human.sorted.bam" | bedtools sort -i > "$BED_FILES_DIR/${base}.sorted.bed"
        echo "✅ BED created."
    done
  
    echo -e "\n=============================================== KARYOTYPE ================================================"
    
    echo "Identifying overlapping regions across samples..."
    mapfile -t SORTED_BEDS < <(find "$BED_FILES_DIR" -name "*.sorted.bed")

    if (( ${#SORTED_BEDS[@]} < 2 )); then
    	echo "⚠️ Only one BED file found – skipping multiinter and using it directly."
    	cp "${SORTED_BEDS[0]}" "$HOST_DNA_ANALYSIS_DIR/common_intervals.bed"
    else
	bedtools multiinter -header -i "${SORTED_BEDS[@]}" > "$HOST_DNA_ANALYSIS_DIR/common_intervals.bed"
    fi

    echo "Generating karyotype plot..."
    Rscript "$ROOT_DIR/homemade_scripts/karyotype.R" "$HOST_DNA_ANALYSIS_DIR/common_intervals.bed" && 
    echo "✅ Karyotype plot generated."
    
    echo -e "\n================================================== BLAST =================================================="
    
    # Prepare BLAST query file
    BLAST_QUERY="$HOST_DNA_ANALYSIS_DIR/blast_query.fasta" 
    > "$BLAST_QUERY"
  
    echo -e "\nExtracting sequences for BLAST query..."
    while IFS=$'\t' read -r chrom start end _; do
        for FILE in "$SORTED_BAM_DIR"/*.sorted.bam; do
            samtools view "$FILE" "$chrom:$start-$end" | awk '{print ">" $1 "\n" $10}' >> "$BLAST_QUERY"
        done
    done < <(tail -n +2 "$BED_FILES_DIR/common_intervals.bed")
    echo "✅ Sequences saved to BLAST query."

    # Skip BLAST if region file is empty
    if [ $(tail -n +2 "$BED_FILES_DIR/common_intervals.bed" | wc -l) -eq 0 ]; then
        echo "❌ No regions found. Skipping BLAST."
    else
        # Deduplicate sequences
        awk '
            NR % 2 == 1 { id = $0; next }
            { seq = $0; pair = id "\n" seq }
            !seen[pair]++ { print id; print seq }
        ' "$BLAST_QUERY" > "$BLAST_QUERY.tmp" && mv "$BLAST_QUERY.tmp" "$BLAST_QUERY"
        echo "✅ Duplicate sequences removed."
        
        # Run BLAST if query size is reasonable
        if [ $(stat -c%s "$BLAST_QUERY") -lt $((100 * 1024)) ]; then
            echo "Running BLAST search..."
            blastn -query "$BLAST_QUERY" -db nt -out "$HOST_DNA_ANALYSIS_DIR/combined_blast_results.txt" \
                   -evalue 1e-5 -max_target_seqs 5 -outfmt "6 qseqid staxids pident evalue bitscore" -remote 
            echo "✅ BLAST completed."

            echo "Generating taxonomy tree..."
            Rscript "$ROOT_DIR/homemade_scripts/human_aligned_tree.R" "$HOST_DNA_ANALYSIS_DIR/combined_blast_results.txt" && 
            echo "✅ Taxonomy tree generated."
        else
            echo "❌ Query file >100KB. Skipping BLAST. Use: https://blast.ncbi.nlm.nih.gov/Blast.cgi"
        fi
    fi 
    
    echo -e "\n=========================================== JACCARD SIMILARITY ==========================================="
  
    if (( ${#SORTED_BEDS[@]} < 2 )); then
    echo "⚠️ Skipping Jaccard similarity – only one sample found."
    else
        # Initialize output
        rm -f "$HOST_DNA_ANALYSIS_DIR/jaccard_results.txt"
        echo -e "Sample1\tSample2\tIntersection\tUnion\tJaccard\tN_Intersections" > "$HOST_DNA_ANALYSIS_DIR/jaccard_results.txt"
        
        # Compare each unique BED file pair
        for ((i=0; i<${#SORTED_BEDS[@]}; i++)); do
            for ((j=i+1; j<${#SORTED_BEDS[@]}; j++)); do
                FILE1="${SORTED_BEDS[i]}"
                FILE2="${SORTED_BEDS[j]}"
                echo -e "$(basename "$FILE1" .sorted.bed)\t$(basename "$FILE2" .sorted.bed)\t$(bedtools jaccard -a "$FILE1" -b "$FILE2" | tail -n1)" \
                    >> "$HOST_DNA_ANALYSIS_DIR/jaccard_results.txt"
            done
        done
        
        echo "Generating Jaccard similarity heatmap..."
        Rscript "$ROOT_DIR/homemade_scripts/jaccard_similarity_heatmap.R" "$HOST_DNA_ANALYSIS_DIR/jaccard_results.txt"
        echo "✅ Jaccard analysis complete."
    fi
fi 

echo -e "\n✅ Pipeline completed successfully."

echo "Storing log file..."
cp -r "$RAW_FASTQ_DIR/../../homemade_scripts/logs" "$RUN_DIR"

echo "SLURM JOB ID: $SLURM_JOB_ID"
