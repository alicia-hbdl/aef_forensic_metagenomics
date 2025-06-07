#!/bin/bash

#SBATCH --job-name=pipeline
#SBATCH --output=logs/pipeline_%j.log
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=60G  
#SBATCH --cpus-per-task=8 
#SBATCH --ntasks=1  
#SBATCH --partition=cpu 

set -e # Exit on error 
set -x  # Print each command for debugging 

# Project root path (user-specific)
ROOT_DIR="/scratch/users/k24087895/final_project"

<< 'COMMENT'
=============================================================================================
Usage:
    ./pipeline.sh --raw-fastq/-f <reads_dir> 
                  [--database/-d <database_path>] 
                  [-t|--trim [adapter_file.fa]] 
                  [-r|--remove-host-dna <index_prefix>] 
                  [-g|--ground-truth <ground_truth.csv>]

Arguments:
    -f, --raw-fastq <path>            Path to the directory containing raw FASTQ reads to process. (Required)
    -d, --database <path>             Path to the Kraken2/Bracken database for classification. (Optional)
    -t, --trim [adapter_file.fa]      Enable adapter and quality trimming. Optionally provide a FASTA adapter file. (Optional)
    -r, --remove-host-dna [prefix]    Enable host DNA removal. Optionally provide a Bowtie2 index prefix. (Optional)
    -g, --ground-truth <file.csv>     CSV file with 'species,abundance' columns for validation. (Optional)
    
Notes:
    - The FASTQ directory must contain at least one file ending in .fq, .fastq, .fq.gz, or .fastq.gz.
    - Ensure all required tools (Kraken2, Bracken, Bowtie2, Trimmomatic, etc.) are installed and accessible.
=============================================================================================
COMMENT

echo -e "\n================================================= ARGUMENT PARSING & VALIDATION =================================================\n"

# Default flags and file paths (can be overridden by arguments)
TRIM=false
REMOVE_HOST_DNA=false
GT_FLAG=false
DATABASE="$ROOT_DIR/data/databases/k2_standard_16gb_20250402"
BOWTIE_PREFIX="$ROOT_DIR/data/bowtie_index/GRCh38_noalt_as/GRCh38_noalt_as"           
ADAPTER_FILE="$ROOT_DIR/data/adapters/TruSeq3-PE-2.fa"            

# Usage message
print_usage() {
  echo "Usage: $0 --raw-fastq/-f <reads_dir> \
  [--database/-d <database_path>] [-t|--trim [adapter_file.fa]] \
  [-r|--remove-host-dna [index_prefix]] [-g|--ground-truth <file.csv>]"
  exit 1
}

# Exit immediately if no arguments are provided
[[ $# -eq 0 ]] && { echo "❌ No arguments provided."; print_usage; }

# Parse command-line arguments
while [[ $# -gt 0 ]]; do 
    case "$1" in
    
        -f|--raw-fastq) # Process raw FASTQ directory
            # Validate input directory and FASTQ file presence
            if [[ -z "$2" || "$2" == -* || ! -d "$2" || -z "$(find "$2" -maxdepth 1 -type f \( -name "*.fq" -o -name "*.fastq" -o -name "*.fq.gz" -o -name "*.fastq.gz" \) 2>/dev/null)" ]]; then
              echo "❌  '$2' is an invalid/missing FASTQ directory. Ensure it contains .fq, .fastq, .fq.gz, or .fastq.gz files."
              print_usage
            fi
            
            # Compress and rename FASTQ files to .fastq.gz
            for file in "$2"/*.{fq,fastq}; do [[ -f "$file" ]] && gzip "$file" && echo "Compressing '$file'..."; done     
            for file in "$2"/*.fq.gz; do [[ -f "$file" ]] && mv "$file" "${file%.*}.fastq.gz" && echo "Renaming '$file'..."; done 

            # Set project paths
            RAW_FASTQ_DIR="$2" 
            PROJ_DIR="$(dirname "$RAW_FASTQ_DIR")"  # Current subproject directory
            shift 2
            ;;
    
        -d|--database) # Process Kraken2/Bracken database directory argument (Optional)
            if [[ -z "$2" || "$2" == -* ]]; then
                echo "⚠️  No database provided. Using default."
                shift 
            elif [[ ! -d "$2" || ! -f "$2/hash.k2d" || ! -f "$2/taxo.k2d" || ! -f "$2/opts.k2d" ]] || ! compgen -G "$2"/*kmer_distrib > /dev/null; then
              echo "⚠️ '$2' is an invalid database. Using default."
              shift 2
            else
                DATABASE="$2"
                shift 2
            fi
            ;;

        -t|--trim) # Enable trimming and set adapter file (Optional)
          TRIM=true
          if [[ -z "$2" || "$2" == -* ]]; then
              echo "⚠️  No adapter file provided. Using default."
              shift
          elif [[ -f "$2" ]] && grep -q "^>" "$2"; then
              ADAPTER_FILE="$2"
              shift 2
          else
            echo "⚠️  '$2' is not a valid FASTA adapter file. Using default."
            shift 2
          fi
          ;;    
          
      -r|--remove-host-dna) # Process host DNA Bowtie2 index prefix  (Optional)
          REMOVE_HOST_DNA=true
          if [[ -z "$2" || "$2" == -* ]]; then
              echo "⚠️  No Bowtie2 index provided. Using default."
              shift 
          elif [[ ! -f "$2.1.bt2" || ! -f "$2.2.bt2" || ! -f "$2.3.bt2" || ! -f "$2.4.bt2" || ! -f "$2.rev.1.bt2" || ! -f "$2.rev.2.bt2" ]]; then
              echo "⚠️  '$2' is an incomplete/ missing Bowtie2 index. Using default."
              shift 2
          else
              BOWTIE_PREFIX="$2"
              shift 2
          fi
          ;;
  
      -g|--ground-truth)  # Ground truth file  (Optional)
        if [[ -z "$2" || "$2" == -* ]]; then
          echo "⚠️  No ground truth file provided. Skipping."
          shift
        elif [[ -f "$2" && "$2" == *.csv ]]; then
          GROUND_TRUTH="$2"
          GT_FLAG=true
          shift 2
        else
          echo "⚠️  '$2' is an invalid or incorrectly formatted ground truth file. Skipping."
          shift 2
        fi
        ;;
  
      # Handle unknown arguments
      *) 
      echo "❌ Unknown argument '$2'"; print_usage        
      ;;   
    esac
done

# Ensure raw FASTQ directory is provided
[[ -z "$RAW_FASTQ_DIR" ]] && { echo "❌  No raw FASTQ directory provided"; print_usage; }

# TODO: Validate the directories and exit 
if  [[ $TRIM == false || $REMOVE_HOST_DNA == false ]]; then
echo "⚠️  Trimming and host DNA removal are optional but must be done once; if skipped, ensure trimmed and filtered files exist in the correct directories.\n"
fi 

# Pipeline configuration summary
echo "===== Pipeline Configuration Summary ====="
echo "Raw FASTQ directory path: $RAW_FASTQ_DIR" 
echo "Kraken2/Bracken Database: $DATABASE"
echo "Ground Truth: $( [[ $GT_FLAG == true ]] && echo "$GROUND_TRUTH ✅" || echo "Not provided ❌" )"
echo "Quality Control & Trimming: $( [[ $TRIM == true ]] && echo "Enabled ✅ (Adapter: $ADAPTER_FILE)" || echo "Disabled ❌" )"
echo "Host DNA Removal: $( [[ $REMOVE_HOST_DNA == true ]] && echo "Enabled ✅ (Index: $BOWTIE_PREFIX)" || echo "Disabled ❌" )"

echo -e "\n================================================= CONDA ENVIRONMENT ACTIVATION ==================================================\n"

source "$ROOT_DIR/scripts/helper_scripts/environment_setup.sh" "$ROOT_DIR/scripts/metagenomics.yml" || { echo "❌ Failed to set up Conda environment."; exit 1; }

echo -e "\n====================================================== PROJECT STRUCTURE ======================================================\n"

# Define output directories
RESULTS_DIR="$PROJ_DIR/results"  # General results of the entire project
HOST_DNA_ANALYSIS_DIR="$RESULTS_DIR/host_dna_analysis"  # Host DNA-related outputs (e.g., BLAST, karyotype)

PROCESSED_DIR="$(dirname "$RAW_FASTQ_DIR")/processed_data"  # Root directory for processed data
TRIMMED_DIR="$PROCESSED_DIR/trimmed"  # Trimmed FASTQ files

# Define human DNA alignment outputs
# Created once per sub-project; only updated when -t or -r is used
HUMAN_DIR="$PROCESSED_DIR/human"  # Human DNA-related outputs
ALIGNED_SAM_DIR="$HUMAN_DIR/aligned_sam"  # SAM files of host-aligned reads
SORTED_BAM_DIR="$HUMAN_DIR/sorted_bam"  # Sorted BAM files of host-alignments
BED_FILES_DIR="$HUMAN_DIR/bed_files"  # BED files for host DNA reads

# Define metagenomic Kraken2/Bracken output directories
# These files are overwritten with each pipeline execution.
METAGENOMIC_DIR="$PROCESSED_DIR/metagenomic"  # Host-filtered reads and metagenomic classification results
FILTERED_FASTQ_DIR="$METAGENOMIC_DIR/filtered_fastq"  # Host-removed paired FASTQ files
KRAKEN2_DIR="$METAGENOMIC_DIR/kraken2"  # Kraken2 classification results
CLASSIFIED_DIR="$METAGENOMIC_DIR/classified"  # Classified reads after Kraken2
UNCLASSIFIED_DIR="$METAGENOMIC_DIR/unclassified"  # Unclassified reads from Kraken2
BRACKEN_DIR="$METAGENOMIC_DIR/bracken"  # Bracken abundance estimation results

# Define run-specific directories
# These directories store final output for each run (e.g., Kraken2, Bracken results).
RUN_DIR="$RESULTS_DIR/runs/$(date +"run_%d%m_%H%M")"  # Unique directory for each pipeline run
REPORTS_DIR="$RUN_DIR/reports"  # Kraken2 and Bracken reports (k2report, breport, etc.)
KRONA_DIR="$RUN_DIR/krona"  # Krona plots for taxonomic classification
DIVERSITY_DIR="$RUN_DIR/diversity"  # Diversity analysis results (e.g., alpha/beta diversity)
LOG_DIR="$RUN_DIR/logs"  # Log files for the pipeline run

# Create all necessary directories
mkdir -p "$RESULTS_DIR" "$HOST_DNA_ANALYSIS_DIR" \
    "$TRIMMED_DIR/paired" "$TRIMMED_DIR/unpaired" "$ALIGNED_SAM_DIR" "$FILTERED_FASTQ_DIR" \
    "$KRAKEN2_DIR" "$CLASSIFIED_DIR" "$UNCLASSIFIED_DIR" "$BRACKEN_DIR" "$REPORTS_DIR" \
    "$KRONA_DIR" "$DIVERSITY_DIR" "$LOG_DIR" "$SORTED_BAM_DIR" "$BED_FILES_DIR"

# Display the directory structure of the project
tree -L 3 -d "$PROJ_DIR"

# Inform the user about the project structure
echo "Human DNA-related files (SAM, BAM, and BED) in the 'human' directory are created once per sub-project and remain unchanged across runs."
echo "Intermediate metagenomic files in the 'metagenomic' directory are overwritten with each run and should be inspected beforehand if desired."
echo "Results are stored in the 'results' directory (including FastQC, host DNA analysis, and run-specific results like Kraken2 and Bracken)."

if [[ "$TRIM" == true ]]; then
    echo -e "\n=================================================== QUALITY CONTROL & TRIMMING ===================================================\n"

    "$ROOT_DIR/scripts/helper_scripts/qc_trim.sh" --reads "$RAW_FASTQ_DIR" --adapters "$ADAPTER_FILE" || { echo "❌ Quality control and trimming failed!"; exit 1; }
    echo -e "✅  QC and trimming completed successfully."

    # Calculate read length stats (raw and trimmed)
    echo -e "\nCalculating read length statistics..."
    files=($RAW_FASTQ_DIR/*.fastq.gz $TRIMMED_DIR/paired/*.fastq.gz) # Expand the file patterns 
    "$ROOT_DIR/scripts/helper_scripts/get_fastq_stats.sh" "${files[@]}" 2>&1 || { echo "❌ Failed to generate read length statistics."; exit 1; }
    echo "✅ Statistics generated."
fi

echo -e "\n================================================= METAGENOMIC ABUNDANCE ESTIMATION ================================================="

echo "⚠️  This step assumes reads have already been trimmed with Trimmomatic and host DNA removed with Bowtie2. Processed reads must be in the correct directory."

# Start timer for metagenomic classification
START_TIME=$(date +%s)

for R1 in "$TRIMMED_DIR"/paired/*_R1_paired.fastq.gz; do  
    base=$(basename "$R1" "_R1_paired.fastq.gz") 
    echo -e "\nProcessing sample: $base"

    # Host DNA removal with Bowtie2
    if [[ "$REMOVE_HOST_DNA" == true ]]; then
        echo -e "\nRemoving human reads with Bowtie2..."
        BOWTIE_CMD="bowtie2 -x \"$BOWTIE_PREFIX\" -p 8 -q --end-to-end --very-sensitive --no-mixed --no-discordant \
	                -1 \"$TRIMMED_DIR/paired/${base}_R1_paired.fastq.gz\" -2 \"$TRIMMED_DIR/paired/${base}_R2_paired.fastq.gz\" \
	                --un-conc \"$FILTERED_FASTQ_DIR/${base}_metagenomic\" -S \"$ALIGNED_SAM_DIR/${base}_human.sam\" 2>&1"
        echo "$BOWTIE_CMD"
        eval "$BOWTIE_CMD"
        echo "✅  Host reads removed."
    fi
    
    # Taxonomic classification with Kraken2
    # Note: Argument order matters—input files must come last, or it may cause a segmentation fault.
    echo -e "\nClassifying metagenomic reads with Kraken2..."
    KRAKEN_CMD="kraken2 --db \"$DATABASE\" --threads 8 --report \"$REPORTS_DIR/${base}.k2report\" \
		--report-minimizer-data --paired --minimum-hit-groups 10 \
		--classified-out \"$CLASSIFIED_DIR/${base}_classified#.fastq\" --unclassified-out \"$UNCLASSIFIED_DIR/${base}_unclassified#.fastq\" \
		--output \"$KRAKEN2_DIR/${base}.kraken2\" --use-names \"$FILTERED_FASTQ_DIR/${base}_metagenomic.1\" \"$FILTERED_FASTQ_DIR/${base}_metagenomic.2\" 2>&1"
    echo "$KRAKEN_CMD"
    eval "$KRAKEN_CMD"
    echo "✅  Classification complete."

    # Abundance estimation with Bracken
    echo -e "\nEstimating species abundance with Bracken..."
    bracken -d "$DATABASE" -i "$REPORTS_DIR/${base}.k2report" -l S \
            -r 100 -t 0 -w "$REPORTS_DIR/${base}.breport" -o "$BRACKEN_DIR/${base}.bracken" 2>&1 
    echo "✅  Abundance estimated."
    
    # Generate Krona interactive plot
    echo -e "\nGenerating Krona visualization..."
    python "$ROOT_DIR/tools/KrakenTools/kreport2krona.py" -r "$REPORTS_DIR/${base}.breport" -o "$KRONA_DIR/${base}.krona.txt" --no-intermediate-ranks && \
    "$ROOT_DIR/tools/Krona/KronaTools/scripts/ImportText.pl" "$KRONA_DIR/${base}.krona.txt" -o "$KRONA_DIR/${base}.krona.html" && rm "$KRONA_DIR/${base}.krona.txt" && \
    echo "✅  Krona plot generated."
done

# Record and display total runtime
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
printf -v FORMATTED_DURATION '%02d:%02d:%02d' $((DURATION/3600)) $(( (DURATION%3600)/60 )) $((DURATION%60)) # Format runtime as HH:MM:SS
echo -e "\nMetagenomic classification completed in: $FORMATTED_DURATION"

# Calculate read length stats (filtered, classified, unclassified FASTQ files)
echo -e "\nCalculating read length statistics..."
files=("$FILTERED_FASTQ_DIR"/*.1 "$FILTERED_FASTQ_DIR"/*.2 "$CLASSIFIED_DIR"/*.fastq "$UNCLASSIFIED_DIR"/*.fastq) # Expand file patterns into actual files
"$ROOT_DIR/scripts/helper_scripts/get_fastq_stats.sh" "${files[@]}" 2>&1 || { echo "❌ Failed to generate read length statistics."; exit 1; }
echo "✅ Statistics generated."

# Combine all Bracken reports
echo -e "\nCombining Bracken reports..."
Rscript "$ROOT_DIR/scripts/helper_scripts/combine_breports.R" "$REPORTS_DIR"/*.breport || { echo "❌ Failed to combine Bracken reports."; exit 1; }
echo "✅ Bracken reports combined."

# Copy SLURM job logs to output directory
echo -e "\nCopying SLURM job logs..."
cp "$ROOT_DIR"/scripts/logs/*_"$SLURM_JOB_ID".* "$LOG_DIR" || { echo "❌ Failed to copy SLURM job logs."; exit 1; }
echo "✅ SLURM job logs copied."

echo -e "\n================================================= COMPARISON TO GROUND TRUTH ================================================="

# Generate the heatmap that looks at how the run differs from the groundn truth 
# Rscript "$ROOT_DIR/scripts/helper_scripts/phylo_classification_comparison.R" -r "$(dirname "$REPORTS_DIR")/combined_breports.csv" -t "$GROUND_TRUTH" 

echo -e "\n================================================= METAGENOMIC DIVERSITY ANALYSIS ================================================="

"$ROOT_DIR/scripts/helper_scripts/diversity_analysis.sh" -b  "$BRACKEN_DIR" -d "$DIVERSITY_DIR" 2>&1 || { echo "❌ Diversity analysis failed!"; exit 1; }
echo -e "✅  Diversity analysis completed successfully."

# Proceed if host DNA removal was performed
if [[ "$REMOVE_HOST_DNA" == true ]]; then    
    echo -e "\n====================================================== HUMAN DNA ANALYSIS ======================================================"
    
    "$ROOT_DIR/scripts/helper_scripts/host_dna_analysis.sh" "$ALIGNED_SAM_DIR" 2>&1 || { echo "❌ Host DNA analysis failed!"; exit 1; }
    echo -e "✅  Host DNA analysis completed successfully."
fi 

echo -e "\n✅  Pipeline completed successfully."

# Move final SLURM logs to run directory and clean up originals
echo -e "\nCopying final SLURM job logs"
if cp "$ROOT_DIR"/scripts/logs/*_"$SLURM_JOB_ID".* "$LOG_DIR"; then

    echo -e "\n====================================================== COMPARISON TO OTHER RUNS ======================================================"
  
    # Generate the summary table and evaluation metrics including this run 
    "$ROOT_DIR/scripts/helper_scripts/runs_summary.sh" "$RESULTS_DIR/runs/" 2>&1 || { echo "❌ Summary generation failed!"; exit 1; }
    
    if [[ "$GT_FLAG" == true ]]; then # Check if ground truth file is set and exists
        BREPORT_FILES=("$RESULTS_DIR"/runs/*/combined_breports.csv)
        if [[ ${#BREPORT_FILES[@]} -gt 0 ]]; then # Then check if breport files exist
            Rscript "$ROOT_DIR/scripts/helper_scripts/downstream_analysis.R" -t "$GROUND_TRUTH" -s "$RESULTS_DIR/runs/runs_summary.csv" \
              "${BREPORT_FILES[@]}" || { echo "❌ Downstream analysis failed."; exit 1; }
        else
            echo "⚠️  No breport files found. Skipping downstream analysis."
        fi
    else
        echo "⚠️  No ground truth provided. Skipping downstream analysis."
    fi    
    
    # Remove the original log file when successful 
    echo -e "\nRemoving original SLURM job logs"
    rm $ROOT_DIR/scripts/logs/*_"$SLURM_JOB_ID".*
else
  echo "❌ Copying log files failed!"
  exit 1
fi
