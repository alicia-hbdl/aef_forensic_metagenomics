#!/bin/bash

#SBATCH --job-name=pipeline
#SBATCH --output=logs/pipeline_%j.log
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=200G  
#SBATCH --cpus-per-task=8 
#SBATCH --ntasks=1  
#SBATCH --partition=cpu 

# TO DO: match the usage with the mac mini and requirements + allow user to specify adaptor file
# maybe don't need to create fastqc directories here since they are created in the sub-script
set -e # Exit on error 
set -x  # Print each command and its arguments as it is executed for debugging 

<< 'COMMENT'
=============================================================================================
Usage:
    ./pipeline.sh --raw-fastq/-f <reads_dir> [--database/-d <database_path>] [-t | --trim] [-r | --remove-host-dna <index_path>]

Arguments:
    -f, --raw-fastq <path>       Path to the directory containing raw FASTQ reads to process.
    -d, --database <path>        (Optional) Path to the Kraken2/Bracken database used for metagenomic classification.
    -t, --trim                    (Optional) Enable adapter and quality trimming of reads.
    -r, --remove-host-dna <path>  (Optional) Remove host DNA contamination.

Notes:
    - Ensure that all necessary dependencies (Kraken2, Trimmomatic, etc.) are installed before running.
=============================================================================================
COMMENT

echo -e "\n================================================= ARGUMENT PARSING & VALIDATION =================================================\n"

# Default values
TRIM=false
REMOVE_HOST_DNA=false

# Parse command-line arguments
while [[ $# -gt 0 ]]; do 
    case "$1" in
    # Process raw FASTQ directory
        -f|--raw-fastq) 
            # Validate FASTQ directory and its contents
            if [[ -z "$2" || "$2" == -* || ! -d "$2" || -z "$(find "$2" -maxdepth 1 -type f \( -name "*.fq" -o -name "*.fastq" -o -name "*.fq.gz" -o -name "*.fastq.gz" \) 2>/dev/null)" ]]; then
              echo "❌  "$2" is an invalid/missing FASTQ directory. Ensure it contains .fq, .fastq, .fq.gz, or .fastq.gz files."
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
            DATABASE="$ROOT_DIR/data/databases/k2_standard_16gb_20250402"
            BOWTIE_PREFIX="$ROOT_DIR/data/bowtie_index/GRCh38_noalt_as/GRCh38_noalt_as"           
            
            shift 2
            ;;
    
        # Process Kraken2/Bracken database directory argument
        -d|--database) 
            if [[ -z "$2" || "$2" == -* || ! -d "$2" || ! -f "$2/hash.k2d" || ! -f "$2/taxo.k2d" || ! -f "$2/opts.k2d" || ! $(ls "$2"/*kmer_distrib 2>/dev/null) ]]; then
                echo "⚠️  "$2" is an invalid/missing database. Using default: $DATABASE"
            else
                DATABASE="$2"
                echo "✅  Database set to: $DATABASE"
            fi
            shift 2
            ;;
            
      # Enable trimming
      -t|--trim)
          TRIM=true
          echo "✅  Trimming enabled."
          shift
          ;;            
        
    # Process host DNA Bowtie2 index prefix
    -r|--remove-host-dna)
        if [[ -z "$2" || "$2" == -* ]]; then
            echo "⚠️  No Bowtie2 index provided. Using default: $BOWTIE_PREFIX"
            shift 
        elif [[ ! -f "$2.1.bt2" || ! -f "$2.2.bt2" || ! -f "$2.3.bt2" || ! -f "$2.4.bt2" || ! -f "$2.rev.1.bt2" || ! -f "$2.rev.2.bt2" ]]; then
            echo "⚠️  "$2" is an incomplete/ missing Bowtie2 index. Using default: $BOWTIE_PREFIX"
            shift 2
        else
            BOWTIE_PREFIX="$2"
            echo "✅  Bowtie2 index set to: $BOWTIE_PREFIX"
            shift 2
        fi
        REMOVE_HOST_DNA=true
        ;;       

      # Handle unknown arguments
      *) 
          echo "❌ Unknown argument: $1. Usage: $0 --raw-fastq/-f <reads_dir> [--database/-d <database_path>] [-t | --trim] [-r | --remove-host-dna <index_path>]"
          exit 1
          ;;    
    esac
done

# Ensure --raw-fastq is provided
if [[ -z "$RAW_FASTQ_DIR" ]]; then  
    echo "❌  --raw-fastq/-f is required. Usage: $0 --raw-fastq/-f <reads_dir> [--database/-d <database_path>] [-t | --trim] [-r | --remove-host-dna <index_path>]"
    exit 1
fi

# Display quality control and trimming status
echo -e "Quality Control & Trimming: $( [[ $TRIM == true ]] && echo 'Enabled ✅' || echo 'Disabled ❌' )"
echo "⚠️  Optional, but must be done at least once on the raw reads for the pipeline to work."

echo ""

# Display host DNA removal status and index
if [[ $REMOVE_HOST_DNA == true ]]; then
    echo -e "Host DNA Removal: Enabled ✅"
    echo "Bowtie2 Index Prefix: $BOWTIE_PREFIX"
else
    echo -e "Host DNA Removal: Disabled ❌"
fi
echo "⚠️  Optional, but must be done at least once on the trimmed reads for the pipeline to work."
echo ""

# Display Kraken2/Bracken database path
echo "Kraken2/Bracken Database: $DATABASE"

echo -e "\n================================================= CONDA ENVIRONMENT ACTIVATION ==================================================\n"

source "$ROOT_DIR/scripts/helper_scripts/environment_setup.sh" "$ROOT_DIR/scripts/metagenomics.yml" || { echo "❌ Failed to set up Conda environment."; exit 1; }

echo -e "\n====================================================== PROJECT STRUCTURE ======================================================\n"

# Define output directories
RESULTS_DIR="$PROJ_DIR/results"  # General results of the entire project
FASTQC_DIR="$RESULTS_DIR/fastqc"  # FastQC reports for quality control analysis
HOST_DNA_ANALYSIS_DIR="$RESULTS_DIR/host_dna_analysis"  # Host DNA-related outputs (e.g., BLAST, karyotype)

PROCESSED_DIR="$(dirname "$RAW_FASTQ_DIR")/processed_data"  # Root directory for processed data
TRIMMED_DIR="$PROCESSED_DIR/trimmed"  # Trimmed FASTQ files

# Define human DNA alignment outputs
# These files are created once for the entire sub-project and do not change with each pipeline run.
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
mkdir -p "$RESULTS_DIR" "$HOST_DNA_ANALYSIS_DIR" "$FASTQC_DIR/pre_trimming" "$FASTQC_DIR/post_trimming" \
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

    "$ROOT_DIR/scripts/helper_scripts/qc_trim.sh" "$RAW_FASTQ_DIR" || { echo "❌ Quality control and trimming failed!"; exit 1; }
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
		--report-minimizer-data --paired --minimum-hit-groups 2 \
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
# Rscript "$ROOT_DIR/scripts/helper_scripts/phylo_classification_comparison.R" -r "$REPORTS_DIR/../combined_breports.csv" -t "$RAW_FASTQ_DIR/ground_truth.csv" 

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
    python "$ROOT_DIR/scripts/helper_scripts/evaluation_metrics.py" --run-dir "$RESULTS_DIR/runs/" --ground-truth "$RAW_FASTQ_DIR/ground_truth.csv" 2>&1 || { echo "❌ Evaluation metrics generation failed!"; exit 1; }
    Rscript "$ROOT_DIR/scripts/helper_scripts/read_progression.R" "$RESULTS_DIR/" 2>&1 || { echo "❌ Read progression plot generation failed!"; exit 1; }
    
    # Remove the original log file when successful 
    echo -e "\nRemoving original SLURM job logs"
    rm $ROOT_DIR/scripts/logs/*_"$SLURM_JOB_ID".*
else
  echo "❌ Copying log files failed!"
  exit 1
fi
