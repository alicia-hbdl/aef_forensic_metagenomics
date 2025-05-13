#!/bin/bash

#This is a pipeline that processes aligned SAM files, converts them to sorted BAM files, and generates BED files for further analysis.
#It also performs karyotype analysis, BLAST searches, and Jaccard similarity calculations.
#Usage: ./host_dna_analysis.sh <ALIGNED_SAM_DIR>
  
# Ensure the SAM directory is provided and valid
if [[ -z "$1" || "$1" == -* || ! -d "$1" || -z "$(find "$1" -maxdepth 1 -name "*.sam" 2>/dev/null)" ]]; then
    echo "❌ Error: Invalid or missing ALIGNED_SAM directory. Ensure it contains .sam files."
    exit 1
fi

ALIGNED_SAM_DIR="$1"  # Set the ALIGNED_SAM_DIR to the directory passed as argument

# Define and create output directories for processed files and results
# These are created in the main script but are also initialized here for standalone use
SORTED_BAM_DIR="$ALIGNED_SAM_DIR/../sorted_bam"  
BED_FILES_DIR="$ALIGNED_SAM_DIR/../bed_files"  
HOST_DNA_ANALYSIS_DIR="$ALIGNED_SAM_DIR/../../../results/host_dna_analysis"  
mkdir -p "$SORTED_BAM_DIR" "$BED_FILES_DIR" "$HOST_DNA_ANALYSIS_DIR"
ROOT_DIR="$ALIGNED_SAM_DIR/../../../.."

# Convert and process each SAM file
for file in "$ALIGNED_SAM_DIR"/*.sam; do  
    base=$(basename "$file" .sam)  

    echo -e "\nProcessing host reads from sample: $base"
    
    # Convert SAM to sorted BAM (excluding unmapped reads)
    samtools view -b -h -F 4 "$file" | samtools sort -@8 -o "$SORTED_BAM_DIR/${base}.sorted.bam" || { echo "❌ Error sorting BAM for $base"; exit 1; }
    echo "✅ BAM sorted."
    
    samtools index "$SORTED_BAM_DIR/${base}.sorted.bam" || { echo "❌ Error indexing BAM for $base"; exit 1; }
    echo "✅ BAM indexed."
    
    # Convert BAM to sorted BED
    bedtools bamtobed -i "$SORTED_BAM_DIR/${base}.sorted.bam" | bedtools sort -i > "$BED_FILES_DIR/${base}.sorted.bed" || { echo "❌ Error creating BED for $base"; exit 1; }
    echo "✅ BED created."
done
  
echo -e "\n=============================================== KARYOTYPE ================================================"
    
# Identify overlapping regions across samples
echo -e "Finding overlapping regions across samples..."
mapfile -t SORTED_BEDS < <(find "$BED_FILES_DIR" -name "*.sorted.bed") || { echo "❌ Failed to find or read BED files."; exit 1; }

if (( ${#SORTED_BEDS[@]} < 2 )); then
    echo "⚠️ Only one BED file found – skipping multiinter and using it directly."
    cp "${SORTED_BEDS[0]}" "$HOST_DNA_ANALYSIS_DIR/intervals.bed" || { echo "❌ Failed to copy BED file."; exit 1; }
else
    bedtools multiinter -header -i "${SORTED_BEDS[@]}" > "$HOST_DNA_ANALYSIS_DIR/intervals.bed" || { echo "❌ Failed to find overlapping regions."; exit 1; }
    echo "✅ Overlapping regions saved."
fi

echo -e "Generating karyotype plot..."
Rscript "$ROOT_DIR/scripts/helper_scripts/karyotype.R" "$HOST_DNA_ANALYSIS_DIR/intervals.bed" || { echo "❌ Failed to generate karyotype plot."; exit 1; }
echo "✅ Karyotype plot generated."

# Filter regions present in more than one sample
awk '$4 > 1' "$HOST_DNA_ANALYSIS_DIR/intervals.bed" > "$HOST_DNA_ANALYSIS_DIR/common_intervals.bed" || { echo "❌ Failed to filter common intervals."; exit 1; }
echo "✅ Filtered common intervals saved."

echo -e "\n================================================== BLAST =================================================="
    
# Prepare BLAST query file
BLAST_QUERY="$HOST_DNA_ANALYSIS_DIR/blast_query.fasta" 
> "$BLAST_QUERY" || { echo "❌ Failed to create BLAST query file."; exit 1; }

echo -e "\nExtracting sequences for BLAST query..."
while IFS=$'\t' read -r chrom start end _; do
    for FILE in "$SORTED_BAM_DIR"/*.sorted.bam; do
        samtools view "$FILE" "$chrom:$start-$end" | awk '{print ">" $1 "\n" $10}' >> "$BLAST_QUERY" || { echo "❌ Failed to extract sequences from $FILE."; exit 1; }
    done
done < <(tail -n +2 "$BED_FILES_DIR/common_intervals.bed") || { echo "❌ Failed to read BED file."; exit 1; }
echo "✅ Sequences saved to BLAST query."

# Skip BLAST if BED region file is empty
if [ "$(tail -n +2 "$BED_FILES_DIR/common_intervals.bed" | wc -l)" -eq 0 ]; then
    echo "❌ No regions found. Skipping BLAST."
    rm -f "$BLAST_QUERY" || { echo "❌ Failed to remove empty BLAST query file."; exit 1; }
else
    # Deduplicate sequences
    awk '
        NR % 2 == 1 { id = $0; next }
        { seq = $0; pair = id "\n" seq }
        !seen[pair]++ { print id; print seq }
    ' "$BLAST_QUERY" > "$BLAST_QUERY.unique" || { echo "❌ Failed to remove duplicates."; exit 1; }
    echo "✅ Duplicate sequences removed."

    # Randomly select 15 unique sequences
    awk 'BEGIN { RS=">"; ORS=""; srand() }
        NR > 1 { seqs[++n] = $0 }
        END {
            for (i = n; i > 1; i--) {
                j = int(rand() * i) + 1
                tmp = seqs[i]; seqs[i] = seqs[j]; seqs[j] = tmp
            }
            for (i = 1; i <= 15 && i <= n; i++) print ">" seqs[i]
        }' "$BLAST_QUERY.unique" > "$BLAST_QUERY.selected" || { echo "❌ Failed to select random sequences."; exit 1; }
    echo "✅ 15 random sequences selected."

    echo "Running BLAST search..."
    blastn -query "$BLAST_QUERY.selected" -db nt -out "$HOST_DNA_ANALYSIS_DIR/combined_blast_results.txt" \
           -evalue 1e-5 -max_target_seqs 3 -outfmt "6 qseqid staxids pident evalue bitscore" -remote || { echo "❌ BLAST failed."; exit 1; }
    echo "✅ BLAST completed."
    
    # Run taxonomy tree only if ≥ 3 unique tax IDs (column 2) exist
    if [ "$(cut -f2 "$HOST_DNA_ANALYSIS_DIR/combined_blast_results.txt" | sort -u | wc -l)" -ge 3 ]; then
       # Generate taxonomy tree from BLAST results
       Rscript "$ROOT_DIR/scripts/helper_scripts/human_aligned_tree.R" "$HOST_DNA_ANALYSIS_DIR/combined_blast_results.txt" || { echo "❌ Failed to generate taxonomy tree."; exit 1; }
        echo "✅ Taxonomy tree generated."
    else
        echo "⚠️ Skipping taxonomy tree: fewer than 3 unique tax IDs."
        echo "Unique tax IDs found:"
        cut -f2 "$HOST_DNA_ANALYSIS_DIR/combined_blast_results.txt" | sort -u || { echo "❌ Failed to list unique tax IDs."; exit 1; }
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
    Rscript "$ROOT_DIR/scripts/helper_scripts/jaccard_similarity_heatmap.R" "$HOST_DNA_ANALYSIS_DIR/jaccard_results.txt"
    echo "✅ Jaccard analysis complete."
fi
