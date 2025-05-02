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
    samtools view -b -h -F 4 "$file" | samtools sort -@8 -o "$SORTED_BAM_DIR/${base}.sorted.bam"
    echo "✅ BAM sorted."
    
    samtools index "$SORTED_BAM_DIR/${base}.sorted.bam"
    echo "✅ BAM indexed."
    
    # Convert BAM to sorted BED
    bedtools bamtobed -i "$SORTED_BAM_DIR/${base}.sorted.bam" | bedtools sort -i > "$BED_FILES_DIR/${base}.sorted.bed"
    echo "✅ BED created."
done
  
echo -e "\n=============================================== KARYOTYPE ================================================"
    
# Identifying overlapping regions across samples
echo "Identifying overlapping regions across samples..."
mapfile -t SORTED_BEDS < <(find "$BED_FILES_DIR" -name "*.sorted.bed")

if (( ${#SORTED_BEDS[@]} < 2 )); then
    echo "⚠️ Only one BED file found – skipping multiinter and using it directly."
    cp "${SORTED_BEDS[0]}" "$BED_FILES_DIR/common_intervals.bed"
else
    bedtools multiinter -header -i "${SORTED_BEDS[@]}" > "$BED_FILES_DIR/common_intervals.bed"
fi

echo "Generating karyotype plot..."
Rscript "$ROOT_DIR/scripts/helper_scripts/karyotype.R" "$BED_FILES_DIR/common_intervals.bed" && \
echo "✅ Karyotype plot generated."

awk '$4 > 1' "$BED_FILES_DIR/common_intervals.bed" > "$BED_FILES_DIR/common_intervals_filtered.bed"

echo -e "\n================================================== BLAST =================================================="
    
# Prepare BLAST query file
BLAST_QUERY="$HOST_DNA_ANALYSIS_DIR/blast_query.fasta" 
> "$BLAST_QUERY"

echo -e "\nExtracting sequences for BLAST query..."
while IFS=$'\t' read -r chrom start end _; do
    for FILE in "$SORTED_BAM_DIR"/*.sorted.bam; do
        samtools view "$FILE" "$chrom:$start-$end" | awk '{print ">" $1 "\n" $10}' >> "$BLAST_QUERY"
    done
done < <(tail -n +2 "$BED_FILES_DIR/common_intervals_filtered.bed")
echo "✅ Sequences saved to BLAST query."

# Skip BLAST if BED region file is empty
if [ $(tail -n +2 "$BED_FILES_DIR/common_intervals_filtered.bed" | wc -l) -eq 0 ]; then
    echo "❌ No regions found. Skipping BLAST."
    rm -f "$BLAST_QUERY"
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
        Rscript "$ROOT_DIR/scripts/helper_scripts/human_aligned_tree.R" "$HOST_DNA_ANALYSIS_DIR/combined_blast_results.txt" && \
        echo "✅ Taxonomy tree generated."
    else
        echo "❌ Query file >100KB. Skipping BLAST."
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
