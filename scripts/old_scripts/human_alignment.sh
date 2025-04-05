#!/bin/bash

# Check if the path to the analysis folder is provided
if [ $# -lt 1 ]; then
  echo -e "❗Usage: $0 <path/to/analysis_folder>\n"
  exit 1
fi

# Define paths to input/output directories
FOLDER="$1"
HUMAN_READS="$FOLDER/data/human_reads"
ALIGNED_SAM="$HUMAN_READS/aligned_sam"
SORTED_BAM="$HUMAN_READS/sorted_bam"
BED_FILES="$HUMAN_READS/bed_files"
ALIGNMENT_SIMILARITY="$FOLDER/outputs/alignment_similarity"
COMMON_REGIONS="$ALIGNMENT_SIMILARITY/common_intersections.bed"
BLAST_QUERY="$HUMAN_READS/blast_query.fasta" 
BLAST_RESULTS="$HUMAN_READS/blast_results"

# Paths to custom R scripts
JACCARD="$FOLDER/../homemade_scripts/jaccard_similarity_heatmap.R"
KARYOTYPE="$FOLDER/../homemade_scripts/karyotype.R"

# Check if SAM files exist before proceeding
if [ ! -d "$ALIGNED_SAM" ] || ! ls "$ALIGNED_SAM"/*.sam >/dev/null 2>&1; then
    echo -e "\n❌ Error: No SAM files found in '$ALIGNED_SAM'."
    exit 1
fi
echo -e "\n✅ Starting analysis of human-aligned genomic sequences in '$ALIGNED_SAM'..."

# Create necessary directories if they don't exist
mkdir -p "$SORTED_BAM" "$BED_FILES" "$ALIGNMENT_SIMILARITY" "$BLAST_RESULTS"

# Process each SAM file
for sam_file in "$ALIGNED_SAM"/*_alignments.sam; do
    sample=$(basename "$sam_file" _alignments.sam)
    echo -e "\nProcessing sample: $sample"

    # Convert SAM to BAM, filter unmapped reads (-F 4), and sort the BAM file
    samtools view -b -h -F 4 "$sam_file" | samtools sort -@8 -o "$SORTED_BAM/${sample}_human.sorted.bam"
    echo " ✅ BAM file created and sorted."
    
    samtools index "$SORTED_BAM/${sample}_human.sorted.bam"
    echo " ✅ BAM file indexed."

    # Convert BAM to BED format and sort the entries
    bedtools bamtobed -i "$SORTED_BAM/${sample}_human.sorted.bam" | sort -k1,1V -k2,2n  > "$BED_FILES/${sample}.sorted.bed"
    echo " ✅ BED file created."
done

echo -e "\nIdentifying overlapping genomic positions across samples..."
mapfile -t sorted_beds < <(find "$BED_FILES" -name "*.sorted.bed")  # Collect all sorted BED file paths
bedtools multiinter -header -i "${sorted_beds[@]}" > "$COMMON_REGIONS.tmp"  # Find common genomic regions

# Generate karyotype plot based on common genomic regions
echo -e "\nGenerating karyotype frequency plot using all identified regions..."
Rscript "$KARYOTYPE" "$COMMON_REGIONS.tmp" 

# Filter out regions where fewer than 3 samples are aligned (retain ≥2 overlapping samples)
awk 'NR==1 || $4 > 2' "$COMMON_REGIONS.tmp" > "$COMMON_REGIONS" && rm "$COMMON_REGIONS.tmp"
echo " ✅ Genomic regions where at least 3 samples align stored in '$COMMON_REGIONS'."

# BLASTing human-aligned sequences
> "$BLAST_QUERY"

# Extract sequences from common genomic regions and store in BLAST query file
echo -e "\n🔹 Extracting sequences for BLAST query..."
while IFS=$'\t' read -r chrom start end _; do  # Read each genomic region (chromosome, start, and end positions)
    for bam_file in "$SORTED_BAM"/*.sorted.bam; do  # Iterate through all sorted BAM files 
        samtools view "$bam_file" "$chrom:$start-$end" | awk '{print ">" $1 "\n" $10}' >> "$BLAST_QUERY"  
    done
done < <(tail -n +2 "$COMMON_REGIONS") # Read from the common regions file, skipping the header row
echo " ✅ Sequences extracted and stored in BLAST query file."

if [ $(tail -n +2 "$COMMON_REGIONS" | wc -l) -eq 0 ]; then
    echo "❌ COMMON_REGIONS file is empty. Skipping sequence extraction and BLAST search."
else
	# Deduplicate (ID, sequence) pairs in BLAST query file.
	awk '
	    NR % 2 == 1 { id = $0; next }  # Store header line (ID)
	    { seq = $0; pair = id "\n" seq }  # Store sequence and create a unique pair
	    !seen[pair]++ { print id; print seq }  # Print only unique (ID, sequence) pairs
	' "$BLAST_QUERY" > "$BLAST_QUERY.tmp" && mv "$BLAST_QUERY.tmp" "$BLAST_QUERY"
	echo " ✅ Duplicate sequences removed from BLAST query."

	# Check BLAST query file size and run BLAST search if valid
	if [ $(stat -c%s "$BLAST_QUERY") -lt $((100 * 1024)) ]; then
	    echo -e "\nRunning BLAST search on extracted sequences..."
	    blastn -query "$BLAST_QUERY" -db nt -out "$BLAST_RESULTS/combined_blast_results.txt" \
	           -evalue 1e-5 -max_target_seqs 5 -outfmt "6 qseqid staxids pident evalue bitscore" -remote
	    echo " ✅ BLAST search completed. Results saved to '$BLAST_RESULTS/combined_blast_results.txt'."
	else
	    echo "❌ BLAST query file exceeds 100KB. Skipping BLAST search."
	    echo "Run selected sequences manually at https://blast.ncbi.nlm.nih.gov/Blast.cgi"
	fi
fi 

# Generate Jaccard similarity heatmap for alignment similarity
echo -e "\nComputing Jaccard similarity across samples..."
JACCARD_RESULTS="$ALIGNMENT_SIMILARITY/jaccard_results.txt"
echo -e "Sample1\tSample2\tIntersection\tUnion\tJaccard\tN_Intersections" > "$JACCARD_RESULTS"

# Compute Jaccard similarity between all pairs of sorted BED files
for file1 in "${sorted_beds[@]}"; do
    for file2 in "${sorted_beds[@]}"; do
        if [[ "$file1" != "$file2" ]]; then
            sample1=$(basename "$file1" .sorted.bed)
            sample2=$(basename "$file2" .sorted.bed)
            jaccard_output=$(bedtools jaccard -a "$file1" -b "$file2" | tail -n1)
            echo -e "$sample1\t$sample2\t$jaccard_output" >> "$JACCARD_RESULTS"
        fi
    done
done

# Generate heatmap visualization of Jaccard similarity
echo -e "\n Creating Jaccard similarity heatmap..."
Rscript "$JACCARD" "$JACCARD_RESULTS"
	
