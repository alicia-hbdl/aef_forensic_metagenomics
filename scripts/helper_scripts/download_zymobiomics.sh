#!/usr/bin/env bash

# This Bash script downloads ZymoBIOMICS reference genomes from NCBI and tags each with its corresponding NCBI taxonomic ID for Kraken2 compatibility.

# Check for required argument
[ -z "$1" ] && { echo "❌  No genomes directory provided."; echo "Usage: $0 <path/to/genomes>"; exit 1; }
GENOMES="$1"
echo "Using genomes directory: $GENOMES"

mkdir -p "$GENOMES" && echo "✅  Created $(basename "$GENOMES") directory."
cd "$GENOMES" || { echo "❌  Failed to enter '$GENOMES'."; exit 1; }
echo "Current directory: $PWD"

# List of genome URLs to download
genomes_urls=(
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/196/035/GCF_000196035.1_ASM19603v1/GCF_000196035.1_ASM19603v1_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/961/225/GCF_029961225.1_ASM2996122v1/GCF_029961225.1_ASM2996122v1_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/393/015/GCF_000393015.1_Ente_faec_T5_V1/GCF_000393015.1_Ente_faec_T5_V1_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/091/045/GCF_000091045.1_ASM9104v1/GCF_000091045.1_ASM9104v1_genomic.fna.gz"
)

# Download and decompress
echo "Downloading genomes..."
for url in "${genomes_urls[@]}"; do
    wget -nc --tries=3 "$url" || { echo "⚠️  Failed download for $url. Skipping."; continue; }
    gunzip -f "$(basename "$url")"
done
echo "✅  Genomes downloaded and decompressed."


# Map filenames to NCBI taxIDs
declare -A taxids=(
    [GCF_000196035.1_ASM19603v1_genomic.fna]=169963     # Listeria monocytogenes EGD-e
    [GCF_000006765.1_ASM676v1_genomic.fna]=208964       # Pseudomonas aeruginosa PAO1
    [GCF_000009045.1_ASM904v1_genomic.fna]=224308       # Bacillus subtilis subsp. subtilis str. 168
    [GCF_000005845.2_ASM584v2_genomic.fna]=562          # Escherichia coli
    [GCF_000008865.2_ASM886v2_genomic.fna]=562          # Escherichia coli
    [GCF_029961225.1_ASM2996122v1_genomic.fna]=1613     # Limosilactobacillus fermentum
    [GCF_000006945.2_ASM694v2_genomic.fna]=28901        # Salmonella enterica
    [GCF_000393015.1_Ente_faec_T5_V1_genomic.fna]=1351  # Enterococcus faecalis
    [GCF_000013425.1_ASM1342v1_genomic.fna]=1280        # Staphylococcus aureus
    [GCF_000146045.2_R64_genomic.fna]=4932              # Saccharomyces cerevisiae
    [GCF_000091045.1_ASM9104v1_genomic.fna]=5207        # Cryptococcus neoformans
)

# Add taxID tags to genome headers
echo "Tagging genome headers..."
for genome in *.fna; do
    taxid=${taxids["$genome"]}
    if [ -z "$taxid" ]; then
        echo "⚠️  No taxid for $genome. Skipping."
        continue
    fi
    if grep -q "kraken:taxid" "$genome"; then
        echo "⚠️  $genome already tagged. Skipping."
    else
        echo "Tagging $genome with taxid $taxid..."
        sed -i -E "s/^(>[^ ]+)/&|kraken:taxid|$taxid/" "$genome"
    fi
done
