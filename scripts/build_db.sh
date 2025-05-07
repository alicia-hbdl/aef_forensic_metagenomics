#!/bin/bash
#SBATCH --job-name=build_db
#SBATCH --output=logs/build_db.log
#SBATCH --error=logs/build_db.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=12:00:00
#SBATCH --partition=cpu

# Usage: ./build_db.sh [-g/--genomes <path/to/genomes>] [-d/--database <path/to/database>]

set -e # Exit on error 
set -x  # Print each command and its arguments as it is executed for debugging 

# This script builds a Kraken2 database.
# Arguments:
#   -g, --genomes    Path to the directory containing genome files (default: /scratch/users/k24087895/final_project/data/genomes)
#   -d, --database   Path to the Kraken2 database directory (default: /scratch/users/k24087895/final_project/data/databases/k2_custom_[date])

# Default values
THREADS=8
KMER_LEN=35
GENOMES="/scratch/users/k24087895/final_project/data/genomes"
DBNAME="/scratch/users/k24087895/final_project/data/databases/k2_custom_$(date +"%Y%m%d")"

# --- ACTIVATE CONDA ENVIRONMENT --- 
source "$(pwd)/helper_scripts/environment_setup.sh" "$(pwd)/metagenomics.yml" \
|| { echo "❌  Failed to set up Conda environment."; exit 1; }

# --- PARSE COMMAND LINE ARGUMENTS ---

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -g|--genomes) GENOMES="$2"; shift 2 ;; # Set GENOMES and shift
        -d|--database) DBNAME="$2"; shift 2 ;; # Set DBNAME and shift
        *) echo "❌  Unknown option: $1"; echo "Usage: ./build_db.sh [-g <genomes>] [-d <database>]"; exit 1 ;;
    esac
done

# Use defaults if arguments are missing
if [ -z "$GENOMES" ] || [ -z "$DBNAME" ]; then
    echo "❌  Missing arguments. Using default values."
    echo "Usage: ./build_db.sh -g <genomes> -d <database>"
else
    echo "✅  Using provided arguments:"
fi
echo "    DBNAME: $DBNAME"
echo "    GENOMES: $GENOMES"

[ -d "$DBNAME" ]  && echo "⚠️  Removing existing $(basename "$DBNAME") directory." && rm -rf "$DBNAME"
mkdir -p  "$DBNAME" && echo "✅  Created $(basename "$DBNAME") directory."

# --- DOWNLOAD AND TAG GENOMES ---
# Optional: comment this block out if $GENOMES already contains the required .fna files.

# Recreate $GENOMES directory
[ -d "$GENOMES" ] && echo "⚠️  Removing existing $(basename "$GENOMES") directory." && rm -rf "$GENOMES"
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
    wget -nc --tries=3 "$url"
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

# --- BUILD KRAKEN DATABASE ---

echo "Downloading NCBI taxonomy..."
# Use --use-ftp to avoid rsync (default)
kraken2-build --download-taxonomy --use-ftp --db "$DBNAME" --threads "$THREADS" \
|| { echo "❌  Failed to download NCBI taxonomy. Exiting."; exit 1; }
echo "✅  NCBI taxonomy downloaded."

# Troubleshooting: check the contents of the taxonomy files
head "$DBNAME"/taxonomy/* || { echo "❌  Failed to read files in '$DBNAME/taxonomy/'."; exit 1; }

# Add human library
echo "Installing human reference library..."
kraken2-build --download-library human --use-ftp --db "$DBNAME" --threads "$THREADS" \
|| echo "⚠️  Failed to install human reference library."
echo "✅  Human reference library installed."

# Troubleshooting: check contents of human library files
head "$DBNAME"/library/human/* || { echo "❌  Failed to read files in '$DBNAME/library/human/'."; exit 1; }

# Add custom genomes sequentially
echo "Adding custom genomes..."
for file in "$GENOMES"/*.fna; do
    echo "Filepath: $file"; head -n 1 "$file"
    kraken2-build --add-to-library "$file" --db "$DBNAME" --threads "$THREADS" \
    || echo "⚠️  '$file' failed to add."
done
echo "✅  Custom genomes added."

# Troubleshooting: check contents of custom library files
head "$DBNAME"/library/added/* || { echo "❌  Failed to read files in '$DBNAME/library/added/'."; exit 1; }

echo "Building Kraken2 database..."
export OMP_NUM_THREADS="$THREADS"  # Set the number of threads for OpenMP - Only the --build step uses OpenMP.
kraken2-build --build --db "$DBNAME" --threads "$THREADS" \
|| { echo "❌  Failed to build Kraken2 database. Exiting."; exit 1; }
echo "✅  Kraken2 database built."

# Troubleshooting: check contents of database files
head "$DBNAME"/*.k2d || { echo "❌  Failed to read files in '$DBNAME/'."; exit 1; }

# --- BUILD BRACKEN DATABASE ---

# Build Bracken databases for different read lengths
for READ_LEN in 50 100 150 200 250 300 ; do 
    bracken-build -d "$DBNAME" -t "$THREADS" -k "$KMER_LEN" -l "$READ_LEN" -x "$(dirname "$(which kraken2)")" \
    || { echo "⚠️  Failed to build Bracken database for read length $READ_LEN. Skipping."; continue; }
    echo "✅  Bracken database built for read length $READ_LEN"
done

echo "All Bracken databases built successfully!"
