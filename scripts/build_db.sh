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
source "$(pwd)/helper_scripts/environment_setup.sh" "$(pwd)/metagenomics.yml" || { echo "❌ Failed to set up Conda environment."; exit 1; }

# --- PARSE COMMAND LINE ARGUMENTS ---

# Check if any arguments are passed
while [[ $# -gt 0 ]]; do
    case "$1" in
        -g|--genomes)
            GENOMES="$2"  # Assign the next argument to GENOMES
            shift 2  # Move past the argument and value
            ;;
        -d|--database)
            DBNAME="$2"  # Assign the next argument to DBNAME
            shift 2  # Move past the argument and value
            ;;
        *)
            echo "❌ Unknown option: $1"
            echo "Usage: ./build_db.sh [-g/--genomes <path/to/genomes>] [-d/--database <path/to/database>]"
            exit 1
            ;;
    esac
done

# Otherwise use default values
if [ -z "$GENOMES" ] || [ -z "$DBNAME" ]; then
    echo "❌ Missing required arguments."
    echo "Usage: ./build_db.sh [-g/--genomes <path/to/genomes>] [-d/--database <path/to/database>]"
    echo "Using default values:"
else
    echo "✅ Arguments provided:"
fi
echo "    DBNAME: $DBNAME"
echo "    GENOMES: $GENOMES"

# Remove existing directories if they exist 
[ -d "$GENOMES" ] && echo "⚠️ Removing existing Genomes directory." && rm -rf "$GENOMES"
[ -d "$DBNAME" ] && echo "⚠️ Removing existing Database directory." && rm -rf "$DBNAME"

# Create new directories 
mkdir -p "$DBNAME" "$GENOMES" && cd "$GENOMES" || { echo "❌ Failed to change directory to '$GENOMES'."; exit 1; }

# --- DOWNLOAD AND TAG GENOMES ---
# This section can be commented out if the $GENOMES directory already contains the genomes in the correct format.

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

# Download and decompress genomes
echo "Downloading genomes..."
for url in "${genomes_urls[@]}"; do
    wget -nc --tries=3 "$url"       # Retry download if needed
    filename=$(basename "$url")
    gunzip -f "$filename"           # Decompress the downloaded genome file
done

# Tag genome headers with appropriate taxIDs from the provided list
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

# Tag genome headers with taxIDs
echo "Tagging genome headers with taxIDs..."
for genome in "$GENOMES"/*.fna; do
    base=$(basename "$genome")  
    taxid=${taxids["$base"]}    
    if [ -z "$taxid" ]; then # Skip if no taxid
        echo "⚠️ Warning: No taxid found for $base. Skipping."
        continue
    fi
    if ! grep -q "kraken:taxid" "$genome"; then # Skip if already tagged
        echo "Tagging $base with taxid $taxid..."
        sed -i -E "s/^(>[^ ]+)/&|kraken:taxid|$taxid/" "$genome"  # Tagging with taxid
    else
        echo "⚠️ Skipping $base: already tagged."
    fi
done

# --- BUILD KRAKEN DATABASE ---

# Download NCBI taxonomyecho "📦 Downloading NCBI taxonomy to: $DBNAME"
if kraken2-build --download-taxonomy --db "$DBNAME" --threads $THREADS; then
    echo "✅ NCBI taxonomy downloaded to: $DBNAME"
else
    echo "❌ Error: Failed to download NCBI taxonomy. Exiting."
    exit 1
fi

# Download human reference library
echo "Downloading human reference library to: $DBNAME"
if kraken2-build --download-library human --db "$DBNAME" --threads $THREADS; then
    echo "✅ Human reference library downloaded."
else
    echo "❌ Error: Failed to download human reference library.  Exiting."
    exit 1
fi

# Add genomes to the Kraken2 database
echo "Adding genomes to Kraken2 DB..."
for genome in "$GENOMES"/*.fna; do
    ADD_TO_LIB="kraken2-build --add-to-library \"$genome\" --db \"$DBNAME\" --threads $THREADS" 
    echo "$ADD_TO_LIB"
    eval "$ADD_TO_LIB" & # Run in background
done
wait  # Wait for all background jobs to finish
echo "✅ All genomes added to Kraken2 DB."

# Build Kraken2 database
echo "Building Kraken2 database..."
if ! kraken2-build --build --db "$DBNAME" --threads $THREADS; then
    echo "❌ Failed to build Kraken2 database."
    exit 1
fi
echo "✅ Kraken2 database built."

tree "$DBNAME"  # Display the directory structure of the database

# --- BUILD BRACKEN DATABASE ---

KRAKEN2_DIR=$(dirname $(which kraken2))  # Get the directory of the Kraken2 executable

# Build Bracken databases for different read lengths
for READ_LEN in 50 100 150 200 250 300; do
    echo "Building Bracken DB for read length $READ_LEN..."
    echo " $KMER_LEN $READ_LEN $DBNAME $KRAKEN2_DIR $THREADS"
    # Try to build Bracken DB using bracken-build
    if ! bracken-build -k "$KMER_LEN" -l "$READ_LEN" -d "$DBNAME" -y kraken2 -t "$THREADS"; then

        # Rebuild Kraken DB and use kmer2read_distr for Bracken database creation
        echo "❌ Bracken DB build failed. Rebuilding Kraken DB for $READ_LEN..."
        
        if ! kraken2 --db="${DBNAME}" --threads=10 <(find -L "${DBNAME}/library" \( -name "*.fna" -o -name "*.fasta" -o -name "*.fa" \) -exec cat {} + ) > database.kraken; then
            echo "❌ Error rebuilding Kraken DB."
            continue  # Skip to the next iteration of the loop
        fi

        # Run kmer2read_distr with the newly rebuilt Kraken DB
        if ! ./kmer2read_distr --seqid2taxid "${DBNAME}/seqid2taxid.map" --taxonomy "${DBNAME}/taxonomy" --kraken database.kraken --output "database${READ_LEN}mers.kraken" -k "${KMER_LEN}" -l "${READ_LEN}" -t "${THREADS}"; then
            echo "❌ Error running kmer2read_distr."
            continue  
        fi

        # Generate kmer distribution using Python script
        if ! python generate_kmer_distribution.py -i "database${READ_LEN}mers.kraken" -o "database${READ_LEN}mers.kmer_distrib"; then
            echo "❌ Error generating kmer distribution."
            continue  
        fi        
        continue
    fi
    echo "✅ Bracken DB for $READ_LEN built."
done

echo "All Bracken databases built successfully!"
