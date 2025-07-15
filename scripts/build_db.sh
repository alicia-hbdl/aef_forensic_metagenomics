#!/usr/bin/env bash

# This script builds a Kraken2 database using genome files from a specified directory.
#SBATCH --job-name=build_db
#SBATCH --output=logs/build_db_%j.log
#SBATCH --error=logs/build_db_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=12:00:00
#SBATCH --partition=cpu

# Usage: ./build_db.sh [-g/--genomes <genomes_dir>] [-d/--database <database_dir>]

set -e  # Exit on error 
set -x  # Print each command and its arguments as it is executed for debugging 

# This script builds a Kraken2 database.
# Arguments:
#   -g, --genomes    Path to the directory containing genome files (default: /scratch/users/k24087895/final_project/data/genomes/zymobiomics)
#   -d, --database   Path to the Kraken2 database directory (default: /scratch/users/k24087895/final_project/data/databases/k2_zymobiomics_[date])

# Default values
THREADS=8
KMER_LEN=35
GENOMES="/scratch/users/k24087895/final_project/data/genomes/zymobiomics"
DBNAME="/scratch/users/k24087895/final_project/data/databases/k2_zymobiomics_$(date +"%Y%m%d")"

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

# Download genomes if missing
if [ ! -d "$GENOMES" ] || ! ls "$GENOMES"/*.fna &> /dev/null; then
    SCRIPT="./helper_scripts/download_$(basename "$GENOMES").sh"
    [ -f "$SCRIPT" ] || { echo "❌  $SCRIPT not found."; exit 1; }
    echo "Downloading genomes..."
    bash "$SCRIPT" "$GENOMES" && echo "✅ Genomes downloaded."
else
    echo "✅ Genomes already present in $GENOMES. Skipping download."
fi

# --- BUILD KRAKEN DATABASE ---

echo "Downloading NCBI taxonomy..."
# Use --use-ftp to avoid rsync (default)
kraken2-build --download-taxonomy --use-ftp --db "$DBNAME" --threads "$THREADS" 2>&1 \
|| { echo "❌  Failed to download NCBI taxonomy. Exiting."; exit 1; }
echo "✅  NCBI taxonomy downloaded."
# head "$DBNAME"/taxonomy/* || { echo "❌  Failed to read files in '$DBNAME/taxonomy/'."; exit 1; } # Troubleshooting

# Add human library
# echo "Installing human reference library..."
# kraken2-build --download-library human --use-ftp --db "$DBNAME" --threads "$THREADS" 2>&1 \
# || echo "⚠️  Failed to install human reference library."
# echo "✅  Human reference library installed."
# head "$DBNAME"/library/human/* || { echo "❌  Failed to read files in '$DBNAME/library/human/'."; exit 1; } # Troubleshooting

# Add custom genomes sequentially
echo "Adding custom genomes..."
for file in "$GENOMES"/*.fna; do
    head -n 1 "$file"
    kraken2-build --add-to-library "$file" --db "$DBNAME" --threads "$THREADS" 2>&1 \
    || echo "⚠️  '$file' failed to add."
done
echo "✅  Custom genomes added."
# head "$DBNAME"/library/added/* || { echo "❌  Failed to read files in '$DBNAME/library/added/'."; exit 1; } # Troubleshooting

echo "Building Kraken2 database..."
export OMP_NUM_THREADS="$THREADS"  # Set OpenMP threads (used only during --build)
kraken2-build --build --db "$DBNAME" --max-db-size 16000000000 --threads "$THREADS" 2>&1 \
|| { echo "❌  Failed to build Kraken2 database. Exiting."; exit 1; } # --max-db-size 16Gb (in bytes)
echo "✅  Kraken2 database built."
# head "$DBNAME"/*.k2d || { echo "❌  Failed to read files in '$DBNAME/'."; exit 1; } # Troubleshooting

# --- BUILD BRACKEN DATABASE ---

# Build Bracken databases for different read lengths
for READ_LEN in 50 75 100 150 200 250 300 ; do 
    bracken-build -d "$DBNAME" -t "$THREADS" -k "$KMER_LEN" -l "$READ_LEN" -x "$(dirname "$(which kraken2)")" 2>&1 \
    || { echo "⚠️  Failed to build Bracken database for read length $READ_LEN. Skipping."; continue; }
    echo "✅  Bracken database built for read length $READ_LEN"
done

echo "All Bracken databases built successfully!"

# --- CLEAN ---
# Remove intermediate files to save space (including taxonomy folder; needed again for future builds)
kraken2-build --clean --db "$DBNAME" 
