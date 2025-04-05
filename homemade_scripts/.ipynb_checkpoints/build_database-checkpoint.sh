#!/bin/bash
#SBATCH --job-name=db_build
#SBATCH --output=logs/db_build.log
#SBATCH --error=logs/db_build.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --partition=cpu
 
# Activate Conda environment
if ! command -v conda &> /dev/null; then
    echo "❌ Error: Conda not found. Exiting."
    exit 1
fi

eval "$(conda shell.bash hook)"
conda activate final_project
echo "✅ Conda environment 'final_project' activated."

# Define database and genome storage paths
DBNAME="/scratch/users/k24087895/final_project/metagenome_analysis/data/databases/custom_db"
GENOMES="/scratch/users/k24087895/final_project/metagenome_analysis/data/databases/genomes"
mkdir -p "$DBNAME" "$GENOMES"
cd "$GENOMES"

# Download taxonomy 
echo "Downloading NCBI taxonomy to: $DBNAME"
kraken2-build --download-taxonomy --db "$DBNAME" --threads 8

# Download human library 
echo "Downloading human reference library to: $DBNAME"
kraken2-build --download-library human --db "$DBNAME" --threads 8

# Download reference genomes (if not already present)
echo "Downloading selected genomes from NCBI..."
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/196/035/GCF_000196035.1_ASM19603v1/GCF_000196035.1_ASM19603v1_genomic.fna.gz
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/961/225/GCF_029961225.1_ASM2996122v1/GCF_029961225.1_ASM2996122v1_genomic.fna.gz
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/393/015/GCF_000393015.1_Ente_faec_T5_V1/GCF_000393015.1_Ente_faec_T5_V1_genomic.fna.gz
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/091/045/GCF_000091045.1_ASM9104v1/GCF_000091045.1_ASM9104v1_genomic.fna.gz

# Decompress downloaded genomes
gunzip -f *.fna.gz

# Map file names to correct NCBI Taxonomy IDs
declare -A taxids=(
    [GCF_000196035.1_ASM19603v1_genomic.fna]=169963
    [GCF_000006765.1_ASM676v1_genomic.fna]=208964
    [GCF_000009045.1_ASM904v1_genomic.fna]=224308
    [GCF_000005845.2_ASM584v2_genomic.fna]=562
    [GCF_000008865.2_ASM886v2_genomic.fna]=562
    [GCF_029961225.1_ASM2996122v1_genomic.fna]=1613
    [GCF_000006945.2_ASM694v2_genomic.fna]=28901
    [GCF_000393015.1_Ente_faec_T5_V1_genomic.fna]=1351
    [GCF_000013425.1_ASM1342v1_genomic.fna]=1280
    [GCF_000146045.2_R64_genomic.fna]=4932
    [GCF_000091045.1_ASM9104v1_genomic.fna]=5207
)

# Rewrite headers to Kraken2 format: >kraken:taxid|TAXID|ACCESSION DESCRIPTION
for genome in "$GENOMES"/*.fna; do
    base=$(basename "$genome")
    taxid=${taxids[$base]}
    if ! grep -q "kraken:taxid" "$genome"; then
          echo "Tagging $base with taxid $taxid..."
          sed -i -E "s/^>([A-Za-z0-9_.]+) (.*)/>kraken:taxid|$taxid|\1 \2/" "$genome"
    fi    
    head -n 1 "$genome"
done

# Add genome sequences to the Kraken2 database
for genome in "$GENOMES"/*.fna; do
    echo "Adding $(basename "$genome") to Kraken2 DB..."
    kraken2-build --add-to-library "$genome" --db "$DBNAME" --threads 8
done

# Build the database with a 16G size limit
echo "Building Kraken2 database (max size = 16G)..."
kraken2-build --build --db "$DBNAME" --max-db-size 16G --threads 8 

# Clean up intermediate files
echo "Cleaning up intermediate files..."
kraken2-build --clean --db "$DBNAME"

echo "✅ Kraken2 database build completed successfully!"