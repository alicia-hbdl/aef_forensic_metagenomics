#!/bin/bash
#SBATCH --job-name=build_db
#SBATCH --output=logs/build_db.log
#SBATCH --error=logs/build_db.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=12:00:00
#SBATCH --partition=cpu

# Default values
GENOMES=/scratch/users/k24087895/final_project/data/genomes
DBNAME=/scratch/users/k24087895/final_project/data/databases/custom_notaxid

THREADS=8
KMER_LEN=35


# Check if Conda is installed and activate the environment
if ! command -v conda &> /dev/null; then
    echo "❌ Conda not found. Exiting."
    exit 1
fi
eval "$(conda shell.bash hook)" || { echo "❌ Conda shell integration not set up. Run 'conda init bash' and restart the shell. Exiting."; exit 1; }
conda activate metagenomics || { echo "❌ Failed to activate Conda environment 'metagenomics'. Exiting."; exit 1; }
echo "✅ Conda environment 'final_project' activated."

# Parse command-line arguments
# Ensure required arguments are provided
if [ -z "$DBNAME" ] || [ -z "$GENOMES" ]; then
    echo "❌ Missing required arguments. Using default values:"
    echo "    DBNAME: $DBNAME"
    echo "    GENOMES: $GENOMES"
else
    echo "✅ Arguments provided:"
    echo "    DBNAME: $DBNAME"
    echo "    GENOMES: $GENOMES"
fi

# Create necessary directories and navigate to GENOMES
mkdir -p "$DBNAME" "$GENOMES"
cd "$GENOMES" || { echo "❌ Failed to change directory to '$GENOMES'. Exiting."; exit 1; }

# Download NCBI Taxonomy
echo "📦 Downloading NCBI taxonomy to: $DBNAME"
if ! kraken2-build --download-taxonomy --db "$DBNAME" --threads $THREADS; then
    echo "❌ Error: Failed to download NCBI taxonomy. Exiting."
    exit 1
fi
echo "✅ NCBI taxonomy downloaded to: $DBNAME"

# Optional: Download Human Reference Library
echo "📦 Downloading human reference library to: $DBNAME"
if ! kraken2-build --download-library human --db "$DBNAME" --threads $THREADS; then
    echo "⚠️ Warning: Failed to download human reference library. Skipping."
else
    echo "✅ Human reference library downloaded."
fi

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
echo "📥 Downloading genomes..."
for url in "${genomes_urls[@]}"; do
    wget -nc --tries=3 "$url"  # Retry download if needed
    filename=$(basename "$url")
    gunzip -f "$filename"  # Decompress downloaded genome
done

# Map file names to correct taxIDs
declare -A taxids=(
    [GCF_000196035.1_ASM19603v1_genomic.fna]=169963  # Listeria monocytogenes EGD-e
    [GCF_000006765.1_ASM676v1_genomic.fna]=208964  # Pseudomonas aeruginosa PAO1
    [GCF_000009045.1_ASM904v1_genomic.fna]=224308  # Bacillus subtilis subsp. subtilis str. 168
    [GCF_000005845.2_ASM584v2_genomic.fna]=562      # Escherichia coli
    [GCF_000008865.2_ASM886v2_genomic.fna]=562      # Escherichia coli
    [GCF_029961225.1_ASM2996122v1_genomic.fna]=1613  # Limosilactobacillus fermentum
    [GCF_000006945.2_ASM694v2_genomic.fna]=28901    # Salmonella enterica
    [GCF_000393015.1_Ente_faec_T5_V1_genomic.fna]=1351 # Enterococcus faecalis
    [GCF_000013425.1_ASM1342v1_genomic.fna]=1280     # Staphylococcus aureus
    [GCF_000146045.2_R64_genomic.fna]=4932           # Saccharomyces cerevisiae
    [GCF_000091045.1_ASM9104v1_genomic.fna]=5207     # Cryptococcus neoformans
)

# Tag genome headers with taxIDs
echo "🧬 Tagging genome headers with taxIDs..."
#for genome in "$GENOMES"/*.fna; do
 #   base=$(basename "$genome")
  #  taxid=${taxids["$base"]}
    
   # if [ -z "$taxid" ]; then
    #    echo "⚠️ Warning: No taxid found for $base. Skipping."
     #   continue
  #  fi
    
   # if ! grep -q "kraken:taxid" "$genome"; then
    #    echo "Tagging $base with taxid $taxid..."
     #   sed -i -E "s/^(>[^ ]+)/&|kraken:taxid|$taxid/" "$genome"
   # else
    #    echo "⚠️ Skipping $base: already tagged."
   # fi
#done

# Add genomes to Kraken2 DB
echo "📒 Adding genomes to Kraken2 DB..."
for genome in "$GENOMES"/*.fna; do
    kraken2-build --add-to-library "$genome" --db "$DBNAME" --threads $THREADS &
done
wait
echo "✅ All genomes added to Kraken2 DB."

# Build Kraken2 database
echo "🧬 Building Kraken2 database..."
if ! kraken2-build --build --db "$DBNAME" --threads $THREADS; then
    echo "❌ Error: Failed to build Kraken2 database. Exiting."
    exit 1
fi
echo "✅ Kraken2 database built successfully!"

KRAKEN2_DIR=$(dirname $(which kraken2))  # Get the directory of the Kraken2 executable

# Build Bracken database for different read lengths
for READ_LEN in 50 100 150 200 250 300; do
    echo "🔄 Building Bracken database for read length $READ_LEN..."
    if ! bracken-build -k "$KMER_LEN" -l "$READ_LEN" -d "$DBNAME" -x "$KRAKEN2_DIR" -y kraken2 -t "$THREADS" ; then
        echo "❌ Error: Failed to build Bracken database for read length $READ_LEN. Exiting."
    fi
    echo "✅ Bracken database for read length $READ_LEN built."
done

echo "🎉 All Bracken databases built successfully!"
