#!/usr/bin/env bash

# This Bash script downloads housepet genomes from NCBI and tags each with its corresponding NCBI taxonomic ID for Kraken2 compatibility.

# Check for required argument
[ -z "$1" ] && { echo "❌  No genomes directory provided."; echo "Usage: $0 <path/to/genomes>"; exit 1; }
GENOMES="$1"
echo "Using genomes directory: $GENOMES"

# Recreate $GENOMES directory
[ -d "$GENOMES" ] && echo "⚠️  Removing existing $(basename "$GENOMES") directory." && rm -rf "$GENOMES"
mkdir -p "$GENOMES" && echo "✅  Created $(basename "$GENOMES") directory."

cd "$GENOMES" || { echo "❌  Failed to enter '$GENOMES'."; exit 1; }
echo "Current directory: $PWD"

# List of genome URLs to download
genomes_urls=(
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/350/175/GCF_018350175.1_F.catus_Fca126_mat1.0/GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/100/685/GCF_011100685.1_UU_Cfam_GSD_1.0/GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/041/296/265/GCF_041296265.1_TB-T2T/GCF_041296265.1_TB-T2T_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/668/045/GCF_003668045.3_CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/190/915/GCF_034190915.1_mCavPor4.1/GCF_034190915.1_mCavPor4.1_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.3_ARS-UCD2.0/GCF_002263795.3_ARS-UCD2.0_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/964/237/555/GCF_964237555.1_mOryCun1.1/GCF_964237555.1_mOryCun1.1_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/323/735/GCF_036323735.1_GRCr8/GCF_036323735.1_GRCr8_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/025/GCF_000003025.6_Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_genomic.fna.gz"
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
    [GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fna.gz]=9685         # Felis catus (domestic cat)
    [GCF_000001635.27_GRCm39_genomic.fna.gz]=10090                      # Mus musculus (house mouse)
    [GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.fna.gz]=9615               # Canis lupus familiaris (dog)
    [GCF_041296265.1_TB-T2T_genomic.fna.gz]=9796                        # Equus caballus (horse)
    [GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna.gz]=10029             # Cricetulus griseus (Chinese hamster)
    [GCF_034190915.1_mCavPor4.1_genomic.fna.gz]=10141                   # Cavia porcellus (domestic guinea pig)
    [GCF_002263795.3_ARS-UCD2.0_genomic.fna.gz]=9913                    # Bos taurus (domestic cattle)
    [GCF_964237555.1_mOryCun1.1_genomic.fna.gz]=9986                    # Oryctolagus cuniculus (rabbit)
    [GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz]=9031   # Gallus gallus (chicken)
    [GCF_036323735.1_GRCr8_genomic.fna.gz]=9823                         # Sus scrofa (pig)
    [GCF_000003025.6_Sscrofa11.1_genomic.fna.gz]=10116                  # Rattus norvegicus (Norway rat)

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
