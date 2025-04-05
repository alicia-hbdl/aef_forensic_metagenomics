#!/bin/bash

# Set up Conda environment
conda create --name final_project -y
conda activate final_project

# Configure Conda channels
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels r

# Install required tools and packages
conda install -c bioconda kraken2 bracken krakentools samtools bowtie2 bedtools blast sra-tools -y
conda install -c r r-essentials -y
conda install -c conda-forge r-dplyr r-pheatmap r-viridis r-tidyverse -y

# Load required modules (for HPC users)
module load samtools/1.17-gcc-13.2.0-python-3.11.6 \
            bowtie2/2.5.1-gcc-13.2.0-python-3.11.6 \
            py-numpy/1.26.1-gcc-13.2.0-python-3.11.6-blas

# Install KrakenTools
mkdir -p /scratch/users/k24087895/final_project/kraken_experimentation/tools
cd /scratch/users/k24087895/final_project/kraken_experimentation/tools
git clone https://github.com/jenniferlu717/KrakenTools

# Verify Python Biopython installation
python -c "import Bio, numpy; print(f'Biopython version: {Bio.__version__}'); print(f'NumPy version: {numpy.__version__}')"

# Install R packages
Rscript -e '
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes");
remotes::install_github("fbreitwieser/pavian");
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager");
BiocManager::install(c("Rsamtools", "karyoploteR", "GenomicRanges"));
'

echo "Setup complete!"
