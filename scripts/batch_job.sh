#!/bin/bash
#SBATCH --job-name=forensics_pipeline
#SBATCH --output=logs/forensics_pipeline.log  
#SBATCH --error=logs/forensics_pipeline.err
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --partition=cpu

# Stop script on any error
set -e

# Ensure Conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: Conda not found. Exiting."
    exit 1
fi

# Initialize Conda (more robust than `source`)
eval "$(conda shell.bash hook)"

# Activate the Conda environment
conda activate final_project

# Force using Conda's Python instead of system Python
export PATH="/users/k24087895/.conda/envs/final_project/bin:$PATH"
export PYTHONPATH="/users/k24087895/.conda/envs/final_project/lib/python3.13/site-packages:$PYTHONPATH"

# Run script
# First perform quality trimming
# bash /scratch/users/k24087895/final_project/metagenome_analysis/homemade_scripts/quality_trimming.sh /scratch/users/k24087895/final_project/metagenome_analysis/zymobiomics_folder/data/raw_fastq
# Then run the pipeline including human DNA 
bash /scratch/users/k24087895/final_project/metagenome_analysis/homemade_scripts/refinement_pipeline.sh zymobiomics -r
