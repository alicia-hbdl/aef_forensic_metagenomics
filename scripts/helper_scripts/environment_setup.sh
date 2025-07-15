#!/bin/bash

# This script sets up a Conda environment based on a provided YAML file.

# Usage:
#   ./helper_scripts/environment_setup.sh <path/to/environment.yml>

# Check if an argument is passed and if it ends with .yml
if [[ -z "$1" || "$1" != *.yml || ! -f "$1" ]]; then
    echo "❌ Error: No valid environment YAML file provided. Usage: $0 <path/to/environment.yml>"
    exit 1
fi

ENV_FILE="$1"
ENV_NAME="$(basename "$ENV_FILE" .yml)" # Get environment name by removing path and .yml extension

# Check for Conda installation and initialize shell integration
if ! command -v conda &> /dev/null; then
    echo "❌ Conda not found. Please install Conda."
    exit 1
fi

eval "$(conda shell.bash hook)" || { echo "❌ Conda shell integration not initialized. Run 'conda init bash'."; exit 1; }

# Create environment if it doesn't exist
if ! conda env list | grep -q "$ENV_NAME"; then
    echo "🔍 '$ENV_NAME' environment not found. Creating it now from '$ENV_FILE'..."
    conda env create -f "$ENV_FILE" || { echo "❌ Failed to create '$ENV_NAME' environment."; exit 1; }
    echo "✅ '$ENV_NAME' environment created."
fi 

# Activate the environment
conda activate "$ENV_NAME" || { echo "❌ Failed to activate '$ENV_NAME' environment."; exit 1; }
echo "✅ '$ENV_NAME' environment activated successfully."

# Ensure that any executable programs installed in the environment are used first when running commands
export PATH="$HOME/.conda/envs/$ENV_NAME/bin:$PATH"
# Ensure that Python will use the modules installed in the Conda environment’s site-packages directory first
export PYTHONPATH="$HOME/.conda/envs/$ENV_NAME/lib/python$(python -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')/site-packages:$PYTHONPATH"
