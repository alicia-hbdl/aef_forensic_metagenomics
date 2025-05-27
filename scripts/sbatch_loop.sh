#!/bin/bash

DB_PATH="/scratch/users/k24087895/final_project/data/databases"

databases=(
  "$DB_PATH/k2_eupathdb48_20230407"
  "$DB_PATH/k2_housepets_250510"
  "$DB_PATH/k2_pluspfp_08gb_20241228"
  "$DB_PATH/k2_pluspfp_16gb_20250402"
  "$DB_PATH/k2_standard_08gb_20241228"
  "$DB_PATH/k2_standard_16gb_20250402"
  "$DB_PATH/k2_zymobiomics_250509"
)

# Run a single pipeline job with trimming and host DNA removal enabled before looping over all databases
jid_first=$(sbatch --parsable pipeline.sh \
  -f /scratch/users/k24087895/final_project/zymobiomics_folder/raw_data \
  -d "$DB_PATH/k2_standard_16gb_20250402" \
  -t -r -g /scratch/users/k24087895/final_project/zymobiomics_folder/raw_data/ground_truth.csv)

# Wait for the first job to finish
while squeue -j "$jid_first" >/dev/null 2>&1; do
    sleep 60
done

# Loop over all databases
for db in "${databases[@]}"; do
    jid=$(sbatch --parsable pipeline.sh \
      -f /scratch/users/k24087895/final_project/zymobiomics_folder/raw_data \
      -d "$db" -g /scratch/users/k24087895/final_project/zymobiomics_folder/raw_data/ground_truth.csv)
    
    # Wait for the submitted job to complete before submitting the next
    while squeue -j "$jid" >/dev/null 2>&1; do
        sleep 60
    done
done