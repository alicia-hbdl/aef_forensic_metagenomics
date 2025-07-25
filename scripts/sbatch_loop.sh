#!/bin/bash

DB_PATH="/scratch/prj/aef_forensic_metagenomics/data/databases"

databases=(
#"$DB_PATH/k2_core_nt_20240904"
#"$DB_PATH/k2_housepets_250510"
#"$DB_PATH/k2_pluspfp_16gb_20250402"
#"$DB_PATH/k2_standard_08gb_20241228"
#"$DB_PATH/k2_standard_20250402"
#"$DB_PATH/k2_eupathdb48_20230407"
#"$DB_PATH/k2_pluspfp_08gb_20241228"
#"$DB_PATH/k2_pluspfp_20250402"
#"$DB_PATH/k2_standard_16gb_20250402"
#"$DB_PATH/k2_zymobiomics_250509"
)

K2_MIN_HIT_VALUES=(5)
#1 2 3 4)
B_THRESHOLD_VALUES=(0 5 10 20 40)

# Run a single pipeline job with trimming and host DNA removal enabled before looping over all databases
jid_first=$(sbatch --parsable pipeline.sh \
  -f /scratch/prj/aef_forensic_metagenomics/zymobiomics_folder/raw_data \
  -d "$DB_PATH/k2_standard_16gb_20250402" \
  -t -r)

# Wait for the first job to finish
#while squeue -j "$jid_first" >/dev/null 2>&1; do
 #   sleep 60
#done

# Loop through all database paths
for db in "${databases[@]}"; do # Loop through all database paths
  for k2_hit in "${K2_MIN_HIT_VALUES[@]}"; do # Loop through all k2-min-hit values
    for b_thresh in "${B_THRESHOLD_VALUES[@]}"; do # Loop through all b-threshold values
      
      echo "Submitting job for DB: $db | k2-min-hit: $k2_hit | b-threshold: $b_thresh"
      
      jid=$(sbatch --parsable pipeline.sh \
        -f /scratch/prj/aef_forensic_metagenomics/zymobiomics_folder/raw_data \
        -d "$db" --k2-min-hit $k2_hit --b-threshold $b_thresh \
        -g /scratch/prj/aef_forensic_metagenomics/zymobiomics_folder/raw_data/ground_truth.csv)
      
      # Wait for job completion before launching the next
      while squeue -j "$jid" >/dev/null 2>&1; do
        sleep 60
      done

    done
  done
done
