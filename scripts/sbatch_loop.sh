#!/bin/bash

databases=(
  /scratch/grp/msc_appbio/recovered/k2_standard_20241228
  /scratch/grp/msc_appbio/recovered/core_nt_database
  /scratch/users/k24087895/final_project/metagenome_analysis/data/databases/k2_standard_08gb_20241228
  /scratch/users/k24087895/final_project/metagenome_analysis/data/databases/k2_pluspfp_16gb_20241228
  /scratch/users/k24087895/final_project/metagenome_analysis/data/databases/k2_pluspfp_08gb_20241228
  /scratch/users/k24087895/final_project/metagenome_analysis/data/databases/k2_standard_16gb_20241228
  /scratch/users/k24087895/final_project/metagenome_analysis/data/databases/k2_eupathdb48_20230407
  /scratch/users/k24087895/final_project/metagenome_analysis/data/databases/custom
)

thresholds=(0 5 10 20 50 100 500)
min_hits=(1 2 3 4 5)
read_lengths=(75 100 150)

for db in "${databases[@]}"; do
  for th in "${thresholds[@]}"; do
    for mh in "${min_hits[@]}"; do
      for len in "${read_lengths[@]}"; do
        jid=$(sbatch --parsable pipeline.sh -fq reads -db "$db" -t "$th" -m "$mh" -l "$len")
        while squeue -j "$jid" >/dev/null 2>&1; do sleep 60; done  # check every 1 minute
      done
    done
  done
done
