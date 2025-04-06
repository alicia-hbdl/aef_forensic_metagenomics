#!/bin/bash

databases=(db1 db2 db3)
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
