#!/bin/bash

# Script to run in the command line to determine the distribution of the FASTQ lengths in the classified and unclassified directories

for file in *.fastq; do
  echo -n "$file: "

  awk 'NR%4==2 {print length($0)}' "$file" | sort -n | awk '
    NR == 1 { min = $1; sum = $1 }
    NR > 1 { sum += $1 }
    { a[NR] = $1; max = $1 }
    END {
      count = NR
      mean = sum / count

      mid = int((count + 1) / 2)
      if (count % 2 == 1) {
        median = a[mid]
      } else {
        median = (a[mid] + a[mid + 1]) / 2
      }

      # Count how many sequences are ≤ median
      below_or_equal = 0
      for (i = 1; i <= count; i++) {
        if (a[i] <= median) {
          below_or_equal++
        } else {
          break
        }
      }
      percent = (below_or_equal / count) * 100

      printf "min = %d, mean = %.1f, median = %.1f, max = %d, ≤ median = %.1f%%\n", min, mean, median, max, percent
    }'
done
