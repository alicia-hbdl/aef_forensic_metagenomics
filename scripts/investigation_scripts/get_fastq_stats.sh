#!/bin/bash

# Check if the argument is a valid directory
if [ ! -d "$1" ]; then
  echo "Error: $1 is not a valid directory."
  echo "Usage: $0 <directory>"
  exit 1
fi

# Set the directory to the argument provided by the user
DIRECTORY=$1

# Extract read lengths from each FASTQ file in the provided directory and calculate mean and median
for file in "$DIRECTORY"/*.fastq; do
  echo -n "$file: "  # Print the filename for reference

  # Use awk to process the FASTQ file:
  # NR%4==2 selects the second line in each 4-line FASTQ entry (the sequence line).
  # We calculate the length of the sequence and store it in an array.
  awk 'NR%4==2 { lengths[NR/4] = length($0) } 
       END {
         n = length(lengths)
         sum = 0
         for (i=1; i<=n; i++) {
           sum += lengths[i]
         }
         mean = sum / n
         # Calculate median
         if (n % 2 == 1) {
           median = lengths[int(n/2) + 1]
         } else {
           median = (lengths[n/2] + lengths[n/2 + 1]) / 2
         }
         # Print the results
         printf "mean = %.1f, median = %.1f\n", mean, median
       }' "$file"
done