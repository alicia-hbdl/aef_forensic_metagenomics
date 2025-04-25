#!/bin/bash

# Check if a directory argument is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <directory>"
  exit 1
fi

# Set the directory to the argument provided by the user
DIRECTORY=$1

# Extract read lengths from each FASTQ file in the provided directory and calculate mean and median
for file in "$DIRECTORY"/*classified*.fastq; do
  echo -n "$file: "  # Print the filename for reference

  # Process the FASTQ file
  # NR%4==2 selects the second line in each 4-line FASTQ entry (the sequence line)
  # Calculate the length of the sequence and store it
  awk 'NR%4==2 {print length($0)}' "$file" | sort -n | awk '
  
  # Initialize variables for the first line
    NR == 1 { 
      sum = $1;          # Initialize sum with the first read length
      count = 1;         # Initialize read count to 1
      a[1] = $1          # Store the first read length in array a[]
    } 

    # For subsequent lines, accumulate the sum, count, and store the lengths
    NR > 1 {   
      sum += $1;         # Add the current read length to sum
      count++;           # Increment the count of reads
      a[count] = $1      # Store the current read length in array a[]
    }

    # Calculate the mean and median of read lengths
    END {
      mean = sum / count  # Mean is the sum of lengths divided by the number of reads

      mid = int(count / 2) + 1  # Calculate the middle index
      if (count % 2 == 1) {
        # If the number of reads is odd, the median is the middle value
        median = a[mid]
      } else {
        # If the number of reads is even, the median is the average of the two middle values
        median = (a[mid - 1] + a[mid]) / 2
      }

      # Output the mean and median values
      printf "mean = %.1f, median = %.1f\n", mean, median
    }'
done