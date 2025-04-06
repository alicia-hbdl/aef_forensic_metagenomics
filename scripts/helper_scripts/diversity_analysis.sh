#!/bin/bash

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    -b|--bracken-dir)
      if [[ -z "$2" || ! -d "$2" ]] || ! ls "$2"/*.bracken &>/dev/null; then
        echo "❌ Error: Invalid or missing Bracken files in $2. Usage: $0 -b/--bracken-dir <path> -d/--diversity-dir <path>"
        exit 1
      fi
      BRACKEN_DIR="$2"
      ROOT_DIR="$(dirname "$(dirname "$(dirname "$BRACKEN_DIR")")")"
      shift 2
      ;;
    -d|--diversity-dir)
      if [[ -z "$2" || ! -d "$2" ]]; then
        echo "❌ Error: Diversity directory $2 does not exist. Usage: $0 -b/--bracken-dir <path> -d/--diversity-dir <path>"
        exit 1
      fi
      DIVERSITY_DIR="$2"
      mkdir -p "$DIVERSITY_DIR"
      shift 2
      ;;
    *)
      echo "❌ Unknown option: $1. Usage: $0 -b/--bracken-dir <path> -d/--diversity-dir <path>"
      exit 1
      ;;
  esac
done


# Define alpha diversity metrics
METRICS=("BP" "Sh" "F" "Si" "ISi") # Shannon's, Berger-Parker's, Simpson's, Inverse Simpson's, Fisher's

# Create output file with header
echo -e "Sample\t${METRICS[*]}" | tr ' ' '\t' > "$DIVERSITY_DIR/alpha_diversity.tsv"

echo "Calculating alpha and beta diversity..."

INPUT_FILES=()

# Loop over Bracken files
for file in "$BRACKEN_DIR"/*.bracken; do  
  base=$(basename "$file" ".bracken")
  DIVERSITY_RESULTS=("$base")

  # Compute each alpha diversity metric
  for METRIC in "${METRICS[@]}"; do
    VALUE=$(python "$ROOT_DIR/tools/KrakenTools/DiversityTools/alpha_diversity.py" -f "$file" -a "$METRIC" | awk -F': ' '{if (NF>1) print $2}')
    DIVERSITY_RESULTS+=("$VALUE")
  done

  # Save alpha diversity results
  echo -e "${DIVERSITY_RESULTS[*]}" | tr ' ' '\t' >> "$DIVERSITY_DIR/alpha_diversity.tsv"
  INPUT_FILES+=("$file")
done

echo "✅ Alpha diversity calculated."

# Beta diversity
if (( ${#INPUT_FILES[@]} < 2 )); then
  echo "⚠️ Skipping beta diversity — fewer than 2 samples."
else
  python "$ROOT_DIR/tools/KrakenTools/DiversityTools/beta_diversity.py" -i "${INPUT_FILES[@]}" --type bracken > "$DIVERSITY_DIR/beta_diversity_matrix.tsv"
  
  Rscript "$ROOT_DIR/scripts/helper_scripts/b_diversity_heatmap.R" "$DIVERSITY_DIR/beta_diversity_matrix.tsv"
  
  echo "✅ Beta diversity heatmap generated."
fi