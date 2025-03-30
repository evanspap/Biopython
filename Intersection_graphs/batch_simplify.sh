#!/bin/bash

# ==============================================================
# Batch CSV Simplifier Script
# Description:
#   Runs the provided CSV_Simplifier.py script on all .csv files
#   in the specified input folder, using the given pocket_base.
#   Results are saved as *_simple.csv in the specified output folder.
#
# Usage:
#   ./batch_simplify.sh <simplifier_script.py> <input_folder> <pocket_base> <output_folder>
#
# Example:
#   ./batch_simplify.sh ./CSV_Simplifier.py ./input Pocket ./output
# ==============================================================

# Check for correct number of arguments
if [ "$#" -ne 4 ]; then
  echo "=============================================================="
  echo "Batch CSV Simplifier Script"
  echo
  echo "Usage:"
  echo "  ./batch_simplify.sh <simplifier_script.py> <input_folder> <pocket_base> <output_folder>"
  echo
  echo "Example:"
  echo "  ./batch_simplify.sh ./CSV_Simplifier.py ./input Pocket ./output"
  echo "=============================================================="
  exit 1
fi

simplifier_script="$1"
input_folder="$2"
pocket_base="$3"
output_folder="$4"

# Check if the simplifier script exists
if [ ! -f "$simplifier_script" ]; then
  echo "Error: simplifier script '$simplifier_script' not found."
  exit 2
fi

# Create output folder if it doesn't exist
mkdir -p "$output_folder"

# Loop over all .csv files in the input folder
for infile in "$input_folder"/*.csv; do
  filename=$(basename "$infile")
  base="${filename%.csv}"
  outfile="${output_folder}/${base}_simple.csv"

  echo "Running: python3 \"$simplifier_script\" \"$infile\" \"$pocket_base\" \"$outfile\""
  python3 "$simplifier_script" "$infile" "$pocket_base" "$outfile"
done
