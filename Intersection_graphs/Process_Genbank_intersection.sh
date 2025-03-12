#!/bin/bash

# ==============================================================
# Script: Process GenBank Files
# Description: This script reads all .gb files from the input folder
#              and processes them using a specified Python script.
#              The output files are saved in the specified output folder
#              with a *_relative.csv suffix.
# Usage: ./Process_Genbank_intersection.sh <python_script> <input_folder> <output_folder>
# Example: ./Process_Genbank_intersection.sh GB_Intersection_Percentage.py genbank_files/ output_results/
# ==============================================================

# Check if correct number of arguments is provided
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <python_script> <input_folder> <output_folder>"
    exit 1
fi

PYTHON_SCRIPT="$1"
INPUT_FOLDER="$2"
OUTPUT_FOLDER="$3"

# Check if the Python script exists
if [[ ! -f "$PYTHON_SCRIPT" ]]; then
    echo "Error: $PYTHON_SCRIPT not found! Please ensure it is in the correct location."
    exit 1
fi

# Ensure the output folder exists
mkdir -p "$OUTPUT_FOLDER"

# Process each .gb file in the input folder
for gb_file in "$INPUT_FOLDER"/*.gb; do
    if [[ -f "$gb_file" ]]; then
        output_file="$OUTPUT_FOLDER/$(basename "${gb_file%.gb}_relative.csv")"
        echo "Processing $gb_file -> $output_file"
        python "$PYTHON_SCRIPT" "$gb_file" "$output_file"
    fi
done

echo "Processing completed!"
