#!/bin/bash

# Script to process and combine GenBank files
# This script takes a Python script for combining GenBank files, two input directories containing GenBank files,
# an output directory, a file extension for GenBank files, and an optional dry-run flag.
# It finds matching files in both input directories based on the given extension and processes them.

# Usage:
# ./script.sh <Combine_Genebank.py> <folder_in1> <folder_in2> <folder_out> <genbank_extension> [--dry-run]
# Example:
# ./script.sh Combine_Genebank.py data/genbank1 data/genbank2 results gb --dry-run

# Check for correct usage
if [ "$#" -lt 5 ]; then
    echo "Usage: $0 <Combine_Genebank.py> <folder_in1> <folder_in2> <folder_out> <genbank_extension> [--dry-run]"
    echo "Example: $0 Combine_Genebank.py data/genbank1 data/genbank2 results gb --dry-run"
    exit 1
fi

# Get arguments
combine_script="$1"
folder_in1="$2"
folder_in2="$3"
folder_out="$4"
gb_extension="$5"
dry_run=false

# Check for optional dry-run flag
if [ "$#" -eq 6 ] && [ "$6" == "--dry-run" ]; then
    dry_run=true
fi

# Check if the Python script exists
if [ ! -f "$combine_script" ]; then
    echo "Error: Python script $combine_script not found."
    exit 1
fi

# Create the output folder if it doesn't exist
mkdir -p "$folder_out"

# Iterate over all GenBank files with the specified extension in folder_in1
for file1 in "$folder_in1"/*.$gb_extension; do
    # Get the base filename
    base_name=$(basename "$file1")
    file2="$folder_in2/$base_name"

    # Check if the corresponding file exists in folder_in2
    if [ -f "$file2" ]; then
        # Define the output file path
        output_file="$folder_out/$base_name"
        
        if [ "$dry_run" = true ]; then
            echo "[DRY RUN] Would process: $file1 and $file2 -> $output_file"
        else
            # Run the Python script
            python "$combine_script" "$file1" "$file2" "$output_file"
            if [ $? -eq 0 ]; then
                echo "Processed: $base_name"
            else
                echo "Failed to process: $base_name"
            fi
        fi
    else
        echo "Warning: Matching file not found for $base_name in $folder_in2"
    fi
done
