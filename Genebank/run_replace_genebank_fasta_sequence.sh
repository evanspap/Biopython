#!/bin/bash

# Ensure script exits on error
set -e

# Check for correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <script_path> <folder1> <folder2>"
    exit 1
fi

# Assign arguments to variables
script_path=$1
folder1=$2
folder2=$3

# Iterate through all subdirectories in folder1 looking for *_pocket.gbk files
find "$folder1" -type f -name '*_pocket.gbk' | while read -r gbk_file; do
    # Extract basename without directory and suffix
    basename=$(basename "$gbk_file" "_pocket.gbk")

    # Define corresponding FASTA file in folder2
    fasta_file="$folder2/${basename}.fasta"

    # Define output file path in the same directory as the GeneBank file
    output_file="$(dirname "$gbk_file")/${basename}_pocket_corr.gbk"

    # Check if the FASTA file exists
    if [ ! -f "$fasta_file" ]; then
        echo "FASTA file not found: $fasta_file"
        continue
    fi

    # Run the Python script
    echo "Processing: $gbk_file with $fasta_file -> $output_file"
    python "$script_path" "$gbk_file" "$fasta_file" "$output_file"
done
