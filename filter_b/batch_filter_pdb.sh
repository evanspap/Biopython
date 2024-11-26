#!/bin/bash

# Usage: ./batch_filter_pdb.sh <filter_pdb_script> <input_directory> <output_directory> <bfactor_threshold> <neighbor_residues>

# Check if the correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <filter_pdb_script> <input_directory> <output_directory> <bfactor_threshold> <neighbor_residues>"
    exit 1
fi

# Assign arguments to variables
filter_pdb_script=$1
input_directory=$2
output_directory=$3
bfactor_threshold=$4
neighbor_residues=$5

# Ensure the output directory exists
mkdir -p "$output_directory"

# Process each PDB file in the input directory
for pdb_file in "$input_directory"/*.pdb; do
    if [ -f "$pdb_file" ]; then
        # Get the basename of the file
        filename=$(basename "$pdb_file")
        output_pdb="$output_directory/$filename"

        # Run the filter_pdb.py script
        echo "Processing $pdb_file..."
        python "$filter_pdb_script" "$pdb_file" "$output_pdb" "$bfactor_threshold" "$neighbor_residues"

        # Check if the command succeeded
        if [ $? -eq 0 ]; then
            echo "Filtered PDB saved to $output_pdb"
        else
            echo "Error processing $pdb_file"
        fi
    else
        echo "No PDB files found in $input_directory."
    fi
done
