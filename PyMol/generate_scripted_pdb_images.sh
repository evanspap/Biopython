#!/bin/bash

# Script: generate_scripted_pdb_images.sh
# Description: Generates images for all PDB files in the current directory using PyMOL.
# Usage:
#   generate_scripted_pdb_images.sh <pymol_command_file> <output_directory>

# Ensure PyMOL is installed
if ! command -v pymol &> /dev/null; then
    echo "Error: PyMOL is not installed or not in PATH."
    exit 1
fi

# Ensure the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <pymol_command_file> <output_directory>"
    echo "Example: $0 pymol_commands.txt ./images"
    exit 1
fi

PYMOL_COMMAND_FILE="$1"
OUTPUT_DIR="$2"

# Ensure the command file exists
if [ ! -f "$PYMOL_COMMAND_FILE" ]; then
    echo "Error: PyMOL command file '$PYMOL_COMMAND_FILE' not found."
    exit 1
fi

# Create and ensure the output directory is writable
mkdir -p "$OUTPUT_DIR"
chmod u+w "$OUTPUT_DIR"

# Process all .pdb files in the current directory
for pdb_file in ./*.pdb; do
    # Ensure we have at least one PDB file
    [ -e "$pdb_file" ] || continue

    # Get the absolute path to avoid PyMOL path issues
    pdb_file="$(realpath "$pdb_file")"

    # Extract the base name (without .pdb extension)
    pdb_basename=$(basename "$pdb_file" .pdb)
    output_png="$(realpath "$OUTPUT_DIR")/${pdb_basename}.png"

    # Create a temporary PyMOL script using the provided command file
    temp_script="render_pdb.pml"
    sed "s|\$pdb_file|$pdb_file|g; s|\$output_png|$output_png|g" "$PYMOL_COMMAND_FILE" > "$temp_script"

    # Run PyMOL with the generated script
    echo "Processing: $pdb_file → $output_png"
    pymol -c "$temp_script" > pymol_log.txt 2>&1

    # Remove temporary script
    rm -f "$temp_script"

    # Check if PNG was created
    if [ -f "$output_png" ]; then
        echo " Saved: $output_png"
    else
        echo " Error: PNG not created for $pdb_file. Check pymol_log.txt."
    fi

done

echo "Processing complete. Images saved in $OUTPUT_DIR"
