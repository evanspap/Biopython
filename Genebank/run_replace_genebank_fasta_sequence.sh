#!/bin/bash

##############################################################
# Script: run_replace_genebank_fasta_sequence.sh
#
# Description:
#   This script replaces the sequence in a GenBank file (.gb)
#   with a provided FASTA sequence (.fa or .fasta). It requires
#   exactly three arguments:
#     1. Path to the Python script performing the replacement.
#     2. GenBank folder containing (.gbk) files.
#     3. FASTA folder containing (.fasta or .fa) files.
#
# Optional argument:
#     --dry-run : Displays commands without executing them.
#
# Usage:
#   ./run_replace_genebank_fasta_sequence.sh script.py GenBank_folder FASTA_folder [--dry-run]
#
# Example:
#   ./run_replace_genebank_fasta_sequence.sh replace_genebank_fasta_sequence.py ./genbank_files ./fasta_files --dry-run
##############################################################

# Ensure script exits on error
set -e

# Check for correct number of arguments (3 or 4 if using --dry-run)
if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 <script_path> <GenBank_folder (contains .gbk files)> <FASTA_folder (contains .fasta files)> [--dry-run]"
    exit 1
fi

# Assign arguments to variables
script_path=$1
GenBank_folder=$2
FASTA_folder=$3
dry_run=false

# Check for the dry-run flag
if [ "$#" -eq 4 ] && [ "$4" == "--dry-run" ]; then
    dry_run=true
fi

# Iterate through all .gbk files in GenBank_folder
find "$GenBank_folder" -type f -name '*.gb' | while read -r gbk_file; do
    # Extract basename without directory and suffix
    basename=$(basename "$gbk_file" ".gb")

    # Define corresponding FASTA file in FASTA_folder
    fasta_file="$FASTA_folder/${basename}.fasta"

    # Define output file path in the same directory as the GenBank file
    output_file="$(dirname "$gbk_file")/${basename}_corr.gb"

    # Check if the FASTA file exists
    if [ ! -f "$fasta_file" ]; then
        echo "FASTA file not found: $fasta_file"
        continue
    fi

    # Print operation details
    echo "Processing: $gbk_file with $fasta_file -> $output_file"

    # Perform dry run or execute the Python script
    if [ "$dry_run" = true ]; then
        echo "Dry run: Would execute -> python \"$script_path\" \"$gbk_file\" \"$fasta_file\" \"$output_file\""
    else
        python "$script_path" "$gbk_file" "$fasta_file" "$output_file"
    fi
done
