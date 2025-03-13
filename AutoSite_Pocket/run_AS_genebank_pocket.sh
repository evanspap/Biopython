#!/bin/bash

# Ensure the script is invoked with the required arguments
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <path_to_genbank_script> <directory_to_search>"
    echo "Example: $0 ~/Biopython/Pocket_Analysis/AS_genbank_pocket.py ./data "
    exit 1
fi

# Parse arguments
GENBANK_SCRIPT_PATH="$1"
SEARCH_DIR="$2"

# Check if the GenBank script exists
if [[ ! -f "$GENBANK_SCRIPT_PATH" ]]; then
    echo "Error: GenBank script not found at $GENBANK_SCRIPT_PATH"
    exit 1
fi

# Check if the search directory exists
if [[ ! -d "$SEARCH_DIR" ]]; then
    echo "Error: Search directory not found at $SEARCH_DIR"
    exit 1
fi

# Iterate over all *_summary.csv files in the specified directory
for summary_file in "$SEARCH_DIR"/*/*_summary.csv; do
    if [[ -f "$summary_file" ]]; then
        # Extract the base name (e.g., ALG3_reduce)
        base_name=$(basename "$summary_file" _summary.csv)
        dir_name=$(dirname "$summary_file")

        # Generate the output GenBank file path
        output_file="${dir_name}/${base_name}_45nm.gb"

        echo "Processing $summary_file..."
        
        # Run the GenBank script
        python "$GENBANK_SCRIPT_PATH" "$summary_file" -o "$output_file"

        if [[ $? -ne 0 ]]; then
            echo "Error: Failed to process $summary_file"
        else
            echo "Output written to $output_file"
        fi
    fi
done

echo "All GenBank files generated!"
