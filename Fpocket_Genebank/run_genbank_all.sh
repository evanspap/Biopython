#!/bin/bash

# Ensure the script is run with at least three arguments
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <python_script> <top_input_folder> <output_folder> [--dry-run]"
    exit 1
fi

# Assign arguments to variables
PYTHON_SCRIPT=$1
INPUT_FOLDER=$2
OUTPUT_FOLDER=$3
DRY_RUN=false

# Check for optional --dry-run flag
if [ "$#" -eq 4 ] && [ "$4" == "--dry-run" ]; then
    DRY_RUN=true
    echo "Dry run mode enabled. No commands will be executed."
fi

# Ensure the output folder exists
mkdir -p "$OUTPUT_FOLDER"

# Loop through each subfolder in the input folder
for SUBFOLDER in "$INPUT_FOLDER"/*/; do
    # Remove trailing slash and extract folder name
    SUBFOLDER_NAME=$(basename "$SUBFOLDER")

    # Check if it's a directory
    if [ -d "$SUBFOLDER" ]; then
        CMD="python \"$PYTHON_SCRIPT\" \"$SUBFOLDER\" \"$OUTPUT_FOLDER\""
        
        if [ "$DRY_RUN" = true ]; then
            echo "Dry Run: $CMD"
        else
            echo "Processing: $SUBFOLDER_NAME"
            echo "Running: $CMD"
            eval $CMD
            echo "Finished processing: $SUBFOLDER_NAME"
        fi
    fi
done

echo "All subfolders processed. GenBank files saved in '$OUTPUT_FOLDER'."
