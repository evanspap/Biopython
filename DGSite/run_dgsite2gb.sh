#!/bin/bash


echo "Script: run_dgsite2gb.sh"
echo "Description: This script scans first-level subfolders and runs the specified Python script on each folder."
echo "It extracts features from PDB files and outputs the results to the specified output folder."
echo "Usage: ./run_dgsite2gb.sh <python_script> <parent_folder> <output_folder> [--dry-run]"
echo "Example: ./run_dgsite2gb.sh /path/to/DGSite2GB_v2.py /data/folders /output/folder"

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <python_script> <parent_folder> <output_folder> [--dry-run]"
    exit 1
fi

PYTHON_SCRIPT="$1"
PARENT_FOLDER="$2"
OUTPUT_FOLDER="$3"
DRY_RUN=false

if [ "$4" == "--dry-run" ]; then
    DRY_RUN=true
    echo "Dry run mode enabled. Commands will be displayed but not executed."
fi

if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "Error: Python script '$PYTHON_SCRIPT' not found."
    exit 1
fi

if [ ! -d "$PARENT_FOLDER" ]; then
    echo "Error: Parent folder '$PARENT_FOLDER' not found."
    exit 1
fi

if [ ! -d "$OUTPUT_FOLDER" ]; then
    echo "Creating output folder: $OUTPUT_FOLDER"
    mkdir -p "$OUTPUT_FOLDER"
fi

for subfolder in "$PARENT_FOLDER"/*/; do
    if [ -d "$subfolder" ]; then
        subfolder_name=$(basename "$subfolder")
        output_file="$OUTPUT_FOLDER/$subfolder_name.gb"
        
        CMD="python3 \"$PYTHON_SCRIPT\" \"$subfolder\" \"$output_file\""
        
        if [ "$DRY_RUN" == true ]; then
            echo "Dry Run: $CMD"
        else
            echo "Running: $CMD"
            eval "$CMD"
        fi
    fi
done

echo "Processing complete."
