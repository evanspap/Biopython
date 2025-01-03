#!/bin/bash

# Check for correct usage
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <Combine_Genebank.py> <folder_in1> <folder_in2> <folder_out>"
    exit 1
fi

# Get arguments
combine_script="$1"
folder_in1="$2"
folder_in2="$3"
folder_out="$4"

# Check if the Python script exists
if [ ! -f "$combine_script" ]; then
    echo "Error: Python script $combine_script not found."
    exit 1
fi

# Create the output folder if it doesn't exist
mkdir -p "$folder_out"

# Iterate over all .gb files in folder_in1
for file1 in "$folder_in1"/*.gb; do
    # Get the base filename
    base_name=$(basename "$file1")
    file2="$folder_in2/$base_name"

    # Check if the corresponding file exists in folder_in2
    if [ -f "$file2" ]; then
        # Define the output file path
        output_file="$folder_out/$base_name"

        # Run the Python script
        python "$combine_script" "$file1" "$file2" "$output_file"

        if [ $? -eq 0 ]; then
            echo "Processed: $base_name"
        else
            echo "Failed to process: $base_name"
        fi
    else
        echo "Warning: Matching file not found for $base_name in $folder_in2"
    fi
done
