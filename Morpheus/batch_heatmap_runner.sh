#!/bin/bash

# Script to generate heatmaps for all CSV files in a directory

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <path_to_csv_heatmap.py> <csv_directory> <config.json> <output_directory>"
    echo "Example: $0 ~/scripts/csv_heatmap.py ./csv_files config.json ./heatmaps"
    exit 1
fi

heatmap_script=$1
csv_dir=$2
config_json=$3
output_dir=$4

# Check if output directory exists, if not create it
mkdir -p "$output_dir"

# Loop through all CSV files in the given directory
for csv_file in "$csv_dir"/*.csv; do
    filename=$(basename "$csv_file" .csv)
    output_file="$output_dir/${filename}_heatmap.pdf"

    echo "Running: python '$heatmap_script' '$csv_file' '$config_json' '$output_file'"
    python "$heatmap_script" "$csv_file" "$config_json" "$output_file"
done


echo "All heatmaps generated in '$output_dir'."
