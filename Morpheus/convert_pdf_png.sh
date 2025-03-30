#!/bin/bash

# =====================================
# Script: Convert all PDFs to PNGs in a given folder
# Usage: ./convert_pdf_png.sh <folder> [--dry-run]
# Example: ./convert_pdf_png.sh /home/user/documents/pdfs
# Requires: ImageMagick installed (for the 'convert' command)
# =====================================

# Check for required argument
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <folder> [--dry-run]"
    echo "Example: $0 /path/to/pdf_folder"
    exit 1
fi

FOLDER="$1"
DRY_RUN=false

# Check for optional dry-run flag
if [ "$2" == "--dry-run" ]; then
    DRY_RUN=true
fi

# Check if folder exists
if [ ! -d "$FOLDER" ]; then
    echo "Error: '$FOLDER' is not a directory"
    exit 1
fi

# Loop through all PDF files in the folder
for pdf in "$FOLDER"/*.pdf; do
    # Skip if no PDF files found
    [ -e "$pdf" ] || continue

    # Get filename without extension
    filename=$(basename "$pdf" .pdf)

    # Output PNG file name
    output="${FOLDER}/${filename}.png"

    # Show the command to be run
    echo "convert -density 150 -depth 8 -quality 85 \"$pdf\" \"$output\""

    # Run the command if not dry-run
    if [ "$DRY_RUN" = false ]; then
        convert -density 150 -depth 8 -quality 85 "$pdf" "$output"
    fi

done

echo "Conversion complete."
