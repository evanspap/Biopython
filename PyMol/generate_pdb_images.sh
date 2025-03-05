#!/bin/bash

# Script: generate_pdb_images.sh
# Description: Generates images for all PDB files in the current directory using PyMOL.
# Usage:
#   chmod +x generate_pdb_images.sh
#   ./generate_pdb_images.sh

# Ensure PyMOL is installed
if ! command -v pymol &> /dev/null; then
    echo "Error: PyMOL is not installed or not in PATH."
    exit 1
fi

# Output directory for images
OUTPUT_DIR="./images"
mkdir -p "$OUTPUT_DIR"

# Ensure the output directory is writable
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

    # Create a temporary PyMOL script with white background and cartoon representation
    temp_script="render_pdb.pml"
    cat > "$temp_script" << EOF
load $pdb_file
hide everything
show cartoon
orient
zoom visible
bg_color white
png $output_png, width=1600, height=900, dpi=300
quit
EOF

    # Run PyMOL with the generated script
    echo "Processing: $pdb_file → $output_png"
    pymol -c "$temp_script" > pymol_log.txt 2>&1

    # Remove temporary script
    rm -f "$temp_script"

    # Check if PNG was created
    if [ -f "$output_png" ]; then
        echo "✅ Saved: $output_png (White Background, Cartoon Representation, 1600x900)"
    else
        echo "❌ Error: PNG not created for $pdb_file. Check pymol_log.txt."
    fi

done

echo "Processing complete. Images saved in $OUTPUT_DIR"
