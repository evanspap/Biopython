#!/bin/bash

# Loop through all *.pdb files in the current directory
for pdb_file in *.pdb; do
    # Check if the file exists and is a regular file
    if [ -f "$pdb_file" ]; then
        # Run the Python script on the current pdb file
        python pdb_to_genebank.py "$pdb_file"
        echo "Processed: $pdb_file"
    fi
done
