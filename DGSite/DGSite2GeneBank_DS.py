#!/usr/bin/env python3
"""
Script: Generate GenBank File with DGPock Features
Description:
    This script processes a `*_desc.txt` file and multiple `*_res_P_*.pdb` files to generate a GenBank file.
    The GenBank file includes `DGPock` features for each pocket, formatted exactly like the provided `ADSS2_DS.txt`.

Usage:
    python3 generate_genbank.py <input_folder> <output_file>

Arguments:
    <input_folder>  Path to the folder containing `*_desc.txt` and `*_res_P_*.pdb` files.
    <output_file>   Path to the output GenBank file.

Example:
    python3 generate_genbank.py ./dogsite_output ./output.gb
"""

import os
import glob
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

def parse_desc_file(desc_file):
    """Parse the *_desc.txt file to extract pocket features."""
    features = {}
    with open(desc_file, 'r') as file:
        lines = file.readlines()
        headers = lines[0].strip().split('\t')
        for line in lines[1:]:
            values = line.strip().split('\t')
            pocket_id = values[0]  # e.g., "P_1"
            features[pocket_id] = dict(zip(headers[1:], values[1:]))
    return features

def parse_pdb_files(pdb_files):
    """Parse the *_res_P_*.pdb files to extract residue numbers."""
    pocket_residues = {}
    for pdb_file in pdb_files:
        # Extract pocket ID from the filename (e.g., "P_1" from "*_res_P_1.pdb")
        pocket_id = os.path.basename(pdb_file).split('_')[-1].split('.')[0]
        residues = []
        with open(pdb_file, 'r') as file:
            for line in file:
                if line.startswith("ATOM"):
                    residue_number = int(line[22:26].strip())
                    if residue_number not in residues:
                        residues.append(residue_number)
        pocket_residues[pocket_id] = sorted(residues)
    return pocket_residues

def create_genbank_record(desc_features, pocket_residues, output_file):
    """Create a GenBank file with DGPock features in the specified format."""
    # Create a dummy sequence (since we don't have a real sequence)
    sequence = Seq("N" * 1000)  # Adjust the length as needed
    record = SeqRecord(sequence, id="ADSS2", name="ADSS2", description="Generated from DoGSite 2.0.0 pockets")

    # Add molecule_type annotation (required for GenBank format)
    record.annotations["molecule_type"] = "protein"

    # Add DGPock features for each pocket
    for pocket_id, residues in pocket_residues.items():
        # Get the corresponding features from the desc file
        features = desc_features.get(pocket_id, {})

        # Create the feature location string
        location = f"join({', '.join(map(str, residues))})"

        # Create the qualifiers dictionary
        qualifiers = {
            "note": [location],
            **features  # Add all other features from the desc file
        }

        # Create the DGPock feature
        feature = SeqFeature(FeatureLocation(0, 0), type="DGPock", qualifiers=qualifiers)
        record.features.append(feature)

    # Write the GenBank file in the specified format
    with open(output_file, "w") as output_handle:
        # Write the header
        output_handle.write("FEATURES             Location/Qualifiers\n")

        # Write each DGPock feature
        for i, feature in enumerate(record.features, start=1):
            pocket_id = f"DGPock_{i:02d}"  # Format as DGPock_01, DGPock_02, etc.
            output_handle.write(f"     {pocket_id}      {feature.qualifiers['note'][0]}\n")
            for key, value in feature.qualifiers.items():
                if key != "note":
                    output_handle.write(f"                    {key}=\"{value[0]}\"\n")

def main(input_folder, output_file):
    """Main function to process input files and generate the GenBank file."""
    desc_file = glob.glob(os.path.join(input_folder, "*_desc.txt"))[0]
    pdb_files = glob.glob(os.path.join(input_folder, "*_res_P_*.pdb"))

    desc_features = parse_desc_file(desc_file)
    pocket_residues = parse_pdb_files(pdb_files)

    create_genbank_record(desc_features, pocket_residues, output_file)

if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    input_folder = sys.argv[1]
    output_file = sys.argv[2]

    # Validate input folder
    if not os.path.isdir(input_folder):
        print(f"Error: Input folder '{input_folder}' does not exist.")
        sys.exit(1)

    # Run the script
    main(input_folder, output_file)
    print(f"GenBank file successfully created: {output_file}")