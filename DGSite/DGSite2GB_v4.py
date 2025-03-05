import os
import glob
import sys
import re
import datetime
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.PDB import PDBParser

"""
Script: DGSite2GB_v2.py
Description: This script extracts residue features from PDB files and matches them with corresponding pocket data from a *_desc.txt file.
             The extracted information is stored in a GenBank file for further analysis.
             It utilizes the output from the Differential of Gaussians (DoG) pocket identification method.
Usage: python script.py <input_folder> <output_file>
Example: python script.py ./input data_output.gb
"""

# Timestamp of script creation
script_written_time = "2025-02-10 14:30:00"  # Update this manually when editing
print(f"Script last modified on: {script_written_time}")

def print_usage():
    print("\nThis script processes PDB files to extract residue features and associates them with pocket properties from a *_desc.txt file.")
    print("The output is stored in a GenBank format file. It utilizes the output from the Differential of Gaussians (DoG) pocket identification method.")
    print("\nUsage: python script.py <input_folder> <output_file>")
    print("Example: python script.py ./input data_output.gb\n")
    sys.exit(1)

# Check for command-line arguments
if len(sys.argv) != 3:
    print("\n[ERROR] Incorrect number of arguments provided.")
    print_usage()

input_folder = sys.argv[1]
output_file = sys.argv[2]

# Function to extract residue numbers from PDB
def extract_residue_numbers_from_pdb(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    
    residue_numbers = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == " ":  # Ignore heteroatoms
                    residue_numbers.append(residue.id[1])
    
    return sorted(set(residue_numbers))  # Return unique sorted residue numbers

# Find the *_desc.txt file
desc_file = glob.glob(os.path.join(input_folder, "*_desc.txt"))
if not desc_file:
    print("Warning: No *_desc.txt file found. Proceeding without additional pocket data.")
    desc_data = None
else:
    desc_file = desc_file[0]
    desc_data = pd.read_csv(desc_file, sep="\t")

# Find all matching PDB files
pdb_files = glob.glob(os.path.join(input_folder, "*_res_P_*.pdb"))

# Create a single GenBank record
record = SeqRecord(
    Seq(""),  # No sequence, only features
    id="Generated_Protein",
    name="Generated_Protein",
    description="GenBank file with extracted residue features from multiple PDB files",
    annotations={"molecule_type": "protein", "script_written_time": script_written_time}  # Required for GenBank format
)

# Add each PDB file's residue numbers as a separate feature
for pdb_file in pdb_files:
    residue_numbers = extract_residue_numbers_from_pdb(pdb_file)
    if not residue_numbers:
        continue  # Skip files with no residues
    
    # Extract pocket number for feature naming
    match = re.search(r'_res_P_(\d+)\.pdb$', pdb_file)
    pocket_number = match.group(1) if match else "XX"
    feature_name = f"DGsite{pocket_number.zfill(2)}"
    
    # Define feature location using CompoundLocation
    location = CompoundLocation([FeatureLocation(pos-1, pos) for pos in residue_numbers])
    
    # Retrieve additional data from *_desc.txt if available
    qualifiers = {"note": feature_name, "source": os.path.basename(pdb_file)}
    if desc_data is not None:
        pocket_row = desc_data[desc_data['name'] == f"P_{pocket_number}"]
        if not pocket_row.empty:
            for col in pocket_row.columns:
                qualifiers[col] = str(pocket_row.iloc[0][col])
    
    feature = SeqFeature(
        location,
        type=feature_name,
        qualifiers=qualifiers
    )
    record.features.append(feature)

# Write to GenBank file
with open(output_file, "w") as output_handle:
    SeqIO.write(record, output_handle, "genbank")

print(f"GenBank file '{output_file}' generated successfully with {len(record.features)} features!")
