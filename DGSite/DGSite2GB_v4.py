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
Script: DGSite2GB_v4.py
Version: 2025-05-20
Last modified: 2025-05-20
Description: This script extracts residue features from only the `_res_P_*` pocket and subpocket PDB files,
             matches them with corresponding descriptor data from a `*_desc.txt` file (if found),
             and writes annotated features into a GenBank file. Features are named `DG_<pocket_number>`
             (e.g., DG_16 or DG_16_2).
"""

# Parse input arguments
pdb_dir = sys.argv[1] if len(sys.argv) > 1 else None
output_file = sys.argv[2] if len(sys.argv) > 2 else None
if not pdb_dir or not output_file:
    print("Usage: python DGSite2GB_v4.py <pdb_dir> <output_genbank_file>")
    print("Example: python DGSite2GB_v4.py ./Human_2e1r/ Human_2e1r.gb")
    sys.exit(1)

# Auto-detect descriptor file based on output filename
record_name = os.path.splitext(os.path.basename(output_file))[0]
possible_desc = [f"{record_name}_desc.txt", os.path.join(pdb_dir, f"{record_name}_desc.txt")]
desc_data = None
for desc_path in possible_desc:
    if os.path.exists(desc_path):
        desc_data = pd.read_csv(desc_path, sep="\t")
        break
if desc_data is None:
    print(f"Descriptor file not found among: {possible_desc}, proceeding without descriptors.")

# Initialize GenBank record
overall_sequence = ''  # Construct or load as needed
record = SeqRecord(Seq(overall_sequence), id=record_name, description='')
record.annotations['molecule_type'] = 'protein'
record.annotations['date'] = datetime.datetime.now().strftime('%Y-%m-%d')

# Process only _res_P_*.pdb pocket and subpocket files
pattern = os.path.join(pdb_dir, '*_res_P_*.pdb')
for pdb_file in sorted(glob.glob(pattern)):
    # Extract pocket/subpocket number for feature naming
    match = re.search(r'_res_P_(\d+(?:_\d+)*)\.pdb$', pdb_file)
    if not match:
        continue
    pocket_number = match.group(1)
    feature_name = f"DG_{pocket_number}"

    # Print annotated file for debugging
    print(f"Annotating {os.path.basename(pdb_file)} as {feature_name}")

    # Extract residue numbers from PDB
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    residue_numbers = [res.id[1]
                       for model in structure
                       for chain in model
                       for res in chain
                       if res.id[0] == ' ']

    # Build location
    parts = [FeatureLocation(pos-1, pos) for pos in residue_numbers]
    if not parts:
        print(f"Warning: no residues found in {pdb_file}, skipping")
        continue
    location = parts[0] if len(parts) == 1 else CompoundLocation(parts)

    # Retrieve descriptor qualifiers if available
    qualifiers = {'note': feature_name, 'source': os.path.basename(pdb_file)}
    if desc_data is not None:
        key = f"P_{pocket_number}"
        pocket_row = desc_data[desc_data['name'] == key]
        if not pocket_row.empty:
            for col in pocket_row.columns:
                qualifiers[col] = str(pocket_row.iloc[0][col])

    # Create and append feature
    feature = SeqFeature(location, type=feature_name, qualifiers=qualifiers)
    record.features.append(feature)

# Write out GenBank file
with open(output_file, 'w') as output_handle:
    SeqIO.write(record, output_handle, 'genbank')

print(f"GenBank file '{output_file}' generated successfully with {len(record.features)} features!")
