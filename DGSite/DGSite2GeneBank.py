#!/usr/bin/env python3

"""
Generate a GenBank file incorporating data from *_desc.txt and *_res_P_*.pdb files.

Usage:
    python DGSite2GeneBank.py <desc_file> <pdb_folder> <output.gb>

Example:
    python DGSite2GeneBank.py ADSS2_desc.txt ./pdb_files/ output.gb
"""

import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

def parse_pdb_residues(pdb_file):
    """Extract residue numbers from the PDB file."""
    residues = set()
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    res_num = int(line[22:26].strip())  # Extracting residue number from PDB column
                    residues.add(res_num)
                except ValueError:
                    continue
    return sorted(residues)

def parse_desc_file(desc_file):
    """Parse the *_desc.txt file to extract relevant pocket feature data."""
    with open(desc_file, 'r') as f:
        headers = f.readline().strip().split("\t")  # Read headers
        data = [line.strip().split("\t") for line in f.readlines()]

    # Convert data into a structured dictionary
    pockets = {}
    for row in data:
        pocket_id = row[0]
        pocket_data = {headers[i]: row[i] for i in range(1, len(headers))}
        pockets[pocket_id] = pocket_data

    return pockets

def create_genbank_record(pockets, pdb_files):
    """Generate a GenBank record incorporating pocket feature data and residue structures."""
    sequence = "NNNNNNNNNN"  # Placeholder sequence as PDB files don't contain sequence data
    record = SeqRecord(Seq(sequence), id="Generated_GB", name="Generated_GB",
                       description="Generated GenBank file with pocket structures.",
                       annotations={"molecule_type": "protein"})

    # Map pockets to their corresponding residue numbers from PDB files
    for index, (pocket_id, pocket_data) in enumerate(pockets.items(), start=1):
        matched_pdbs = [pdb for pdb in pdb_files if pocket_id in os.path.basename(pdb)]
        residues = []
        for pdb_file in matched_pdbs:
            residues.extend(parse_pdb_residues(pdb_file))
        residues = sorted(set(residues))

        if residues:
            qualifiers = {k: v for k, v in pocket_data.items() if k not in ["name"]}
            feature = SeqFeature(FeatureLocation(1, 2), type=f"DGsite{index:02}", qualifiers=qualifiers)
            feature.qualifiers["join"] = f"join({', '.join(map(str, residues))})"
            record.features.append(feature)

    return record

def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)

    desc_file = sys.argv[1]
    pdb_folder = sys.argv[2]
    output_file = sys.argv[3]

    print(f"Running: python {' '.join(sys.argv)}")

    if not os.path.exists(desc_file):
        print(f"Error: Description file {desc_file} not found.")
        print(__doc__)
        sys.exit(1)

    if not os.path.exists(pdb_folder):
        print(f"Error: PDB folder {pdb_folder} not found.")
        print(__doc__)
        sys.exit(1)

    pdb_files = [os.path.join(pdb_folder, f) for f in os.listdir(pdb_folder) if f.endswith(".pdb")]
    if not pdb_files:
        print("Warning: No PDB files found in the specified folder.")

    # Parse desc file and create GenBank record
    pockets = parse_desc_file(desc_file)
    record = create_genbank_record(pockets, pdb_files)

    # Save to GenBank file
    SeqIO.write(record, output_file, "genbank")
    print(f"GenBank file saved as {output_file}")

if __name__ == "__main__":
    main()
