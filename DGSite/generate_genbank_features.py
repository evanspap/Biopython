"""
Generate GenBank File with Pocket Features

This script reads all PDB files matching the pattern *_res_P_*.pdb, extracts their
residue numbers using Biopython, and generates a GenBank file where each pocket
is represented as a feature with joined residue numbers.

Usage:
    python generate_genbank_features.py <output_genbank_file>

Example:
    python generate_genbank_features.py ADSS2_pockets.gb
"""

import os
import sys
import glob
from Bio import PDB, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

def extract_residue_numbers(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    residue_numbers = set()
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue):  # Check if it's an amino acid
                    residue_numbers.add(residue.get_id()[1])
    
    return sorted(residue_numbers)

def create_genbank(pdb_files, output_file):
    record = SeqRecord(Seq(""), id="Pockets", description="Extracted pockets from PDB files",
                        annotations={"molecule_type": "protein"})
    
    features = []
    for pdb_file in pdb_files:
        pocket_id = os.path.basename(pdb_file).split("_res_P_")[-1].split(".")[0]
        pocket_name = f"DGpock_{pocket_id.zfill(2)}"
        residues = extract_residue_numbers(pdb_file)
        
        if residues:
            location = FeatureLocation(0, 0)  # Dummy location to prevent Biopython errors
            qualifiers = {"label": pocket_name, "residues": f"join({','.join(map(str, residues))})"}
            feature = SeqFeature(location, type="misc_feature", qualifiers=qualifiers)
            features.append(feature)
    
    record.features = features
    SeqIO.write(record, output_file, "genbank")
    print(f"GenBank file saved as {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("\nGenerate GenBank File with Pocket Features")
        print("\nUsage: python generate_genbank_features.py <output_genbank_file>\n")
        print("Example: python generate_genbank_features.py ADSS2_pockets.gb\n")
        sys.exit(1)
    
    output_file = sys.argv[1]
    pdb_files = glob.glob("*_res_P_*.pdb")
    
    if not pdb_files:
        print("No matching PDB files found.")
        sys.exit(1)
    
    create_genbank(pdb_files, output_file)
