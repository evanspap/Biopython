import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
import re

def parse_pdb(pdb_file):
    """Extracts unique residues from a PDB file and header information."""
    residues = set()
    header_info = {}
    
    def clean_header_key(key):
        key = key.replace("HEADER", "").strip()  # Remove "HEADER"
        key = re.sub(r'^\d+\s*-\s*', '', key)  # Remove leading numbers and dash
        key = re.sub(r'\s+\d+$', '', key)  # Remove trailing numbers
        return key

    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("HEADER") or line.startswith("REMARK"):
                parts = line.strip().split(':', 1)
                if len(parts) == 2:
                    raw_key, value = parts[0].strip(), parts[1].strip()
                    key = clean_header_key(raw_key)
                    if key.lower() not in ["information about the pocket", ""]:
                        header_info[key] = value
            if line.startswith("ATOM"):
                chain_id = line[21]
                res_name = line[17:20].strip()
                res_seq = line[22:26].strip()
                residues.add((chain_id, res_name, int(res_seq)))
    return sorted(residues, key=lambda x: (x[0], x[2])), header_info

def create_genbank(folder_name, output_folder):
    """Generates a GenBank file from PDB pocket files and saves it to a specific folder."""

    # Normalize paths
    folder_name = os.path.abspath(os.path.normpath(folder_name))
    output_folder = os.path.abspath(os.path.normpath(output_folder))

    print(f"Normalized input folder: {folder_name}")
    print(f"Normalized output folder: {output_folder}")

    features = []
    sequence = "NNNNNNNNNNNNNNNNNNNN"  # Placeholder sequence
    record = SeqRecord(Seq(sequence), id=os.path.basename(folder_name), name=os.path.basename(folder_name),
                       description=f"Generated from {folder_name} PDB pockets")
    
    record.annotations["molecule_type"] = "protein"
    record.annotations["source"] = folder_name
    
    pocket_folder = os.path.join(folder_name, "pockets")
    if not os.path.isdir(pocket_folder):
        print(f"Error: Pocket folder '{pocket_folder}' not found.")
        sys.exit(1)
    
    print(f"Looking for pocket files in folder: {pocket_folder}")
    for pdb_file in sorted(os.listdir(pocket_folder)):
        if pdb_file.startswith("pocket") and pdb_file.endswith("_atm.pdb"):
            full_path = os.path.join(pocket_folder, pdb_file)
            print(f"Processing file: {full_path}")
            residues, header_info = parse_pdb(full_path)
            
            if residues:
                pocket_number = re.search(r'pocket(\d+)_atm.pdb', pdb_file).group(1)
                feature_key = f"Pocket{pocket_number.zfill(2)}"
                residue_positions = [res[2] for res in residues]
                location = CompoundLocation([FeatureLocation(pos, pos + 1) for pos in residue_positions])
                
                cleaned_headers = {k: v for k, v in header_info.items() if k and not k.isspace()}
                
                feature = SeqFeature(location, type=feature_key, qualifiers=cleaned_headers)
                features.append(feature)
            else:
                print(f"No residues found in {pdb_file}")
    
    record.features = features

    # Ensure output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Save GenBank file in the output folder
    output_file = os.path.join(output_folder, f"{os.path.basename(folder_name)}.gb")
    with open(output_file, "w") as gb_file:
        SeqIO.write(record, gb_file, "genbank")
    
    print(f"GenBank file created: {output_file} with {len(features)} features.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_folder> <output_folder>")
        print("Example: python script.py ADSS2 GeneBank")
        sys.exit(1)
    
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]

    create_genbank(input_folder, output_folder)
