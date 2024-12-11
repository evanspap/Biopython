import os
import math
import argparse
import pandas as pd
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, CompoundLocation, FeatureLocation
from Bio.Seq import Seq
from Bio import SeqIO

# Mapping from three-letter codes to single-letter amino acid codes
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}

def parse_fasta(fasta_dir, base_name):
    """Parse the FASTA sequence from a file."""
    fasta_file = os.path.join(fasta_dir, f"{base_name}.fasta")
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"FASTA file '{fasta_file}' not found!")
    
    with open(fasta_file, 'r') as file:
        fasta_record = next(SeqIO.parse(file, "fasta"))
        return str(fasta_record.seq)

def generate_residues_from_fasta(sequence):
    """Generate residue data from the FASTA sequence."""
    residues = []
    for i, amino_acid in enumerate(sequence, start=1):
        res_name = [k for k, v in THREE_TO_ONE.items() if v == amino_acid]
        if not res_name:
            res_name = "UNK"  # Unknown residue
        else:
            res_name = res_name[0]
        # Dummy coordinates for compatibility
        residues.append(("ATOM", res_name, "A", str(i), 0.0, 0.0, 0.0))
    return residues

def parse_pocket_pdb(file_path):
    """Parse pocket ATOM coordinates from PDB file."""
    atoms = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                atoms.append((x, y, z))
    return atoms

def calculate_distance(coord1, coord2):
    """Calculate Euclidean distance between two 3D points."""
    return math.sqrt(sum((c1 - c2) ** 2 for c1, c2 in zip(coord1, coord2)))

def find_close_residues(residues, pocket_atoms, threshold):
    """Find unique residues within a threshold distance from pocket atoms."""
    close_residues = {}
    for residue in residues:
        res_seq = residue[3]  # Residue sequence number
        res_coord = residue[4:7]  # (x, y, z)
        for atom_coord in pocket_atoms:
            if calculate_distance(res_coord, atom_coord) < threshold:
                close_residues[res_seq] = residue  # Use res_seq as a key to ensure uniqueness
                break  # No need to check further for this residue
    return list(close_residues.values())  # Return unique residues for the pocket

def combine_annotations(annotations):
    """Combine ranges and attributes for each pocket, supporting non-contiguous ranges."""
    combined = {}
    for start, end, name, cluster, extras in annotations:
        if name not in combined:
            combined[name] = {
                "ranges": [],
                "cluster": cluster,
                "extras": extras
            }
        combined[name]["ranges"].append(FeatureLocation(start, end, strand=None))  # No direction

    # Convert combined ranges to GenBank features
    features = []
    for name, data in combined.items():
        # Combine ranges into a CompoundLocation if they are non-contiguous
        ranges = sorted(data["ranges"], key=lambda r: r.start)
        location = CompoundLocation(ranges) if len(ranges) > 1 else ranges[0]

        feature = {
            "location": location,
            "name": name,
            "cluster": data["cluster"],
            "extras": data["extras"]
        }
        features.append(feature)
    return features

def write_genbank(output_file, combined_annotations, sequence, base_name):
    """Write annotations and sequence to a GenBank file."""
    # Create a SeqRecord
    record = SeqRecord(
        Seq(sequence),
        id=base_name,
        name=base_name,
        description="Pocket annotations for protein",
        annotations={"molecule_type": "protein"}  # Add molecule_type
    )

    # Add combined features
    for annotation in combined_annotations:
        location = annotation["location"]
        name = annotation["name"]
        cluster = annotation["cluster"]
        extras = annotation["extras"]

        # Split extras into multiple fields
        qualifiers = {"name": name}
        for attr in extras.split(";"):
            key, value = attr.split("=")
            qualifiers[key.strip()] = value.strip()

        feature = SeqFeature(
            location=location,
            type=name,  # Use the pocket name as the type
            qualifiers=qualifiers
        )
        record.features.append(feature)

    # Write to GenBank file
    with open(output_file, "w") as file:
        SeqIO.write(record, file, "genbank")

def main():
    parser = argparse.ArgumentParser(description="Generate GenBank annotations for pocket residues.")
    parser.add_argument(
        "summary_csv",
        help="Full path to the summary CSV file (e.g., './TLE3_reduce_summary.csv')."
    )
    parser.add_argument(
        "-f", "--fasta-dir", required=True,
        help="Directory containing the FASTA file (e.g., './fasta_dir')."
    )
    parser.add_argument(
        "-t", "--threshold", type=float, default=3.0,
        help="Threshold distance for residue-pocket proximity (default: 3.0)"
    )
    parser.add_argument(
        "-o", "--output-file", default=None,
        help="Optional: Output GenBank file name. If not provided, will use '<base-name>_pockets.gbk'."
    )
    args = parser.parse_args()

    # Extract input directory and base name from summary_csv argument
    summary_csv = args.summary_csv
    input_dir = os.path.dirname(summary_csv)
    base_name = os.path.basename(summary_csv).replace("_summary.csv", "")
    threshold = args.threshold

    fasta_dir = args.fasta_dir
    output_file = args.output_file or os.path.join(input_dir, f"{base_name}_pockets.gbk")

    # Check input files
    if not os.path.exists(summary_csv):
        print(f"Error: Summary CSV file '{summary_csv}' not found!")
        return

    # Step 1: Read the FASTA sequence
    sequence = parse_fasta(fasta_dir, base_name)

    # Step 2: Generate residue data from FASTA
    residues = generate_residues_from_fasta(sequence)

    # Step 3: Read and clean pocket cluster data from CSV
    summary_data = pd.read_csv(summary_csv)

    # Ensure the column is treated as a string
    summary_data['Cluster # | Energy | #points | Radius ofN'] = summary_data['Cluster # | Energy | #points | Radius ofN'].astype(str)

    # Extract cluster numbers from the string column
    summary_data['cluster#'] = summary_data['Cluster # | Energy | #points | Radius ofN'].str.extract(r'(\d+)').astype(int)

    # Step 4: Parse and process pocket PDB files
    annotations = []
    for index, row in summary_data.iterrows():
        cluster_num = row['cluster#']
        pocket_file = os.path.join(input_dir, f"{base_name}_fp_{cluster_num:03d}.pdb")
        
        if not Path(pocket_file).exists():
            print(f"Warning: Pocket file '{pocket_file}' not found, skipping...")
            continue

        pocket_atoms = parse_pocket_pdb(pocket_file)
        close_residues = find_close_residues(residues, pocket_atoms, threshold)

        for residue in close_residues:
            _, _, _, res_seq, x, y, z = residue
            annotations.append((
                int(res_seq) - 1,  # Convert to 0-based index
                int(res_seq),      # End position
                f"Pocket{cluster_num:02d}",
                f"Cluster{cluster_num}",
                f"e={row['e']};v={row['v']};rg={row['rg']};epv={row['epv']};buriedness={row['buriedness']};v*buriedness^2/rg={row['v*buriedness^2/rg']}"
            ))

    # Step 5: Combine and write the GenBank file
    combined_annotations = combine_annotations(annotations)
    write_genbank(output_file, combined_annotations, sequence, base_name)
    print(f"GenBank file written to {output_file}")

if __name__ == "__main__":
    main()
