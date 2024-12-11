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
