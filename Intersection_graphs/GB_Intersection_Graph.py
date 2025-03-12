import pandas as pd
import argparse
import os

# ==============================================================
# Script: Intersection Graph Table Generator
# Description: This script processes a GenBank (.gb) file, extracts
#              residue positions from all annotations, and generates
#              an intersection graph table showing shared residues.
# Usage: python GB_Intersection_Graph.py <input_file> <output_file>
# Example: python GB_Intersection_Graph.py ADSS2.gb results/Intersection_Graph_Table.csv
# ==============================================================

# Function to extract residue positions from a join() field
def extract_positions(line):
    positions = line.split("join(")[-1].replace(")", "").split(",")
    return set(map(int, filter(str.isdigit, positions)))

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Generate an intersection graph table for all annotation residues.")
parser.add_argument("input_file", help="Path to the GenBank input file")
parser.add_argument("output_file", help="Full path to the output CSV file (including filename)")
args = parser.parse_args()

genbank_file = args.input_file
output_file = args.output_file


pockets = {}
print("Processing file:", genbank_file)

with open(genbank_file, "r") as file:
    for line in file:
        line = line.strip()
        if "join(" in line:
            pocket_name = line.split()[0]  # First word is the pocket name
            positions = extract_positions(line)
            pockets[pocket_name] = positions

# Create an intersection graph table with total residue counts
pocket_names = list(pockets.keys())
total_residues = {p: len(pockets[p]) for p in pocket_names}

intersection_matrix = {p1: {p2: 0 for p2 in pocket_names} for p1 in pocket_names}

# Compute intersection counts
for p1 in pocket_names:
    for p2 in pocket_names:
        if p1 == p2:
            intersection_matrix[p1][p2] = total_residues[p1]  # Self-comparison should be total residues
        else:
            common_residues = pockets[p1] & pockets[p2]
            intersection_matrix[p1][p2] = len(common_residues)

# Convert matrix into a DataFrame and add total residue counts
intersection_df = pd.DataFrame(intersection_matrix).T
intersection_df.insert(0, "Total Residues", [total_residues[p] for p in pocket_names])

# Add the "Total Residues" row at the top
total_residue_row = ["-"] + [total_residues[p] for p in pocket_names]
intersection_df.loc["Total Residues"] = total_residue_row

# Reorder to move "Total Residues" to the first row
intersection_df = intersection_df.iloc[[-1] + list(range(len(intersection_df) - 1))]

# Save the output to the specified file
intersection_df.to_csv(output_file)

print(f"Intersection graph table saved to {output_file}")
