#!/usr/bin/env python3

"""
Triple Pocket Residue Intersection Analysis
------------------------------------------

This script reads a GenBank file with pocket annotations from Fpocket (Fpockt#),
DGSite (DGsite#), and ASpock (ASpock#). It performs a triple nested loop over the
indices of these pocket annotations and for each combination (Fpockt_i,
DGsite_j, ASpock_k), it calculates how many and which residues are common among all
three pockets and writes the results to an output file:
  - Basename of the input GenBank file (Protein)
  - Maximum percentage overlap among the three pockets (max%)
  - Number of residues in each pocket (Nr_Fp, Nr_DG, Nr_Asp)
  - Score keys and values (score1_key, score1, score2_key, score2) for each pocket
  - Number of common residues
  - Percentage of common residues relative to each pocket's size
  - List of common residues at the end, formatted without commas
  - Write output in output_folder file 

Printouts occur when any percentage exceeds 0.5.

Usage:
    python triple_pocket_analysis.py input.gb output_folder

Example:
    python triple_pocket_analysis.py ADSS2.gb results/

Requirements:
    - biopython

Author: EP
Date: 2025-06-03
Version: 1.12
"""

import sys
import os
from Pocket_List_Builder_score_key import count_annotations, build_pocket_list, pocket_list_score


def main():
    # Verify arguments: input GB file and output folder
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    gb_file = sys.argv[1]
    out_dir = sys.argv[2]
    # Ensure output directory exists
    if not os.path.isdir(out_dir):
        try:
            os.makedirs(out_dir)
        except Exception as e:
            sys.stderr.write(f"Error: Could not create output directory '{out_dir}': {e}\n")
            sys.exit(1)

    # Extract basename without extension
    protein_name = os.path.splitext(os.path.basename(gb_file))[0]
    # Construct output file path
    out_path = os.path.join(out_dir, f"{protein_name}.out")

    # Define patterns for pocket annotation feature types
    patterns = ["Fpockt#", "DGsite#", "ASpock#"]

    # Count annotation features to determine max index for each pocket source
    counts = count_annotations(gb_file, patterns)
    fp_max = counts.get("Fpockt#", 0)
    dg_max = counts.get("DGsite#", 0)
    asp_max = counts.get("ASpock#", 0)

    # Build initial dictionary of pocket indices to residue lists
    pocket_list = build_pocket_list(gb_file, patterns)

    # Enrich pocket_list with scores for each source:
    pocket_list = pocket_list_score(gb_file, pocket_list, ["Fpockt#"], "Drug Score", "Pocket Score")
    pocket_list = pocket_list_score(gb_file, pocket_list, ["DGsite#"], "drugScore", "depth")
    pocket_list = pocket_list_score(gb_file, pocket_list, ["ASpock#"], "e", "v*buriedness^2/rg")

    # Open output file for writing
    try:
        out_f = open(out_path, 'w')
    except Exception as e:
        sys.stderr.write(f"Error: Could not open output file '{out_path}' for writing: {e}\n")
        sys.exit(1)

    # Write header
    header = (
        "Protein,max%,Pocket_A,Nr_Fp,Drug Score,Pocket Score,Pocket_B,Nr_DG,"
        "drugScore,depth,Pocket_C,Nr_Asp,e,v*buriedness^2/rg,common residues,%Fpockt,%DGsite,%ASpock,list"
    )
    out_f.write(header + "\n")

    # Triple nested loop over each pocket index combination
    for i in range(1, fp_max + 1):
        fp_data = pocket_list.get("Fpockt", {}).get(i, {})
        residues_fp = set(fp_data.get("residues", []))
        len_fp = len(residues_fp)
        fp_score1 = fp_data.get("score1")
        fp_score2 = fp_data.get("score2")
        
        for j in range(1, dg_max + 1):
            dg_data = pocket_list.get("DGsite", {}).get(j, {})
            residues_dg = set(dg_data.get("residues", []))
            len_dg = len(residues_dg)
            dg_score1 = dg_data.get("score1")
            dg_score2 = dg_data.get("score2")

            for k in range(1, asp_max + 1):
                asp_data = pocket_list.get("ASpock", {}).get(k, {})
                residues_asp = set(asp_data.get("residues", []))
                len_asp = len(residues_asp)
                asp_score1 = asp_data.get("score1")
                asp_score2 = asp_data.get("score2")

                # Compute common residues among the three sets
                common_residues = residues_fp & residues_dg & residues_asp
                num_common = len(common_residues)

                # Compute percentage relative to each pocket
                pct_fp = (num_common / len_fp) if len_fp > 0 else 0
                pct_dg = (num_common / len_dg) if len_dg > 0 else 0
                pct_asp = (num_common / len_asp) if len_asp > 0 else 0

                # Determine maximum percentage among the three
                max_pct = max(pct_fp, pct_dg, pct_asp)

                # Only write if any percentage > 0.5
                if pct_fp > 0.5 or pct_dg > 0.5 or pct_asp > 0.5:
                    sorted_res = sorted(common_residues)
                    res_str = "[" + " ".join(str(r) for r in sorted_res) + "]"

                    line = (
                        f"{protein_name},{max_pct:.3f},Fpockt_{i},{len_fp},{fp_score1},{fp_score2},"
                        f"DGsite_{j},{len_dg},{dg_score1},{dg_score2},"
                        f"ASpock_{k},{len_asp},{asp_score1},{asp_score2},"
                        f"{num_common},{pct_fp:.3f},{pct_dg:.3f},{pct_asp:.3f},{res_str}"
                    )
                    out_f.write(line + "\n")

    out_f.close()

if __name__ == "__main__":
    main()
