#!/usr/bin/env python3
"""
parse_vina_out.py

Extracts docking mode data from AutoDock Vina output files and prints it in a CSV-style format.

Usage:
    python parse_vina_out.py <vina_output_file> [--header] [--top X]

Example:
    python parse_vina_out.py CHEMBL3680879_BACE1_2ZHR_aln6UJ0.out --header --top 3

Outputs:
    Ligand, avg_affinity_all, avg_affinity_topX, mode_1, affinity_1(Kcal/mol), rmsd_lb_1, rmsd_ub_1, ...
    CHEMBL3680879, -9.89, -10.24, 1, -10.37, 0.0, 0.0, 2, -10.3, 3.469, 4.814, ...

Author: Evangelos Papadopoulos
Version: 1.1
Date: 2025-07-16
"""

import sys
import os

def parse_vina_output(filename):
    ligand_id = os.path.basename(filename).split('_')[0]
    data_started = False
    results = []

    with open(filename, 'r') as f:
        for line in f:
            if line.strip().startswith("mode |"):
                data_started = True
                continue
            if data_started:
                if line.strip() == "" or line.startswith("-----"):
                    continue
                parts = line.split()
                if len(parts) == 4:
                    mode, affinity, rmsd_lb, rmsd_ub = parts
                    results.append((int(mode), float(affinity), float(rmsd_lb), float(rmsd_ub)))
    return ligand_id, results

def format_output(ligand_id, results, top_n):
    affinities = [affinity for _, affinity, _, _ in results]
    avg_all = sum(affinities) / len(affinities) if affinities else 0
    top_avg = sum(sorted(affinities)[:top_n]) / min(top_n, len(affinities)) if affinities else 0

    flat_values = [ligand_id, round(avg_all, 3), round(top_avg, 3)]
    for mode, affinity, rmsd_lb, rmsd_ub in results:
        flat_values.extend([mode, affinity, rmsd_lb, rmsd_ub])
    return ', '.join(map(str, flat_values))

def print_header(max_modes=9):
    header = ["Ligand", "avg_affinity_all", "avg_affinity_topX"]
    for i in range(1, max_modes + 1):
        header.extend([
            f"mode_{i}",
            f"affinity_{i}(Kcal/mol)",
            f"rmsd_lb_{i}",
            f"rmsd_ub_{i}"
        ])
    print(', '.join(header))

def main():
    if len(sys.argv) < 2:
        print("Usage: parse_vina_out.py <vina_output_file> [--header] [--top X]")
        sys.exit(1)

    show_header = '--header' in sys.argv
    top_n = 3  # default value
    for i, arg in enumerate(sys.argv):
        if arg == '--top' and i + 1 < len(sys.argv):
            try:
                top_n = int(sys.argv[i + 1])
            except ValueError:
                print("Error: Invalid value for --top. Must be an integer.")
                sys.exit(1)

    vina_file = next((arg for arg in sys.argv[1:] if not arg.startswith('--')), None)

    if not vina_file or not os.path.isfile(vina_file):
        print(f"Error: file not found or not specified: {vina_file}")
        sys.exit(1)

    ligand_id, results = parse_vina_output(vina_file)

    if show_header:
        print_header(max_modes=len(results))

    print(format_output(ligand_id, results, top_n))

if __name__ == '__main__':
    main()
