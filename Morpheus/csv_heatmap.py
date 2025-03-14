#!/usr/bin/env python3

"""
Create a heatmap from a CSV file with plot styling from a JSON configuration.

Usage:
    python csv_heatmap.py <input.csv> <config.json> <output.pdf>

Example:
    python csv_heatmap.py data.csv config.json heatmap.pdf
"""

import sys
import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def usage():
    print(__doc__)
    sys.exit(1)

def main():
    if len(sys.argv) != 4:
        usage()

    csv_file = sys.argv[1]
    json_file = sys.argv[2]
    output_file = sys.argv[3]

    print(f"Running: python {sys.argv[0]} {csv_file} {json_file} {output_file}")

    # Load CSV data, handling non-numeric entries as NaN
    data = pd.read_csv(csv_file, index_col=0, na_values=['-', 'NA', ''])
    data = data.astype(float)

    # Load JSON configuration
    with open(json_file, 'r') as jf:
        config = json.load(jf)

    # Extract JSON configuration (example structure)
    cmap = config.get('colormap', 'viridis')
    title = config.get('title', 'Heatmap')
    xlabel = config.get('xlabel', '')
    ylabel = config.get('ylabel', '')
    annot = config.get('annot', False)
    figsize = tuple(config.get('figsize', [10, 8]))

    # Create heatmap
    plt.figure(figsize=figsize)
    sns.heatmap(data, cmap=cmap, annot=annot)

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.tight_layout()
    plt.savefig(output_file)
    print(f"Heatmap successfully saved to '{output_file}'.")

if __name__ == "__main__":
    main()
