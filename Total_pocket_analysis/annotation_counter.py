#!/usr/bin/env python3

"""
Annotation Counter Script
-------------------------

This script reads a GenBank file and counts the number of features
whose feature name (type) matches given prefix patterns followed by numbers
(e.g., Fpockt#, DGsite#, ASpock#).

Usage:
    python annotation_counter.py input.gb Fpockt# DGsite# ASpock#

Example:
    python annotation_counter.py ADSS2.gb Fpockt# DGsite# ASpock#

Requirements:
    - biopython

Author: EP
Date: 2025-06-03
Version: 1.3
"""

import sys
import re
from Bio import SeqIO

def count_annotations(file_path, patterns):
    record = SeqIO.read(file_path, "genbank")
    counts = {p: 0 for p in patterns}
    for feature in record.features:
        feature_name = feature.type
        for pattern in patterns:
            # replace '#' with '\d+' to match any number
            pattern_regex = pattern.replace('#', r'\d+')
            regex = re.compile(r'^' + pattern_regex + r'$')
            if regex.match(feature_name):
                counts[pattern] += 1
    return counts

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    gb_file = sys.argv[1]
    annotation_types = sys.argv[2:]

    # Print running and keep cursor on same line
    print(f"Running: python annotation_counter.py {gb_file} {' '.join(annotation_types)}", end=',')

    result = count_annotations(gb_file, annotation_types)

    # Print all counts immediately after
    counts_line = ", ".join(f"{k}: {v}" for k, v in result.items())
    print(counts_line)
