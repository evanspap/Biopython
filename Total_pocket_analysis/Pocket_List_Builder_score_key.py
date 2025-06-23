#!/usr/bin/env python3

"""
Pocket List Builder Script
-------------------------

This script reads a GenBank file and counts the number of features
whose feature name (type) matches given prefix patterns followed by numbers
(e.g., Fpockt#, DGsite#, ASpock#). It also builds a data structure
containing the residues for each pocket annotation, and can enrich
with two score values and record the score keys for clarity.

Usage:
    python Pocket_List_Builder.py input.gb Fpockt# DGsite# ASpock#

Example:
    python Pocket_List_Builder.py ADSS2.gb Fpockt# DGsite# ASpock#

Requirements:
    - biopython

Author: EP
Date: 2025-06-03
Version: 1.5
"""

import sys
import re
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, FeatureLocation


def count_annotations(file_path, patterns):
    record = SeqIO.read(file_path, "genbank")
    counts = {p: 0 for p in patterns}
    for feature in record.features:
        feature_name = feature.type
        for pattern in patterns:
            pattern_regex = pattern.replace('#', r'\d+')
            regex = re.compile(r'^' + pattern_regex + r'$')
            if regex.match(feature_name):
                counts[pattern] += 1
    return counts


def build_pocket_list(file_path, patterns):
    record = SeqIO.read(file_path, "genbank")
    pocket_list = {}
    for pattern in patterns:
        source = pattern.rstrip('#')
        pocket_list[source] = {}

    for feature in record.features:
        fname = feature.type
        for pattern in patterns:
            source = pattern.rstrip('#')
            pattern_regex = pattern.replace('#', r'\d+')
            regex = re.compile(r'^' + pattern_regex + r'$')
            if regex.match(fname):
                idx_match = re.search(r'(\d+)$', fname)
                if not idx_match:
                    continue
                idx = int(idx_match.group(1))
                residues = []
                loc = feature.location
                parts = []
                if isinstance(loc, CompoundLocation):
                    parts = loc.parts
                elif isinstance(loc, FeatureLocation):
                    parts = [loc]
                for part in parts:
                    residues.append(int(part.start) + 1)
                pocket_list[source][idx] = {"residues": residues}
                break
    return pocket_list


def pocket_list_score(file_path, pocket_list, patterns, score1_key, score2_key):
    """
    Enrich pocket_list with score1 and score2 for each pocket, and record keys.

    Parameters:
        file_path: path to GenBank file
        pocket_list: dict returned by build_pocket_list
        patterns: list of patterns (e.g., ["Fpockt#"])
        score1_key: qualifier key for primary score (e.g., 'Drug Score')
        score2_key: qualifier key for secondary score (e.g., 'Pocket Score')

    Returns:
        Updated pocket_list where each pocket entry includes:
            'score1', 'score2', 'score1_key', 'score2_key'
    """
    record = SeqIO.read(file_path, "genbank")
    for feature in record.features:
        fname = feature.type
        for pattern in patterns:
            source = pattern.rstrip('#')
            pattern_regex = pattern.replace('#', r'\d+')
            regex = re.compile(r'^' + pattern_regex + r'$')
            if regex.match(fname):
                idx_match = re.search(r'(\d+)$', fname)
                if not idx_match:
                    continue
                idx = int(idx_match.group(1))
                s1 = None
                s2 = None
                quals = feature.qualifiers
                if score1_key in quals:
                    try:
                        s1 = float(quals[score1_key][0])
                    except (ValueError, TypeError):
                        s1 = None
                if score2_key in quals:
                    try:
                        s2 = float(quals[score2_key][0])
                    except (ValueError, TypeError):
                        s2 = None
                if source in pocket_list and idx in pocket_list[source]:
                    pocket_list[source][idx]["score1"] = s1
                    pocket_list[source][idx]["score2"] = s2
                    pocket_list[source][idx]["score1_key"] = score1_key
                    pocket_list[source][idx]["score2_key"] = score2_key
                break
    return pocket_list


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    gb_file = sys.argv[1]
    annotation_types = sys.argv[2:]

    print(f"Running: python Pocket_List_Builder.py {gb_file} {' '.join(annotation_types)}", end=',')

    result = count_annotations(gb_file, annotation_types)
    counts_line = ", ".join(f"{k}: {v}" for k, v in result.items())
    print(counts_line)

    pocket_list = build_pocket_list(gb_file, annotation_types)
    pocket_list = pocket_list_score(gb_file, pocket_list, annotation_types, annotation_types[0].rstrip('#'), annotation_types[0].rstrip('#'))  # example usage
    # Uncomment to pretty-print:
    # import pprint
    # pprint.PrettyPrinter(width=200, compact=True).pprint(pocket_list)
