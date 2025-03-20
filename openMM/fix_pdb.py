#!/usr/bin/env python3
"""
Script Name: fix_pdb.py

Description:
  This script uses PDBFixer to fix a given PDB file. It can:
    1) Identify missing residues
    2) Identify missing atoms
    3) Add missing atoms
    4) Add missing hydrogens
  and then write out a new "fixed" PDB file for subsequent usage.

Usage:
  ./fix_pdb.py input_file.pdb output_file.pdb

Example:
  ./fix_pdb.py 4TPW_Full_noTER.pdb 4TPW_fixed.pdb
"""

import sys
import os
from pdbfixer import PDBFixer
from openmm.app import PDBFile


def main():
    # Check arguments
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    input_pdb = sys.argv[1]
    output_pdb = sys.argv[2]

    # Check for file existence
    if not os.path.isfile(input_pdb):
        print(f"Error: Input file '{input_pdb}' not found.")
        sys.exit(1)

    print(f"Running: Fixing {input_pdb} -> {output_pdb}")

    # Load your PDB
    fixer = PDBFixer(filename=input_pdb)

    # Add missing residues, atoms, and hydrogens (as needed)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)  # or adjust pH as desired

    # Write out a fixed PDB
    with open(output_pdb, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)


if __name__ == "__main__":
    main()
