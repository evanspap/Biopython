'''
smiles_to_sdf.py

Author: Evangelos Papadopoulos
Version: 1.1
Date: 2025-06-18

Summary:
  This script reads a text file of SMILES strings, generates 3D coordinates for each molecule using RDKit,
  and writes the resulting molecules to an output SDF file. Hydrogens are added, 3D coordinates are embedded,
  and energy minimization is performed with UFF.

Usage:
  python smiles_to_sdf.py input_smiles.txt output.sdf
  python smiles_to_sdf.py input_smiles.txt --basename
  python smiles_to_sdf.py input_smiles.txt --std-out

Example:
  python smiles_to_sdf.py CHEMBL506258.smi --basename
  python smiles_to_sdf.py CHEMBL506258.smi --std-out
'''

import sys
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import SDWriter

if len(sys.argv) < 2 or len(sys.argv) > 3:
    print("Usage: python smiles_to_sdf.py input_smiles.txt [output.sdf | --basename | --std-out]")
    sys.exit(1)

input_file = sys.argv[1]
output_mode = sys.argv[2] if len(sys.argv) == 3 else None

# Determine output destination
if output_mode == '--basename':
    base = os.path.splitext(os.path.basename(input_file))[0]
    output_file = f"{base}.sdf"
    writer = SDWriter(output_file)
    output_type = "file"
elif output_mode == '--std-out':
    import sys
    writer = SDWriter(sys.stdout)
    output_type = "stdout"
elif output_mode:  # assumed to be a filename
    output_file = output_mode
    writer = SDWriter(output_file)
    output_type = "file"
else:
    print("Output mode not specified correctly. Use output filename, --basename, or --std-out.")
    sys.exit(1)

with open(input_file, 'r') as f:
    for i, line in enumerate(f):
        smi = line.strip()
        if not smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"Skipping invalid SMILES on line {i+1}: {smi}", file=sys.stderr)
            continue

        mol = Chem.AddHs(mol)
        success = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        if success != 0:
            print(f"Embedding failed for line {i+1}: {smi}", file=sys.stderr)
            continue
        AllChem.UFFOptimizeMolecule(mol)

        mol.SetProp("_Name", f"mol_{i+1}")
        writer.write(mol)

writer.close()
if output_type == "file":
    print(f"Finished writing all molecules to {output_file}")
