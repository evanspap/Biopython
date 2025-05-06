#!/usr/bin/env python
"""
Convert an SDF file to a canonical SMILES file using RDKit.

Usage:
    python sdf_to_smi.py -i input.sdf -o output.smi [--dry-run]

Arguments:
    -i, --input      Input SDF file path
    -o, --output     Output SMILES file path
    --dry-run        Perform a dry run without actual file writing

Example:
    python sdf_to_smi.py -i molecules.sdf -o smiles.smi
"""

import argparse
from rdkit import Chem


def sdf_to_smiles(input_file, output_file, dry_run=False):
    suppl = Chem.SDMolSupplier(input_file)

    if dry_run:
        print(f"[Dry Run] Would convert '{input_file}' to canonical SMILES in '{output_file}'")
        return

    with open(output_file, 'w') as f:
        for mol in suppl:
            if mol is None:
                continue  # Skip invalid molecules
            smi = Chem.MolToSmiles(mol, canonical=True)
            mol_name = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
            f.write(f"{smi} {mol_name}\n")

    print(f"Conversion complete: '{output_file}' generated.")


def main():
    parser = argparse.ArgumentParser(description="Convert SDF to canonical SMILES format.")
    parser.add_argument('-i', '--input', required=True, help="Input SDF file")
    parser.add_argument('-o', '--output', required=True, help="Output SMILES file")
    parser.add_argument('--dry-run', action='store_true', help="Perform dry run")

    args = parser.parse_args()

    print(f"Running: sdf_to_smiles(input='{args.input}', output='{args.output}', dry_run={args.dry_run})")

    sdf_to_smiles(args.input, args.output, args.dry_run)


if __name__ == "__main__":
    main()
