"""
Script: sdf_2d_to_3d.py
Description: Converts 2D molecules in an SDF file to 3D structures using RDKit.
             Uses the ETKDGv3 method for realistic conformer generation, with UFF/MMFF optimization fallback.
Version: 1.1
Author: Evangelos Papadopoulos
Date: 2025-07-17

Example usage:
    python sdf_2d_to_3d.py input_2d.sdf output_3d.sdf
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import SDWriter
from rdkit import RDLogger

# Suppress RDKit warnings like UFFTYPER messages
RDLogger.DisableLog('rdApp.warning')

def generate_3d_sdf(input_sdf, output_sdf, max_confs=1):
    suppl = Chem.SDMolSupplier(input_sdf, removeHs=False)
    writer = SDWriter(output_sdf)

    for i, mol in enumerate(suppl):
        if mol is None:
            print(f"Skipping molecule {i} (could not be read)")
            continue

        mol = Chem.AddHs(mol)

        params = AllChem.ETKDGv3()
        params.randomSeed = 0xf00d
        success = AllChem.EmbedMolecule(mol, params)
        if success != 0:
            print(f"Embedding failed for molecule {i}")
            continue

        try:
            AllChem.UFFOptimizeMolecule(mol)
        except Exception as e:
            print(f"UFF optimization failed for molecule {i}: {e}\nTrying MMFF...")
            try:
                props = AllChem.MMFFGetMoleculeProperties(mol)
                ff = AllChem.MMFFGetMoleculeForceField(mol, props)
                ff.Minimize()
            except Exception as mmff_e:
                print(f"MMFF optimization also failed for molecule {i}: {mmff_e}")
                continue

        Chem.SanitizeMol(mol)
        writer.write(mol)

    writer.close()
    print(f"3D SDF written to: {output_sdf}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python sdf_2d_to_3d.py input_2d.sdf output_3d.sdf")
        sys.exit(1)

    input_sdf = sys.argv[1]
    output_sdf = sys.argv[2]
    generate_3d_sdf(input_sdf, output_sdf)
