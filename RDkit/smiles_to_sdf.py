import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import SDWriter

if len(sys.argv) != 3:
    print("Usage: python smiles_to_sdf.py input_smiles.txt output.sdf")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

writer = SDWriter(output_file)

with open(input_file, 'r') as f:
    for i, line in enumerate(f):
        smi = line.strip()
        if not smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"Skipping invalid SMILES on line {i+1}: {smi}")
            continue

        mol = Chem.AddHs(mol)
        success = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        if success != 0:
            print(f"Embedding failed for line {i+1}: {smi}")
            continue
        AllChem.UFFOptimizeMolecule(mol)

        mol.SetProp("_Name", f"mol_{i+1}")
        writer.write(mol)

writer.close()
print(f"Finished writing all molecules to {output_file}")
