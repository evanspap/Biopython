from Bio.PDB import PDBParser

parser = PDBParser()

distance_cutoff = 3.5

target_structure = parser.get_structure("target", "target.pdb")

reference_structure = parser.get_structure("reference", "reference.pdb")

selected_residues = []

for tresidue in target_structure[0]['A']:
	print(tresidue)
	for tatom in tresidue:
		print(tatom)
		print(tatom.coord)

			
