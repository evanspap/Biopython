from Bio.PDB import PDBParser

parser = PDBParser()

distance_cutoff = 3.5

target_structure = parser.get_structure("target", "target.pdb")

reference_structure = parser.get_structure("reference", "reference.pdb")

selected_residues = []

for tmodel in target_structure:
	print(tmodel)
	for tchain in tmodel:
		print(tchain)
		for tresidue in tchain:
			print(tresidue)
			for tatom in tresidue:
				print(tatom)
				print(tatom.coord)

			
