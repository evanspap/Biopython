from Bio.PDB import PDBParser

parser = PDBParser()

distance_cutoff = 3

target_structure = parser.get_structure("target", "target.pdb")

reference_structure = parser.get_structure("reference", "reference.pdb")

selected_residues = []

for tmodel in target_structure:
	for tchain in tmodel:
		for tresidue in tchain:
			print(tresidue)
			for tatom in tresidue:
#				print(model, chain,residue,atom)
#				print(atom.coord)
				for rmodel in reference_structure:
					for rchain in rmodel:
						for rresidue in tchain:
							for ratom in rresidue:
								if tatom-ratom <= distance_cutoff:
									print(tatom-ratom)
									print(tresidue)
									selected_residues.append(tresidue)
									# The argument to the function may be any descriptive text 
									input("Press the Enter key to continue: ") 

									break

print(selected_residues)

