from Bio.PDB import PDBParser
import sys

parser = PDBParser()

distance_cutoff = 3

target_structure = parser.get_structure("target", sys.argv[1])

reference_structure = parser.get_structure("reference", sys.argv[2])

selected_residues = []

def Check_Residue(Res):
	for tatom in Res:
		#print(tatom)
		#print(tatom.coord)
		for ratom in reference_structure[0][' '][1]:
			#print(ratom)
			#print(ratom.coord)
			if tatom-ratom <= distance_cutoff:
				#print(tatom,tatom.coord)
				#print(ratom,ratom.coord)
				#print(tatom-ratom)
				#print(tresidue)
				selected_residues.append(Res)
				return
		

for tresidue in target_structure[0]['A']:
	print(tresidue)
	Check_Residue(tresidue)
	

for i in range(len(selected_residues)):
	print(selected_residues[i].resname,selected_residues[i]._id[1])
	
for i in range(len(selected_residues)):
	print(selected_residues[i]._id[1], end="+")

print()
	
		
"""or ratom in reference_structure[0][' '][1]:
	print(ratom)
	print(ratom.coord)"""
	

			
