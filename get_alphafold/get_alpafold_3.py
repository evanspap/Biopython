import sys


print("#!/bin/sh")

f = open(sys.argv[1], "r")

for ln in f:
  #print(x)
	lnargs=ln.split()
	out="curl -o "+lnargs[1]+"_v1.pdb https://alphafold.ebi.ac.uk/files/AF-"+lnargs[0]+"-F1-model_v1.pdb"
	print(out)
	out="curl -o "+lnargs[1]+"_v2.pdb https://alphafold.ebi.ac.uk/files/AF-"+lnargs[0]+"-F1-model_v2.pdb"
	print(out)
	out="curl -o "+lnargs[1]+"_v3.pdb https://alphafold.ebi.ac.uk/files/AF-"+lnargs[0]+"-F1-model_v3.pdb"
	print(out)
	out="curl -o "+lnargs[1]+"_v4.pdb https://alphafold.ebi.ac.uk/files/AF-"+lnargs[0]+"-F1-model_v4.pdb"
	print(out)
	out="curl -o "+lnargs[1]+"_v5.pdb https://alphafold.ebi.ac.uk/files/AF-"+lnargs[0]+"-F1-model_v5.pdb"
	print(out)
