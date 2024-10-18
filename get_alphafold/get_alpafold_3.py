import sys

f = open(sys.argv[1], "r")

for ln in f:
  #print(x)
	lnargs=ln.split()
	out="curl -o "+lnargs[1]+".pdb https://alphafold.ebi.ac.uk/files/AF-"+lnargs[0]+"-F1-model_v1.pdb"
	out="curl -o "+lnargs[1]+".pdb https://alphafold.ebi.ac.uk/files/AF-"+lnargs[0]+"-F1-model_v2.pdb"
	out="curl -o "+lnargs[1]+".pdb https://alphafold.ebi.ac.uk/files/AF-"+lnargs[0]+"-F1-model_v3.pdb"
	out="curl -o "+lnargs[1]+".pdb https://alphafold.ebi.ac.uk/files/AF-"+lnargs[0]+"-F1-model_v4.pdb"
	out="curl -o "+lnargs[1]+".pdb https://alphafold.ebi.ac.uk/files/AF-"+lnargs[0]+"-F1-model_v5.pdb"
	print(out)
