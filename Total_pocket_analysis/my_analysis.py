# my_analysis.py

import sys
# Αν χρησιμοποιείς Windows CMD:
# sys.path.append(r"G:\My Drive\VSCode_Github\BioPython\Total_pocket_analysis")
# Αν χρησιμοποιείς WSL/Linux:
# sys.path.append("/mnt/g/My Drive/VSCode_Github/BioPython/Total_pocket_analysis")

from Pocket_List_Builder import count_annotations, build_pocket_list
import pprint

gb_file = "ADSS2.gb"
patterns = ["Fpockt#", "DGsite#", "ASpock#"]

counts = count_annotations(gb_file, patterns)
print("Counts:", counts)

pocket_list = build_pocket_list(gb_file, patterns)
pp = pprint.PrettyPrinter(width=120, compact=True)

pp.pprint(pocket_list)
