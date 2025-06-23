# my_analysis_2.py

import sys
# Αν χρησιμοποιείς Windows CMD:
# sys.path.append(r"G:\My Drive\VSCode_Github\BioPython\Total_pocket_analysis")
# Αν χρησιμοποιείς WSL/Linux:
# sys.path.append("/mnt/g/My Drive/VSCode_Github/BioPython/Total_pocket_analysis")

from Pocket_List_Builder_score_key import count_annotations, build_pocket_list, pocket_list_score
import pprint

# Αρχείο GenBank
gb_file = "ADSS2.gb"
# Patterns για κάθε πηγή
patterns = ["Fpockt#", "DGsite#", "ASpock#"]

# Πρώτα μετράμε πόσα έχει κάθε κατηγορία
counts = count_annotations(gb_file, patterns)
print("Counts:", counts)

# Δημιουργούμε το pocket_list με τα residues
pocket_list = build_pocket_list(gb_file, patterns)

# Εφαρμόζουμε τα score1/score2 ανά πηγή
# - Fpocket: score1 = 'Drug Score', score2 = 'Pocket Score'
# - DGSite: score1 = 'drugScore', score2 = 'depth'
# - ASpock: score1 = 'e', score2 = 'v*buriedness^2/rg'

# Βήμα 1: προσθήκη score1/score2 για Fpocket
pocket_list = pocket_list_score(gb_file, pocket_list, ["Fpockt#"], "Drug Score", "Pocket Score")
# Βήμα 2: προσθήκη score1/score2 για DGSite
pocket_list = pocket_list_score(gb_file, pocket_list, ["DGsite#"], "drugScore", "depth")
# Βήμα 3: προσθήκη score1/score2 για ASpock απευθείας από qualifiers 'e' και 'v*buriedness^2/rg'
# (χρησιμοποιούμε την ίδια pocket_list_score)
pocket_list = pocket_list_score(gb_file, pocket_list, ["ASpock#"], "e", "v*buriedness^2/rg")

# Εμφανίζουμε τελικά το pocket_list με scores
pp = pprint.PrettyPrinter(width=170, compact=True)
pp.pprint(pocket_list)
