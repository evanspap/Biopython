from Bio.PDB import PDBParser, Superimposer



def select_residues_by_distance(target_pdb_file, reference_pdb_file, distance_cutoff):

    parser = PDBParser()

    target_structure = parser.get_structure("target", target_pdb_file)

    reference_structure = parser.get_structure("reference", reference_pdb_file)



    # Select atoms from reference structure

    reference_selection = [atom in reference_structure ]
    print (reference_selection)



    selected_residues = []

    for atom in target_structure:

        for selection_atom in reference_selection:


            if atom.get_distance(selection_atom) <= distance_cutoff:

                        selected_residues.append(residue)

                        break  # Move to next residue if a match is found



    return selected_residues



# Example usage

target_pdb = "target.pdb"

reference_pdb = "reference.pdb"

distance_threshold = 3.0  # Angstroms



selected_residues = select_residues_by_distance(target_pdb, reference_pdb, distance_threshold)



print("Selected residues:", [residue.get_id() for residue in selected_residues]) 
