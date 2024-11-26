import sys
from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain, Residue


def filter_atoms_with_neighbors(input_pdb, output_pdb, bfactor_threshold, neighbor_residues):
    """
    Filters atoms based on a B-factor threshold and includes neighboring residues.

    Args:
        input_pdb (str): Path to the input PDB file.
        output_pdb (str): Path to the output PDB file containing filtered atoms.
        bfactor_threshold (float): The B-factor threshold.
        neighbor_residues (int): Number of residues to include around selected residues.
    """
    # Parse the input PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", input_pdb)

    # Create a new structure to hold filtered atoms
    filtered_structure = Structure.Structure("filtered_structure")

    # Loop through each atom and filter by B-factor
    for model in structure:
        filtered_model = Model.Model(model.id)
        for chain in model:
            filtered_chain = Chain.Chain(chain.id)
            residues = list(chain)  # Get all residues in the chain
            selected_residue_indices = []

            # Find residues with at least one atom meeting the B-factor threshold
            for i, residue in enumerate(residues):
                if any(atom.get_bfactor() > bfactor_threshold for atom in residue):
                    selected_residue_indices.append(i)

            # Add neighboring residues to the selection
            extended_indices = set()
            for index in selected_residue_indices:
                start = max(0, index - neighbor_residues)
                end = min(len(residues), index + neighbor_residues + 1)
                extended_indices.update(range(start, end))

            # Add filtered residues to the chain
            for index in sorted(extended_indices):
                filtered_chain.add(residues[index])

            if len(filtered_chain):
                filtered_model.add(filtered_chain)
        if len(filtered_model):
            filtered_structure.add(filtered_model)

    # Save the filtered structure to the output PDB file
    io = PDBIO()
    io.set_structure(filtered_structure)
    io.save(output_pdb)
    print(f"Filtered PDB file saved as {output_pdb}")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <input_pdb> <output_pdb> <bfactor_threshold> <neighbor_residues>")
        sys.exit(1)

    # Get arguments from the command line
    input_pdb = sys.argv[1]
    output_pdb = sys.argv[2]
    bfactor_threshold = float(sys.argv[3])
    neighbor_residues = int(sys.argv[4])

    # Check if input and output filenames are the same
    if input_pdb == output_pdb:
        print("Error: The output PDB file cannot have the same name as the input PDB file.")
        sys.exit(1)

    # Run the function
    filter_atoms_with_neighbors(input_pdb, output_pdb, bfactor_threshold, neighbor_residues)
