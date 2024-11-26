import os
import sys
from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain


def filter_atoms_with_neighbors(input_pdb, output_pdb, bfactor_threshold, neighbor_residues):
    """
    Filters atoms based on a B-factor threshold and includes neighboring residues.

    Args:
        input_pdb (str): Path to the input PDB file.
        output_pdb (str): Path to the output PDB file containing filtered atoms.
        bfactor_threshold (float): The B-factor threshold.
        neighbor_residues (int): Number of residues to include around selected residues.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", input_pdb)

    # Create a new structure to hold filtered atoms
    filtered_structure = Structure.Structure("filtered_structure")

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


def batch_filter_pdb(input_directory, output_directory, bfactor_threshold, neighbor_residues):
    """
    Processes all PDB files in a directory, filtering them based on B-factor and neighbors,
    and saves the results in an output directory.

    Args:
        input_directory (str): Path to the directory containing PDB files.
        output_directory (str): Path to the directory where filtered PDB files will be saved.
        bfactor_threshold (float): The B-factor threshold.
        neighbor_residues (int): Number of residues to include around selected residues.
    """
    # Ensure the output directory exists
    os.makedirs(output_directory, exist_ok=True)

    # Find all PDB files in the input directory
    pdb_files = [f for f in os.listdir(input_directory) if f.endswith(".pdb")]
    if not pdb_files:
        print(f"No PDB files found in directory: {input_directory}")
        sys.exit(1)

    for pdb_file in pdb_files:
        input_pdb = os.path.join(input_directory, pdb_file)
        output_pdb = os.path.join(output_directory, pdb_file)

        print(f"Processing {input_pdb}...")
        filter_atoms_with_neighbors(input_pdb, output_pdb, bfactor_threshold, neighbor_residues)
        print(f"Filtered PDB saved to {output_pdb}")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python batch_filter_pdb.py <input_directory> <output_directory> <bfactor_threshold> <neighbor_residues>")
        sys.exit(1)

    # Get arguments from the command line
    input_directory = sys.argv[1]
    output_directory = sys.argv[2]
    bfactor_threshold = float(sys.argv[3])
    neighbor_residues = int(sys.argv[4])

    # Run the batch processing
    batch_filter_pdb(input_directory, output_directory, bfactor_threshold, neighbor_residues)
