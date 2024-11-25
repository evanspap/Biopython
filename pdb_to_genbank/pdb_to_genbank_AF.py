import sys
import os
from Bio.PDB import PDBParser, PPBuilder
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio import SeqIO


def pdb_to_genbank_with_bfactor_categories(input_pdb):
    """
    Converts a PDB file into a GenBank file with annotations based on B-factor categories.

    Args:
        input_pdb (str): Path to the input PDB file.
    """
    # Generate output GenBank filename and record name
    basename = os.path.splitext(os.path.basename(input_pdb))[0]
    output_genbank = f"{basename}.gb"

    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(basename, input_pdb)

    # Extract the amino acid sequence
    ppb = PPBuilder()
    seq = ""
    for pp in ppb.build_peptides(structure):
        seq += str(pp.get_sequence())  # Combine sequences from all polypeptides

    # Define B-factor categories
    bfactor_ranges = {
        "b_50": lambda b: b < 50,
        "50_b_70": lambda b: 50 <= b < 70,
        "70_b_90": lambda b: 70 <= b < 90,
        "90_b": lambda b: b >= 90,
    }

    # Initialize annotations for each category
    bfactor_annotations = {key: [] for key in bfactor_ranges}

    # Extract residues and categorize by B-factor
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:  # Use the alpha carbon atom as a reference for residue
                    atom = residue["CA"]
                    bfactor = atom.get_bfactor()

                    for category, condition in bfactor_ranges.items():
                        if condition(bfactor):
                            # Use residue sequence ID
                            res_id = residue.id[1]
                            bfactor_annotations[category].append(res_id)

    # Create SeqRecord for GenBank file
    record = SeqRecord(
        Seq(seq),  # Use the actual sequence extracted from the PDB
        id=basename[:10],  # GenBank `id` must be short and alphanumeric (max 10 chars)
        name=basename[:10],  # GenBank `name` must be short and alphanumeric (max 10 chars)
        description=f"B-factor annotations for {basename}",
        annotations={"molecule_type": "protein"},  # Required field for GenBank
    )

    # Add annotations to SeqRecord
    for category, residues in bfactor_annotations.items():
        if residues:
            # Join contiguous residues as a single feature
            locations = [FeatureLocation(start - 1, start, strand=1) for start in residues]
            record.features.append(
                SeqFeature(
                    location=sum(locations),
                    type=category,  # Use the revised naming convention for feature type
                )
            )

    # Write GenBank file
    SeqIO.write(record, output_genbank, "genbank")
    print(f"GenBank file saved as {output_genbank}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_pdb_file>")
        sys.exit(1)

    input_pdb = sys.argv[1]
    pdb_to_genbank_with_bfactor_categories(input_pdb)
