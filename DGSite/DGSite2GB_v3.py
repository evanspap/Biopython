import os
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.PDB import PDBParser

def parse_desc_file(desc_file):
    """Parse the *_desc.txt file to extract selected feature data."""
    df = pd.read_csv(desc_file, sep='\t')
    selected_features = [
        "volume", "hull", "surface", "lid", "depth", "surf/vol", "lid/hull", "ellVol", "ell c/a", "ell b/a", 
        "siteAtms", "accept", "donor", "aromat", "hydrophobicity", "metal", "Cs", "Ns", "Os", "Ss", "Xs", 
        "negAA", "posAA", "polarAA", "apolarAA", "sumAA", "ALA", "ARG", "ASN", "ASP", "CSO", "CYS", "GLN", 
        "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", 
        "A", "C", "G", "U", "I", "N", "DA", "DC", "DG", "DT", "DN", "UNK", "simpleScore", "drugScore"
    ]
    df = df[["name"] + selected_features]
    return df.set_index("name").to_dict(orient='index')

def extract_residue_numbers_from_pdb(pdb_file):
    """Extract residue numbers from a PDB file using Bio.PDB."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    
    residue_numbers = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == " ":  # Ignore heteroatoms
                    residue_numbers.append(residue.id[1])
    
    return sorted(set(residue_numbers))

def create_genbank_record(pockets, pdb_files):
    """Generate a GenBank record incorporating pocket feature data and residue structures."""
    sequence = "NNNNNNNNNN"  # Placeholder sequence as PDB files don't contain sequence data
    record = SeqRecord(Seq(sequence), id="Generated_GB", name="Generated_GB",
                       description="Generated GenBank file with pocket structures.",
                       annotations={"molecule_type": "protein"})
    
    for index, (pocket_id, pocket_data) in enumerate(pockets.items(), start=1):
        matched_pdbs = [pdb for pdb in pdb_files if pocket_id in os.path.basename(pdb)]
        residues = []
        for pdb_file in matched_pdbs:
            residues.extend(extract_residue_numbers_from_pdb(pdb_file))
        residues = sorted(set(residues))
        
        if residues:
            qualifiers = {k: v for k, v in pocket_data.items() if k in pocket_data}
            location = CompoundLocation([FeatureLocation(pos - 1, pos) for pos in residues])
            feature = SeqFeature(location, type=f"DGsite{index:02}", qualifiers=qualifiers)
            record.features.append(feature)
    
    return record

def main():
    if len(sys.argv) < 4:
        print("Usage: python DGSite2GB_v3.py <desc_file> <pdb_folder> <output.gb>")
        sys.exit(1)
    
    desc_file = sys.argv[1]
    pdb_folder = sys.argv[2]
    output_file = sys.argv[3]
    
    print(f"Running: python {' '.join(sys.argv)}")
    
    if not os.path.exists(desc_file):
        print(f"Error: Description file {desc_file} not found.")
        sys.exit(1)
    
    if not os.path.exists(pdb_folder):
        print(f"Error: PDB folder {pdb_folder} not found.")
        sys.exit(1)
    
    pdb_files = [os.path.join(pdb_folder, f) for f in os.listdir(pdb_folder) if f.endswith(".pdb")]
    if not pdb_files:
        print("Warning: No PDB files found in the specified folder.")
    
    pockets = parse_desc_file(desc_file)
    record = create_genbank_record(pockets, pdb_files)
    
    SeqIO.write(record, output_file, "genbank")
    print(f"GenBank file saved as {output_file}")

if __name__ == "__main__":
    main()
