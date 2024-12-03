import xml.etree.ElementTree as ET
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio import SeqIO
import os
import argparse

# Parse the UNIPROT XML file
def parse_uniprot_xml(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()

    namespace = {'ns': "http://uniprot.org/uniprot"}
    proteins = []

    for entry in root.findall("ns:entry", namespace):
        # Get LOCUS name from the <name> field
        locus_name = entry.find("ns:name", namespace)
        locus_name = locus_name.text if locus_name is not None else "Unknown"

        # Get the protein full name
        protein_name = entry.find("ns:protein/ns:recommendedName/ns:fullName", namespace)
        protein_name = protein_name.text if protein_name is not None else "Unknown"

        # Get sequence and its length
        sequence_element = entry.find("ns:sequence", namespace)
        sequence = sequence_element.text.strip() if sequence_element is not None else ""
        length = int(sequence_element.attrib.get("length", "0"))

        # Extract regions from <dbReference>
        features = []
        for db_ref in entry.findall("ns:dbReference", namespace):
            method = db_ref.attrib.get("type")
            if method == "PDB":
                for prop in db_ref.findall("ns:property", namespace):
                    if prop.attrib.get("type") == "method":
                        structure_method = prop.attrib.get("value")
                    elif prop.attrib.get("type") == "chains":
                        chains = prop.attrib.get("value")
                        for chain in chains.split(","):
                            parts = chain.split("=")
                            if len(parts) == 2:
                                residues = parts[1]
                                start, end = map(int, residues.split("-"))
                                features.append({
                                    "method": structure_method,
                                    "start": start,
                                    "end": end,
                                    "evidence": f"PBD:{db_ref.attrib.get('id')}",
                                    "label": db_ref.attrib.get("id")
                                })

        proteins.append({
            "locus_name": locus_name,
            "protein_name": protein_name,
            "sequence": sequence,
            "length": length,
            "features": features
        })

    return proteins

# Generate GenBank files
def generate_genbank_files(proteins, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for protein in proteins:
        # Create a GenBank record
        record = SeqRecord(
            Seq(protein["sequence"]),
            id=protein["locus_name"],
            name=protein["locus_name"],
            description=protein["protein_name"]
        )
        record.annotations["molecule_type"] = "protein"

        # Add features for regions
        for feature in protein["features"]:
            genbank_feature = SeqFeature(
                FeatureLocation(feature["start"] - 1, feature["end"]),
                type=feature["method"],
                qualifiers={
                    "evidence": feature["evidence"],
                    "label": feature["label"]
                }
            )
            record.features.append(genbank_feature)

        # Set annotations
        record.annotations["length"] = f"{protein['length']}aa"

        # Save to GenBank file
        output_path = os.path.join(output_dir, f"{protein['locus_name']}.gb")
        with open(output_path, "w") as output_handle:
            SeqIO.write(record, output_handle, "genbank")

# Main execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse UNIPROT XML and generate GenBank files.")
    parser.add_argument("input_file", help="Path to the UNIPROT XML file")
    parser.add_argument("--output_dir", default="genbank_files", help="Directory to save GenBank files")

    args = parser.parse_args()

    proteins = parse_uniprot_xml(args.input_file)
    generate_genbank_files(proteins, args.output_dir)

    print(f"GenBank files have been generated in the directory: {args.output_dir}")
