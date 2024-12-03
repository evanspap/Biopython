import xml.etree.ElementTree as ET
import os
import argparse

# Parse the UNIPROT XML file
def parse_uniprot_xml(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()

    namespace = {'ns': "http://uniprot.org/uniprot"}
    proteins = []

    for entry in root.findall("ns:entry", namespace):
        # Get the <name> field for seqid
        seqid = entry.find("ns:name", namespace)
        seqid = seqid.text if seqid is not None else "Unknown"

        # Get the protein sequence
        sequence_element = entry.find("ns:sequence", namespace)
        sequence = sequence_element.text.strip() if sequence_element is not None else ""

        # Extract features
        features = []
        for db_ref in entry.findall("ns:dbReference", namespace):
            if db_ref.attrib.get("type") == "PDB":
                method = None
                chains = None
                for prop in db_ref.findall("ns:property", namespace):
                    if prop.attrib.get("type") == "method":
                        method = prop.attrib.get("value")
                    elif prop.attrib.get("type") == "chains":
                        chains = prop.attrib.get("value")
                if method and chains:
                    for chain in chains.split(","):
                        parts = chain.split("=")
                        if len(parts) == 2:
                            residues = parts[1]
                            start, end = map(int, residues.split("-"))
                            features.append({
                                "seqid": seqid,
                                "source": "UniProt",
                                "type": method,
                                "start": start,
                                "end": end,
                                "score": ".",
                                "strand": ".",
                                "phase": ".",
                                "attributes": f"ID={method}_{start}_{end};Evidence=PBD:{db_ref.attrib.get('id')}"
                            })

        proteins.append({"seqid": seqid, "features": features, "sequence": sequence})

    return proteins

# Generate GFF3 files
def generate_gff3(proteins, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for protein in proteins:
        seqid = protein["seqid"]
        features = protein["features"]
        sequence = protein["sequence"]

        output_path = os.path.join(output_dir, f"{seqid}.gff")
        with open(output_path, "w") as gff_file:
            # Write GFF3 header
            gff_file.write("##gff-version 3\n")
            if features:
                gff_file.write(f"##sequence-region {seqid} 1 {max([f['end'] for f in features])}\n")

            # Write features
            for feature in features:
                gff_file.write(
                    f"{feature['seqid']}\t{feature['source']}\t{feature['type']}\t"
                    f"{feature['start']}\t{feature['end']}\t{feature['score']}\t"
                    f"{feature['strand']}\t{feature['phase']}\t{feature['attributes']}\n"
                )

            # Append sequence in FASTA format
            if sequence:
                gff_file.write("##FASTA\n")
                gff_file.write(f">{seqid}\n")
                # Split sequence into 60-character lines
                for i in range(0, len(sequence), 60):
                    gff_file.write(sequence[i:i + 60] + "\n")

# Main execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse UNIPROT XML and generate GFF3 annotation files.")
    parser.add_argument("input_file", help="Path to the UNIPROT XML file")
    parser.add_argument("--output_dir", default="gff3_files", help="Directory to save GFF3 files")

    args = parser.parse_args()

    proteins = parse_uniprot_xml(args.input_file)
    generate_gff3(proteins, args.output_dir)

    print(f"GFF3 files have been generated in the directory: {args.output_dir}")
