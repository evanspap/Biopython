import xml.etree.ElementTree as ET
import csv
import argparse

# Parse the UNIPROT XML file
def parse_uniprot_xml(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()

    namespace = {'ns': "http://uniprot.org/uniprot"}
    proteins = []

    for entry in root.findall("ns:entry", namespace):
        # Get LOCUS name (for protein_name in CSV)
        protein_name = entry.find("ns:name", namespace)
        protein_name = protein_name.text if protein_name is not None else "Unknown"

        # Count the number of X-ray, NMR, and EM structures
        xray_count = 0
        nmr_count = 0
        em_count = 0

        for db_ref in entry.findall("ns:dbReference", namespace):
            if db_ref.attrib.get("type") == "PDB":
                for prop in db_ref.findall("ns:property", namespace):
                    if prop.attrib.get("type") == "method":
                        method = prop.attrib.get("value")
                        if method == "X-ray":
                            xray_count += 1
                        elif method == "NMR":
                            nmr_count += 1
                        elif method == "EM":
                            em_count += 1

        proteins.append({
            "protein_name": protein_name,
            "Nr_X-ray": xray_count,
            "Nr_NMR": nmr_count,
            "Nr_EM": em_count
        })

    return proteins

# Generate CSV summary
def generate_csv(proteins, output_csv):
    with open(output_csv, mode='w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["protein_name", "Nr_X-ray", "Nr_NMR", "Nr_EM"])
        writer.writeheader()
        writer.writerows(proteins)

# Main execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse UNIPROT XML and generate CSV summary.")
    parser.add_argument("input_file", help="Path to the UNIPROT XML file")
    parser.add_argument("--output_csv", default="protein_summary.csv", help="Output CSV file")

    args = parser.parse_args()

    proteins = parse_uniprot_xml(args.input_file)
    generate_csv(proteins, args.output_csv)

    print(f"CSV summary has been generated: {args.output_csv}")
