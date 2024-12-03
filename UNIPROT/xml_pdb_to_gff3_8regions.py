import xml.etree.ElementTree as ET
import os
import argparse
from collections import defaultdict

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
        sequence_length = len(sequence)

        # Extract residue ranges for X-ray, NMR, and EM
        xray_ranges = []
        nmr_ranges = []
        em_ranges = []

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
                            if method == "X-ray":
                                xray_ranges.append((start, end))
                            elif method == "NMR":
                                nmr_ranges.append((start, end))
                            elif method == "EM":
                                em_ranges.append((start, end))

        # Assign each residue to one of the 8 regions
        region_labels = ["No PDB structure", "X-ray", "NMR", "EM", "X-ray and NMR", "X-ray and EM", "NMR and EM", "X-ray and NMR and EM"]
        residue_annotations = [0] * sequence_length

        def mark_regions(ranges, value):
            for start, end in ranges:
                for i in range(start - 1, end):
                    residue_annotations[i] |= value

        mark_regions(xray_ranges, 1)
        mark_regions(nmr_ranges, 2)
        mark_regions(em_ranges, 4)

        # Consolidate annotations into regions
        annotated_regions = defaultdict(list)
        current_region = None
        region_start = None

        for i, annotation in enumerate(residue_annotations):
            if annotation != current_region:
                if current_region is not None:
                    region_label = region_labels[current_region]
                    annotated_regions[region_label].append((region_start + 1, i))
                current_region = annotation
                region_start = i

        # Finalize the last region
        if current_region is not None:
            region_label = region_labels[current_region]
            annotated_regions[region_label].append((region_start + 1, sequence_length))

        proteins.append({"seqid": seqid, "sequence": sequence, "regions": annotated_regions})

    return proteins

# Generate GFF3 files
def generate_gff3(proteins, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for protein in proteins:
        seqid = protein["seqid"]
        sequence = protein["sequence"]
        regions = protein["regions"]

        output_path = os.path.join(output_dir, f"{seqid}.gff")
        with open(output_path, "w") as gff_file:
            # Write GFF3 header
            gff_file.write("##gff-version 3\n")
            gff_file.write(f"##sequence-region {seqid} 1 {len(sequence)}\n")

            # Write annotated regions
            for region_label, ranges in regions.items():
                for start, end in ranges:
                    gff_file.write(
                        f"{seqid}\tUniProt\t{region_label}\t{start}\t{end}\t.\t.\t.\t"
                        f"ID={region_label.replace(' ', '_')}_{start}_{end}\n"
                    )

            # Append sequence in FASTA format
            gff_file.write("##FASTA\n")
            gff_file.write(f">{seqid}\n")
            for i in range(0, len(sequence), 100):
                gff_file.write(sequence[i:i + 100] + "\n")

# Main execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse UNIPROT XML and generate GFF3 annotation files.")
    parser.add_argument("input_file", help="Path to the UNIPROT XML file")
    parser.add_argument("--output_dir", default="gff3_files", help="Directory to save GFF3 files")

    args = parser.parse_args()

    proteins = parse_uniprot_xml(args.input_file)
    generate_gff3(proteins, args.output_dir)

    print(f"GFF3 files have been generated in the directory: {args.output_dir}")
