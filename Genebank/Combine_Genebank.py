import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def combine_genbank_files(file1, file2, output_file):
    try:
        # Read the sequences from the input GenBank files
        records1 = list(SeqIO.parse(file1, "genbank"))
        records2 = list(SeqIO.parse(file2, "genbank"))

        if not records1 or not records2:
            raise ValueError("One or both input files are empty or not in GenBank format.")

        # Combine sequences and their annotations
        combined_features = []
        for record in records1 + records2:
            combined_features.extend(record.features)

        # Use the sequence only from file1
        combined_sequence = str(records1[0].seq)

        # Get fields from file1
        locus = records1[0].id if records1[0].id else "Combined"
        accession = records1[0].annotations.get("accessions", ["Unknown"])[0]
        version = records1[0].annotations.get("sequence_version", "1")

        # Merge DEFINITION fields
        definition1 = records1[0].description if records1[0].description else ""
        definition2 = records2[0].description if records2[0].description else ""
        merged_definition = f"{definition1} {definition2}".strip()

        # Create a new SeqRecord object with combined data
        combined_record = SeqRecord(
            Seq(combined_sequence),
            id=locus,
            name=locus,
            description=merged_definition,
            annotations={**records1[0].annotations, **records2[0].annotations},
        )

        # Update specific annotations
        combined_record.annotations["accessions"] = [accession]
        combined_record.annotations["sequence_version"] = version

        # Add all features to the new record
        combined_record.features = combined_features

        # Write the combined record to an output GenBank file
        SeqIO.write(combined_record, output_file, "genbank")

    except Exception as e:
        print(f"Error during processing: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <file1> <file2> <output_file>")
        sys.exit(1)

    # Get file paths from command line arguments
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    output_file = sys.argv[3]

    print("Starting to combine GenBank files...")
    combine_genbank_files(file1, file2, output_file)
    print(f"Combined GenBank file saved to {output_file}")
