import sys
from Bio import SeqIO
from Bio.Seq import Seq

def replace_sequence_in_genebank(genebank_file, fasta_file, output_file):
    # Read the FASTA sequence
    fasta_seq = None
    with open(fasta_file, "r") as fasta_handle:
        fasta_seq = str(next(SeqIO.parse(fasta_handle, "fasta")).seq)
    
    # Parse the GeneBank file
    genbank_record = None
    with open(genebank_file, "r") as genebank_handle:
        genbank_record = SeqIO.read(genebank_handle, "genbank")
    
    # Replace the sequence in the GeneBank record
    genbank_record.seq = Seq(fasta_seq)  # Ensure sequence is a Seq object
    
    # Write the updated GeneBank file
    with open(output_file, "w") as output_handle:
        SeqIO.write(genbank_record, output_handle, "genbank")
    print(f"Updated GeneBank file written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python replace_sequence.py <gene_bank_file> <fasta_file> <output_file>")
        sys.exit(1)
    
    genebank_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3]

    replace_sequence_in_genebank(genebank_file, fasta_file, output_file)
