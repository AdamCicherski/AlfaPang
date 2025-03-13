from Bio import SeqIO
import sys

def filter_fasta(fasta_file, names_file, output_file):
    # Read sequence names from file
    with open(names_file, "r") as f:
        sequence_names = set(line.strip() for line in f)
    
    # Filter and write sequences to new FASTA file
    with open(output_file, "w") as out_f:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in sequence_names:
                SeqIO.write(record, out_f, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input.fasta> <names.txt> <output.fasta>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    names_file = sys.argv[2]
    output_file = sys.argv[3]
    
    filter_fasta(fasta_file, names_file, output_file)

