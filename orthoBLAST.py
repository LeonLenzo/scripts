import argparse
import subprocess
from Bio import SeqIO

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("genes1", help="Path to the FASTA file containing known genes from genome 1.")
parser.add_argument("genes2", help="Path to the FASTA file containing known genes from genome 2.")
parser.add_argument("-e", "--evalue", type=float, default=1e-6, help="E-value threshold for considering matches as significant. Default is 1e-6.")
args = parser.parse_args()

# Run BLAST from genes1 to genes2
def run_blast(query_file, subject_file, evalue):
    cmd = f"blastn -query {query_file} -subject {subject_file} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'"
    blast_output = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return blast_output.stdout.splitlines()

blast_output = run_blast(args.genes1, args.genes2, args.evalue)

# Filter BLAST output to include results with alignment length > 250
filtered_output = [line for line in blast_output if int(line.split('\t')[3]) > 250]

# Sort filtered output by % Identity in descending order
sorted_output = sorted(filtered_output, key=lambda x: float(x.split('\t')[2]), reverse=True)

# Write column headings to sorted BLAST output
column_headings = "Query ID\tSubject ID\t% Identity\tAlignment Length\tMismatches\tGap Opens\tQuery Start\tQuery End\tSubject Start\tSubject End\tE-value\tBit Score\n"
sorted_output.insert(0, column_headings)

# Count the number of blast hits in the final output
num_blast_hits = len(sorted_output) - 1  # Exclude the header line

print(f"Potential Targets: {num_blast_hits}")

# Write sorted BLAST output to a file with the name of the input files
output_file = f"{args.genes1}_vs_{args.genes2}_orthologs.txt"
with open(output_file, "w") as file:
    file.write("\n".join(sorted_output))

print(f"Output written to {output_file}")

# Get the sequences of queries and targets
queries = set()
targets = set()

for line in sorted_output[1:]:
    query_id, subject_id = line.split('\t')[:2]
    queries.add(query_id)
    targets.add(subject_id)

# Write sequences of queries to a FASTA file with the name of the input file
query_sequences = SeqIO.to_dict(SeqIO.parse(args.genes1, "fasta"))
query_output_file = f"{args.genes1}.targets.fasta"
with open(query_output_file, "w") as f:
    for query_id in queries:
        SeqIO.write(query_sequences[query_id], f, "fasta")

print(f"Filtered queries written to {query_output_file}")

# Write sequences of targets to a FASTA file with the name of the input file
target_sequences = SeqIO.to_dict(SeqIO.parse(args.genes2, "fasta"))
target_output_file = f"{args.genes2}.targets.fasta"
with open(target_output_file, "w") as f:
    for target_id in targets:
        SeqIO.write(target_sequences[target_id], f, "fasta")

print(f"Filtered targets written to {target_output_file}")
