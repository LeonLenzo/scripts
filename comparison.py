import os
from collections import defaultdict
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.motifs import create
import subprocess

def regroup_fasta_files(working_directory):
    """Regroup sequences in FASTA files by gene instead of isolate."""
    gene_dict = defaultdict(list)

    # Read all consensus files in the working directory
    for file_name in os.listdir(working_directory):
        if file_name.endswith('.consensus'):
            file_path = os.path.join(working_directory, file_name)
            for record in SeqIO.parse(file_path, 'fasta'):
                header_parts = record.id.split('.')
                if len(header_parts) == 1:
                    gene = header_parts[0]
                else:
                    gene, isolate = header_parts
                gene_dict[gene].append(record)

    # Write the regrouped sequences to new FASTA files in the working directory
    for gene, sequences in gene_dict.items():
        output_file_path = os.path.join(working_directory, f'{gene}.fasta')
        SeqIO.write(sequences, output_file_path, 'fasta')

def align_sequences(input_file, output_file):
    """Perform multiple sequence alignment using MUSCLE."""
    with open(output_file, 'w') as out_handle:
        subprocess.run(['muscle', '-in', input_file, '-out', output_file, '-quiet'], stdout=out_handle, check=True)

def generate_consensus(alignment_file):
    """Generate a consensus sequence from the alignment."""
    alignment = AlignIO.read(alignment_file, 'fasta')
    motif = create(alignment)
    consensus = motif.degenerate_consensus
    return str(consensus)

def align_and_generate_consensus(working_directory):
    """Align sequences in FASTA files and generate consensus sequences."""
    consensus_sequences = []

    for file_name in os.listdir(working_directory):
        if file_name.endswith('.fasta'):
            input_file = os.path.join(working_directory, file_name)
            alignment_file = input_file.replace('.fasta', '_aligned.fasta')
            align_sequences(input_file, alignment_file)
            consensus = generate_consensus(alignment_file)
            gene_name = os.path.splitext(file_name)[0]
            consensus_sequences.append((gene_name, consensus))
            os.remove(alignment_file)  # Clean up intermediate file

    # Write all consensus sequences to a single file
    output_file_path = os.path.join(working_directory, 'all_consensus_sequences.fasta')
    with open(output_file_path, 'w') as output_file:
        for gene_name, consensus in consensus_sequences:
            output_file.write(f'>{gene_name}\n{consensus}\n')

def main():
    working_directory = os.getcwd()
    
    print('Regrouping FASTA files by gene...')
    regroup_fasta_files(working_directory)

    print('Aligning sequences and generating consensus sequences...')
    align_and_generate_consensus(working_directory)

    print('Done!')

if __name__ == "__main__":
    main()
