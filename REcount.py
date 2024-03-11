#!/usr/bin/env python3

import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

def count_restriction_sites(fasta_file, enzyme_data):
    # Initialize a dictionary to store counts for each enzyme
    enzyme_counts = {enzyme: 0 for enzyme in enzyme_data.iloc[:, 0]}

    # Step 1: Read the FASTA file
    records = SeqIO.parse(fasta_file, "fasta")

    # Step 2: Search for recognition sites
    for record in records:
        sequence = str(record.seq)

        # Step 3: Count occurrences for each enzyme sequence
        for _, row in enzyme_data.iterrows():
            enzyme_sequence = str(row.iloc[1])
            occurrences = str(sequence).count(enzyme_sequence)
            enzyme_counts[row.iloc[0]] += occurrences

    return enzyme_counts

def read_enzyme_data(file_path):
    # Read enzyme data from a CSV file (no headers)
    enzyme_data = pd.read_csv(file_path, header=None)
    return enzyme_data

def process_directory(directory_path, enzyme_data):
    # Create a DataFrame to store results
    results_df = pd.DataFrame(index=enzyme_data.iloc[:, 0])

    # Process each FASTA file in the directory
    for file_name in os.listdir(directory_path):
        if file_name.endswith(".fasta"):
            fasta_file_path = os.path.join(directory_path, file_name)

            enzyme_counts = count_restriction_sites(fasta_file_path, enzyme_data)

            # Add results to DataFrame
            results_df[file_name[:-6]] = pd.Series(enzyme_counts)

    return results_df

def main():
    # Prompt the user for input directory and enzyme file
    input_directory = input("Enter the path to the input directory: ")
    enzyme_file_path = input("Enter the path to the enzyme data file (.csv): ")

    # Read enzyme data from the file
    enzyme_data = read_enzyme_data(enzyme_file_path)

    # Process each FASTA file in the directory and compile results into a table
    results_df = process_directory(input_directory, enzyme_data)

    # Print the results DataFrame
    print("\nResults Table:")
    print(results_df)

    # Save the results DataFrame to a CSV file
    results_df.to_csv("REcount_table.csv")

if __name__ == "__main__":
    main()

