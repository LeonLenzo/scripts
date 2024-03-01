#!/bin/bash

working_dir="."

for reads in $working_dir/*.fastq.gz; do

    # Extract sample IDs without extensions
    sample_ID=$(basename -s _1.fastq.gz "$reads")

    #Trim Raw Reads
	bbduk.sh \
	in1="${sample_ID}_1.fastq.gz" \
	in2="${sample_ID}_2.fastq.gz" \
	out1="${sample_ID}_1.trimmed.fastq.gz" \
	out2="${sample_ID}_2.trimmed.fastq.gz" \
	ref="/onedrive/Omics/Metagenomics/ll_adapt.fa" \
	
done
