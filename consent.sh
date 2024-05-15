#!/bin/bash

# Define arguments
INDEX="$1" # Contains the path to the directory containing your index
REFERENCE="$2" # Contains the path to your reference .fasta synthetic genome

# Iterate each paired .fastq.gz files in the directory
for file in *_1.trimmed.fastq.gz; do

    # Extract sample name
    sample=$(basename "$file" _1.trimmed.fastq.gz)
    
    # Define input files
    read1="$sample"_1.trimmed.fastq.gz
    read2="$sample"_2.trimmed.fastq.gz

    # Map reads using HISAT2
    echo "Mapping reads for: $sample"
    hisat2 -p 16 -t -x $INDEX -1 $read1 -2 $read2 -S "$sample.sam"

    # Convert .sam to .bam
    echo "Converting $sample.sam to $sample.bam"
    samtools view -bS "$sample.sam" > "$sample.bam"
    rm -r "$sample.sam"

    # Sort and index .bam
    echo "Sorting: $sample.bam"
    samtools sort "$sample".bam -o "$sample".sorted.bam

    echo "Indexing: $sample.bam"
    samtools index "$sample".sorted.bam
    rm -r "$sample.bam"

    # Generate VCF file
    echo "Generating VCF for $sample.sorted.bam"
    bcftools mpileup -Ou -f "$REFERENCE" "$sample".sorted.bam | bcftools call -mv -Ov -o "$sample".vcf
    rm -r "$sample.sorted.bam"

    # Compress and index VCF file
    echo "Compressing $sample.vcf"
    bgzip "$sample".vcf

    echo "Indexing $sample.vcf"
    tabix -p vcf "$sample".vcf.gz

    # Generate consensus sequence
    echo "Generating consensus for $sample"
    bcftools consensus -f "$REFERENCE" "$sample".vcf.gz > "$sample".consensus.fasta
    rm -r *vcf*

    # Append sample name to sequence names
    echo "Appending $sample to sequence names"
    awk -v sample="$sample" '/^>/{print $0"_"sample;next}{print}' "$sample".consensus.fasta > "$sample".consensus_temp.fasta
    mv "$sample".consensus_temp.fasta "$sample".consensus.fasta

done

echo "Mapping and consensus sequence generation complete."
