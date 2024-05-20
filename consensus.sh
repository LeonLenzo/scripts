#!/bin/bash

# Define arguments
INDEX="$1" # Contains the path to the directory containing your index
REFERENCE="$2" # Contains the path to your reference .fasta synthetic genome

# Function to display elapsed time
function display_elapsed_time {
    local current_time=$(date +%s)
    local elapsed_seconds=$((current_time - start_time))
    local hours=$((elapsed_seconds / 3600))
    local minutes=$((elapsed_seconds % 3600 / 60))
    local seconds=$((elapsed_seconds % 60))
    printf "Elapsed time for $sample: %02d:%02d:%02d\n" $hours $minutes $seconds
}

# Iterate each paired .fastq.gz files in the directory
for file in *_1.trimmed.fastq.gz; do

    # Extract sample name
    sample=$(basename "$file" _1.trimmed.fastq.gz)
    
    # Define input files
    read1="$sample"_1.trimmed.fastq.gz
    read2="$sample"_2.trimmed.fastq.gz

    # Start timer for sample processing
    start_time=$(date +%s)

    # Map reads using HISAT2
    echo "Mapping reads for: $sample"
    hisat2 -p 16 -t -x $INDEX -1 $read1 -2 $read2 -S "$sample.sam"
    display_elapsed_time

    # Convert .sam to .bam
    echo "Converting $sample.sam to $sample.bam"
    samtools view -@ 8 -bS "$sample.sam" > "$sample.bam"  # Use 8 threads for compression
    rm -r "$sample.sam"
    display_elapsed_time

    # Sort and index .bam
    echo "Sorting: $sample.bam"
    samtools sort "$sample".bam -o "$sample".sorted.bam
    echo "Indexing: $sample.bam"
    samtools index "$sample".sorted.bam
    rm -r "$sample.bam"
    display_elapsed_time

    # Generate VCF file
    echo "Generating VCF for $sample.sorted.bam"
    bcftools mpileup -Ou -f "$REFERENCE" "$sample".sorted.bam | bcftools call -mv -Ov -o "$sample".vcf
    rm -r "$sample.sorted.bam"
    display_elapsed_time

    # Compress and index VCF file
    echo "Compressing $sample.vcf"
    bgzip "$sample".vcf
    echo "Indexing $sample.vcf"
    tabix -p vcf "$sample".vcf.gz
    display_elapsed_time

    # Generate consensus sequence
    echo "Generating consensus for $sample"
    bcftools consensus -f "$REFERENCE" "$sample".vcf.gz > "$sample".consensus.fasta
    rm -r *vcf*
    display_elapsed_time

    # Append sample name to sequence names
    awk -v sample="$sample" '/^>/{print $0"_"sample;next}{print}' "$sample".consensus.fasta > "$sample".consensus_temp.fasta
    mv "$sample".consensus_temp.fasta "$sample".consensus.fasta
done

echo "Mapping and consensus sequence generation complete."
