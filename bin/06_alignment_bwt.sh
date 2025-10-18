#!/bin/bash


BASE_DIR=~/pheonix/MEF2C_ChIPseq
TRIMMED_DIR="$BASE_DIR/data/processed/trimmed"
ALIGNED_DIR="$BASE_DIR/data/processed/aligned"
GENOME_INDEX="$BASE_DIR/resources/genomes/hg38/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as"

mkdir -p "$ALIGNED_DIR"

for file in "$TRIMMED_DIR"/*_1_val_1.fq; do
    base=$(basename "$file" _1_val_1.fq)
    echo "Aligning $base..."
    
    bowtie2 \
        -x "$GENOME_INDEX" \
        -1 "$TRIMMED_DIR/${base}_1_val_1.fq" \
        -2 "$TRIMMED_DIR/${base}_2_val_2.fq" \
        --threads 8 \
        --very-sensitive \
        --no-mixed \
        --no-discordant \
        -S "$ALIGNED_DIR/${base}.sam" 2> "$ALIGNED_DIR/${base}_alignment_stats.txt"
    
    # Convert to BAM and sort
    samtools view -bS "$ALIGNED_DIR/${base}.sam" | samtools sort -o "$ALIGNED_DIR/${base}.bam"
    samtools index "$ALIGNED_DIR/${base}.bam"
    
    # Clean up SAM file
    rm "$ALIGNED_DIR/${base}.sam"
    
    echo "Finished $base"
done

echo "Alignment complete for all samples"
