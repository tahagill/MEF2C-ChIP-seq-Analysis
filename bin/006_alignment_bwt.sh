#!/bin/bash


BASE_DIR=~/pheonix/MEF2C_ChIPseq
TRIMMED_DIR="$BASE_DIR/data/processed/trimmed"
ALIGNED_DIR="$BASE_DIR/data/processed/aligned"
GENOME_INDEX="$BASE_DIR/resources/genomes/hg38/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as"

mkdir -p "$ALIGNED_DIR"


NEW_SAMPLES=("SRR35220292" "SRR35220288" "SRR35220289" "SRR35220284")

for SAMPLE in "${NEW_SAMPLES[@]}"; do
    echo "Aligning NEW sample: $SAMPLE..."
    
    bowtie2 \
        -x "$GENOME_INDEX" \
        -1 "$TRIMMED_DIR/${SAMPLE}_1_val_1.fq" \
        -2 "$TRIMMED_DIR/${SAMPLE}_2_val_2.fq" \
        --threads 8 \
        --very-sensitive \
        --no-mixed \
        --no-discordant \
        -S "$ALIGNED_DIR/${SAMPLE}.sam" 2> "$ALIGNED_DIR/${SAMPLE}_alignment_stats.txt"
    
    # Convert to BAM and sort
    samtools view -bS "$ALIGNED_DIR/${SAMPLE}.sam" | samtools sort -o "$ALIGNED_DIR/${SAMPLE}.bam"
    samtools index "$ALIGNED_DIR/${SAMPLE}.bam"
    
    # Clean up SAM file
    rm "$ALIGNED_DIR/${SAMPLE}.sam"
    
    echo "Finished NEW sample: $SAMPLE"
done

echo "Alignment complete for NEW samples only"
