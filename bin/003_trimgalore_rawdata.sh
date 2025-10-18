#!/bin/bash


BASE_DIR=~/pheonix/MEF2C_ChIPseq
RAW_DIR="$BASE_DIR/data/raw"
TRIMMED_DIR="$BASE_DIR/data/processed/trimmed"
QC_DIR="$BASE_DIR/results/qc_new"

mkdir -p "$TRIMMED_DIR"
mkdir -p "$QC_DIR/trimmed"

NEW_SAMPLES=("SRR35220292" "SRR35220288" "SRR35220289" "SRR35220284")

for SAMPLE in "${NEW_SAMPLES[@]}"; do
    echo "Trimming NEW sample: $SAMPLE..."
    
    trim_galore \
        --paired \
        --quality 20 \
        --length 20 \
        --output_dir "$TRIMMED_DIR" \
        --fastqc \
        --fastqc_args "-o $QC_DIR/trimmed" \
        "$RAW_DIR/${SAMPLE}_1.fastq" \
        "$RAW_DIR/${SAMPLE}_2.fastq"
done

echo "Trimming complete for NEW samples only"
