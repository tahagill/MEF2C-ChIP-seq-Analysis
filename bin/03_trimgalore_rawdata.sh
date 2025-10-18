#!/bin/bash

BASE_DIR=~/pheonix/MEF2C_ChIPseq
RAW_DIR="$BASE_DIR/data/raw"
TRIMMED_DIR="$BASE_DIR/data/processed/trimmed"
QC_DIR="$BASE_DIR/results/qc"

mkdir -p "$TRIMMED_DIR"
mkdir -p "$QC_DIR/trimmed"

for file in "$RAW_DIR"/*_1.fastq; do
    base=$(basename "$file" _1.fastq)
    echo "Trimming $base..."
    
    trim_galore \
        --paired \
        --quality 20 \
        --length 20 \
        --output_dir "$TRIMMED_DIR" \
        --fastqc \
        --fastqc_args "-o $QC_DIR/trimmed" \
        "$RAW_DIR/${base}_1.fastq" \
        "$RAW_DIR/${base}_2.fastq"
done
