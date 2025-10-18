#!/bin/bash

BASE_DIR=~/pheonix/MEF2C_ChIPseq
FASTQ_DIR="$BASE_DIR/data/raw"
QC_DIR="$BASE_DIR/results/qc_new"

mkdir -p "$QC_DIR"

NEW_SAMPLES=("SRR35220292" "SRR35220288" "SRR35220289" "SRR35220284")

for SAMPLE in "${NEW_SAMPLES[@]}"; do
    echo "Running FastQC on $SAMPLE..."
    fastqc "$FASTQ_DIR"/"$SAMPLE"*.fastq -o "$QC_DIR" -t 4
done

multiqc "$QC_DIR" -o "$QC_DIR"

echo "QC complete for NEW samples only"
