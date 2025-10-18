#!/bin/bash

BASE_DIR=~/pheonix/MEF2C_ChIPseq
RAW_SRA_DIR="$BASE_DIR/data/raw_sra"
FASTQ_DIR="$BASE_DIR/data/raw"
LOG_DIR="$BASE_DIR/logs"

mkdir -p "$RAW_SRA_DIR" "$FASTQ_DIR" "$LOG_DIR"

SAMPLES=("SRR35220290" "SRR35220291" "SRR35220286" "SRR35220287" "SRR35220282" "SRR35220283")

cd "$RAW_SRA_DIR" || exit

for SRR in "${SAMPLES[@]}"; do
    echo "Processing $SRR..."
    echo "Processing $SRR..." >> "$LOG_DIR/download.log"

    prefetch "$SRR"

    echo "Converting $SRR to FASTQ..."
    echo "Converting $SRR to FASTQ..." >> "$LOG_DIR/download.log"
    fasterq-dump "$SRR" -O "$FASTQ_DIR" --split-files

    echo "$SRR completed."
    echo "$SRR completed." >> "$LOG_DIR/download.log"
done

echo "All samples downloaded and converted successfully."
echo "All samples downloaded and converted successfully." >> "$LOG_DIR/download.log"
