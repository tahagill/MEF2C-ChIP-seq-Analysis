#!/bin/bash

BASE_DIR=~/pheonix/MEF2C_ChIPseq
RAW_SRA_DIR="$BASE_DIR/data/raw_sra"
FASTQ_DIR="$BASE_DIR/data/raw"
LOG_DIR="$BASE_DIR/logs"

mkdir -p "$RAW_SRA_DIR" "$FASTQ_DIR" "$LOG_DIR"

SAMPLES=("SRR35220292" "SRR35220288" "SRR35220289" "SRR35220284")

cd "$RAW_SRA_DIR" || exit

for SRR in "${SAMPLES[@]}"; do
    echo "Processing $SRR..."
    echo "Processing $SRR..." >> "$LOG_DIR/download_new.log"

    prefetch "$SRR"

    echo "Converting $SRR to FASTQ..."
    echo "Converting $SRR to FASTQ..." >> "$LOG_DIR/download_new.log"
    fasterq-dump "$SRR" -O "$FASTQ_DIR" --split-files

    echo "$SRR finished."
    echo "$SRR finished." >> "$LOG_DIR/download_new.log"
done

echo "All new samples downloaded and converted."
echo "All new samples downloaded and converted." >> "$LOG_DIR/download_new.log"
