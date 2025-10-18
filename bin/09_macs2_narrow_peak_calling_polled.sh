#!/bin/bash

BASE_DIR=~/pheonix/MEF2C_ChIPseq
ALIGNED_DIR="$BASE_DIR/data/processed/aligned"
PEAKS_DIR="$BASE_DIR/results/peaks"

mkdir -p "$PEAKS_DIR"

echo "Peak calling for pooled replicates"
echo "Started: $(date)"

# WT Peaks (pooled replicates)
echo "Calling peaks for WT samples..."
macs2 callpeak \
    -t "$ALIGNED_DIR/SRR35220290.bam" "$ALIGNED_DIR/SRR35220291.bam" "$ALIGNED_DIR/SRR35220292.bam" \
    -f BAMPE \
    -g hs \
    -n "MEF2C_WT_pooled" \
    --outdir "$PEAKS_DIR" \
    -q 0.05

# HET Peaks (pooled replicates)
echo "Calling peaks for HET samples..."
macs2 callpeak \
    -t "$ALIGNED_DIR/SRR35220287.bam" "$ALIGNED_DIR/SRR35220288.bam" "$ALIGNED_DIR/SRR35220289.bam" \
    -f BAMPE \
    -g hs \
    -n "MEF2C_HET_pooled" \
    --outdir "$PEAKS_DIR" \
    -q 0.05

# KO Peaks (pooled replicates)
echo "Calling peaks for KO samples..."
macs2 callpeak \
    -t "$ALIGNED_DIR/SRR35220282.bam" "$ALIGNED_DIR/SRR35220283.bam" "$ALIGNED_DIR/SRR35220284.bam" \
    -f BAMPE \
    -g hs \
    -n "MEF2C_KO_pooled" \
    --outdir "$PEAKS_DIR" \
    -q 0.05

echo "Peak calling finished: $(date)"
echo "Output directory: $PEAKS_DIR"
