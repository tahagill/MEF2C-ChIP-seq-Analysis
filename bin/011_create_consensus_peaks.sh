#!/bin/bash

BASE_DIR=~/pheonix/MEF2C_ChIPseq
PEAKS_DIR="$BASE_DIR/results/peaks"

echo "=== CREATING CONSENSUS PEAK SET ==="
echo "Started: $(date)"

# Combine all peaks and create consensus set
cat \
    "$PEAKS_DIR/MEF2C_WT_pooled_peaks.narrowPeak" \
    "$PEAKS_DIR/MEF2C_HET_pooled_peaks.narrowPeak" \
    "$PEAKS_DIR/MEF2C_KO_pooled_peaks.narrowPeak" | \
    bedtools sort | \
    bedtools merge -i - > "$PEAKS_DIR/consensus_peaks.bed"

echo "Consensus peaks created: $PEAKS_DIR/consensus_peaks.bed"
echo "Number of consensus peaks: $(wc -l < "$PEAKS_DIR/consensus_peaks.bed")"
echo "Finished: $(date)"
