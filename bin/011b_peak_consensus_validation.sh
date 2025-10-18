#!/bin/bash

BASE_DIR=~/pheonix/MEF2C_ChIPseq
PEAKS_DIR="$BASE_DIR/results/peaks"

echo "Consensus peaks validation"
echo "Started: $(date)"

# Check if consensus file exists
if [[ -f "$PEAKS_DIR/consensus_peaks.bed" ]]; then
    consensus_count=$(wc -l < "$PEAKS_DIR/consensus_peaks.bed")
    echo "Consensus peaks file exists"
    echo "Total consensus peaks: $consensus_count"
else
    echo "Consensus peaks file missing!"
    exit 1
fi

# Check file format
echo ""
echo "File format check:"
echo "First 5 lines of consensus peaks:"
head -5 "$PEAKS_DIR/consensus_peaks.bed"

# Check chromosome distribution
echo ""
echo "Chromosome distribution (top 10):"
cut -f1 "$PEAKS_DIR/consensus_peaks.bed" | sort | uniq -c | sort -rn | head -10

# Check peak widths
echo ""
echo "Peak width statistics:"
awk '{print $3-$2}' "$PEAKS_DIR/consensus_peaks.bed" | \
awk '{sum+=$1; array[NR]=$1} END {
    for(x=1;x<=NR;x++){sumsq+=((array[x]-(sum/NR))^2)}
    print "Mean: " sum/NR " bp"
    print "StdDev: " sqrt(sumsq/NR) " bp"
    print "Min: " min " bp"
    print "Max: " max " bp"
}' min=$(awk '{print $3-$2}' "$PEAKS_DIR/consensus_peaks.bed" | sort -n | head -1) \
   max=$(awk '{print $3-$2}' "$PEAKS_DIR/consensus_peaks.bed" | sort -n | tail -1)

# Check overlap with original peaks
echo ""
echo "Overlap with original peak sets:"
for genotype in WT HET KO; do
    original_peaks="$PEAKS_DIR/MEF2C_${genotype}_pooled_peaks.narrowPeak"
    if [[ -f "$original_peaks" ]]; then
        original_count=$(wc -l < "$original_peaks")
        overlap_count=$(bedtools intersect -a "$PEAKS_DIR/consensus_peaks.bed" -b "$original_peaks" -u | wc -l)
        overlap_percent=$(echo "scale=2; $overlap_count * 100 / $original_count" | bc)
        echo "$genotype: $overlap_count/$original_count peaks in consensus ($overlap_percent%)"
    fi
done

echo ""
echo "Expected results:"
echo "1. Consensus peaks should be 100,000-150,000 regions"
echo "2. Majority should be on main chromosomes (chr1-22, X, Y)"
echo "3. Mean peak width should be 200-1000 bp"
echo "4. High overlap with WT peaks (>90%)"
echo "5. Lower overlap with KO peaks (<50%)"
