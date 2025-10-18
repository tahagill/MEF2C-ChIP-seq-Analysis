#!/bin/bash

BASE_DIR=~/pheonix/MEF2C_ChIPseq
PEAKS_DIR="$BASE_DIR/results/peaks"

echo "Peak calling quality check"
echo "Started: $(date)"

# Check if peak files exist
echo ""
echo "File existence check:"
for genotype in WT HET KO; do
    peak_file="$PEAKS_DIR/MEF2C_${genotype}_pooled_peaks.narrowPeak"
    summit_file="$PEAKS_DIR/MEF2C_${genotype}_pooled_summits.bed"
    
    if [[ -f "$peak_file" ]]; then
        peak_count=$(wc -l < "$peak_file")
        echo "$genotype: $peak_count peaks"
    else
        echo "$genotype: peaks file missing!"
    fi
    
    if [[ -f "$summit_file" ]]; then
        summit_count=$(wc -l < "$summit_file")
        echo "   Summits: $summit_count"
    else
        echo "   Summits file missing!"
    fi
done

# Check MACS2 model files
echo ""
echo "MACS2 statistics check:"
for genotype in WT HET KO; do
    model_file="$PEAKS_DIR/MEF2C_${genotype}_pooled_model.r"
    if [[ -f "$model_file" ]]; then
        echo "$genotype: Model file present"
    else
        echo "$genotype: Model file missing"
    fi
done

# Preview first few lines of peak files
echo ""
echo "Peak file preview:"
for genotype in WT HET KO; do
    peak_file="$PEAKS_DIR/MEF2C_${genotype}_pooled_peaks.narrowPeak"
    if [[ -f "$peak_file" ]]; then
        echo ""
        echo "--- $genotype peaks (first 3 lines) ---"
        head -3 "$peak_file"
    fi
done

echo ""
echo "Quality metrics to verify:"
echo "1. WT should have the most peaks (MEF2C binding sites)"
echo "2. HET should have fewer peaks than WT"  
echo "3. KO should have the fewest peaks (background noise)"
echo "4. All files should have proper chromosome format (chr1, chr2, etc.)"
echo "5. Peak counts should be reasonable (5,000-50,000 for TF ChIP-seq)"
