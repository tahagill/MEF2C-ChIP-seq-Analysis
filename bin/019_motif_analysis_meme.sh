#!/bin/bash
# 019_motif_analysis_meme.sh
# MEF2C Motif Analysis with MEME-CHIP (chr-prefixed genome)

echo "=== MEF2C MOTIF ANALYSIS ==="

# Directories and files
MOTIF_DIR="results/motifs"
ORIG_GENOME="resources/genomes/hg38/GRCh38_noalt_as.fa"
GENOME_FASTA="resources/genomes/hg38/GRCh38_noalt_as_chr.fa"
CHROM_SIZES="resources/genomes/hg38/chrom.sizes_chr"
PEAK_FILE="results/differential_final_20k/WT_vs_KO_significant_peaks.csv"

mkdir -p "$MOTIF_DIR"

# 1. Create chr-prefixed genome files for consistency
echo "Preparing chr-prefixed genome files..."
awk 'BEGIN {OFS="\t"} {print "chr"$1, $2}' resources/genomes/hg38/chrom.sizes > "$CHROM_SIZES"

if [ ! -f "$GENOME_FASTA" ]; then
    echo "Creating chr-prefixed genome FASTA..."
    sed 's/^>/>chr/' "$ORIG_GENOME" > "$GENOME_FASTA"
fi

# 2. Extract sequences under peaks for motif analysis
echo "Extracting peak sequences for MEME analysis..."
BEDTOOLS_FAIDX="$GENOME_FASTA"
PEAK_FASTA="$MOTIF_DIR/MEF2C_peak_sequences.fasta"
bedtools getfasta -fi "$BEDTOOLS_FAIDX" -bed <(awk 'NR>1 {print $1,$2,$3}' "$PEAK_FILE") -fo "$PEAK_FASTA"

# 3. Run MEME for motif discovery
echo "Running MEME on extracted peak sequences..."
MEME_OUT="$MOTIF_DIR/meme_final"
mkdir -p "$MEME_OUT"

meme "$PEAK_FASTA" \
     -oc "$MEME_OUT" \
     -dna -mod anr -nmotifs 10 -minw 6 -maxw 20 -revcomp

echo "Motif analysis complete."
echo "Results available at: $MEME_OUT/meme.html"
