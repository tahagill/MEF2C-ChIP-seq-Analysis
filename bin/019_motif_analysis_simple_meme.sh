#!/bin/bash

echo "=== SIMPLE MEF2C MOTIF ANALYSIS ==="

meme results/motifs/MEF2C_peak_sequences.fasta \
    -oc results/motifs/meme_final \
    -dna -mod anr -nmotifs 10 -minw 6 -maxw 20 -revcomp

echo "MOTIF ANALYSIS COMPLETE!"
echo "Check results: results/motifs/meme_final/meme.html"
