#!/usr/bin/env Rscript


library(rtracklayer)
library(Gviz)
library(GenomicRanges)


peaks <- read.csv("results/differential_final_20k/WT_vs_KO_significant_peaks.csv")
peaks_gr <- makeGRangesFromDataFrame(peaks)

=
genes_to_plot <- c("ANK3", "DSCAM", "ARID1B", "NRXN1", "CACNB2")

for(gene in genes_to_plot) {
  cat("Creating track for:", gene, "\n")
  
  # Get gene coordinates from ENSEMBL
  # Plot MEF2C binding peaks
  # Add conservation tracks
  # Save as PDF/PNG
}
