#!/usr/bin/env Rscript

# Create smaller, high-confidence peak set

cat("CREATING SMALLER PEAK SET\n")

base_dir <- "~/pheonix/MEF2C_ChIPseq"
peaks_dir <- file.path(base_dir, "results/peaks")

# Read WT peaks (most confident MEF2C binding sites)
wt_peaks <- read.table(file.path(peaks_dir, "MEF2C_WT_pooled_peaks.narrowPeak"), 
                       sep = "\t", header = FALSE)

cat("Original WT peaks:", nrow(wt_peaks), "\n")

# Sort by significance (column 9 = -log10(q-value))
wt_peaks_sorted <- wt_peaks[order(-wt_peaks[,9]), ]

# Take top 30,000 most significant peaks
top_peaks <- wt_peaks_sorted[1:30000, 1:3]

# Save as BED format
write.table(top_peaks, 
            file.path(peaks_dir, "consensus_peaks_small.bed"),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

cat("Created smaller peak set with", nrow(top_peaks), "peaks\n")
