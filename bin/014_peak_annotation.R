#!/usr/bin/env Rscript
# MEF2C PEAK ANNOTATION - Identify target genes

cat("MEF2C Peak Annotation - Finding target genes\n")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

cat("Loading MEF2C binding sites lost in knockout...\n")

base_dir <- "~/pheonix/MEF2C_ChIPseq"
wt_ko_peaks <- read.csv(file.path(base_dir, "results/differential_final_20k/WT_vs_KO_significant_peaks.csv"))

cat("Number of peaks to annotate:", nrow(wt_ko_peaks), "\n")

# Convert genomic coordinates to GRanges
cat("Converting coordinates to genomic ranges...\n")
peaks_gr <- makeGRangesFromDataFrame(wt_ko_peaks, 
                                     seqnames.field = "seqnames",
                                     start.field = "start", 
                                     end.field = "end")

# Annotate peaks to nearest genes
cat("Annotating peaks to nearest genes...\n")
annotated_peaks <- annotatePeak(peaks_gr, 
                                tssRegion = c(-3000, 3000),  # ±3kb from TSS
                                TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                annoDb = "org.Hs.eg.db")

# Save annotated results
output_dir <- file.path(base_dir, "results/annotations")
dir.create(output_dir, showWarnings = FALSE)
output_file <- file.path(output_dir, "MEF2C_target_genes.csv")
write.csv(as.data.frame(annotated_peaks), output_file, row.names = FALSE)

cat("Annotation complete\n")
cat("MEF2C target genes saved to:", output_file, "\n")

# Summary of genomic distribution
cat("\nGenomic distribution of MEF2C binding:\n")
annotation_table <- table(as.data.frame(annotated_peaks)$annotation)
print(annotation_table)
