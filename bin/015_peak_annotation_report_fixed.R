#!/usr/bin/env Rscript
# Generate a corrected MEF2C target gene report using existing SYMBOL column

cat("MEF2C Target Gene Report\n")

base_dir <- "~/pheonix/MEF2C_ChIPseq"
annotated_file <- file.path(base_dir, "results/annotations/MEF2C_target_genes.csv")
report_file <- file.path(base_dir, "results/annotations/CORRECT_ANNOTATION_REPORT.txt")

# Read annotated peaks
data <- read.csv(annotated_file)
cat("Analyzing", nrow(data), "MEF2C binding sites\n")

# Start writing to report
sink(report_file)

cat("MEF2C ChIP-seq Target Genes Report\n")
cat("=================================\n")
cat("Generated:", date(), "\n\n")

cat("Experiment Summary:\n")
cat("-------------------\n")
cat("Total MEF2C binding sites lost in KO:", nrow(data), "\n")
cat("Unique genes targeted:", length(unique(data$SYMBOL[!is.na(data$SYMBOL)])), "\n\n")

cat("Genomic Distribution:\n")
cat("---------------------\n")
# Clean annotation categories
clean_anno <- gsub(" \\(.*\\)", "", data$annotation)
clean_anno <- gsub("Exon.*", "Exon", clean_anno)
clean_anno <- gsub("Intron.*", "Intron", clean_anno)
anno_table <- table(clean_anno)

for (item in names(anno_table)) {
  cat(sprintf("%-25s: %4d peaks (%5.1f%%)\n", 
              item, anno_table[item], 
              anno_table[item]/nrow(data)*100))
}
cat("\n")

cat("Top 20 MEF2C Target Genes:\n")
cat("---------------------------\n")
valid_genes <- data$SYMBOL[!is.na(data$SYMBOL) & data$SYMBOL != ""]
gene_counts <- table(valid_genes)
top_genes <- head(sort(gene_counts, decreasing = TRUE), 20)

for (i in 1:length(top_genes)) {
  cat(sprintf("%-20s: %2d binding sites\n", names(top_genes)[i], top_genes[i]))
}
cat("\n")

cat("Key Biological Insights:\n")
cat("------------------------\n")
promoter_peaks <- sum(grepl("Promoter", clean_anno))
cat("•", promoter_peaks, "promoter peaks → direct transcriptional regulation\n")

intronic_peaks <- sum(grepl("Intron", clean_anno))
cat("•", intronic_peaks, "intronic peaks → potential enhancer regulation\n")

intergenic_peaks <- sum(grepl("Intergenic", clean_anno))
cat("•", intergenic_peaks, "intergenic peaks → long-range regulation\n\n")

cat("Sample of MEF2C Target Genes:\n")
cat("-----------------------------\n")
sample_genes <- head(unique(valid_genes), 15)
cat(paste(sample_genes, collapse=", "), "\n")

sink()

cat("Report generated successfully.\n")
cat("Report file:", report_file, "\n")

# Quick terminal summary
cat("\nQuick Summary:\n")
cat("Unique genes:", length(unique(valid_genes)), "\n")
cat("Top gene:", names(which.max(gene_counts)), "with", max(gene_counts), "binding sites\n")

