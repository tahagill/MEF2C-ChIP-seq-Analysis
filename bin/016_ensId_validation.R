#!/usr/bin/env Rscript
# Final check of annotated MEF2C target data

cat("Checking final annotation data\n")

base_dir <- "~/pheonix/MEF2C_ChIPseq"
final_file <- file.path(base_dir, "results/annotations/MEF2C_target_genes_FINAL.csv")
report_file <- file.path(base_dir, "results/annotations/ANNOTATION_REPORT_FINAL.txt")

# Read the final annotated data
final_data <- read.csv(final_file)

cat("Final data overview:\n")
cat("-------------------\n")
cat("Total rows:", nrow(final_data), "\n")
cat("Columns:", paste(colnames(final_data), collapse=", "), "\n\n")

cat("Gene symbol summary:\n")
cat("-------------------\n")
gene_symbols <- final_data$hgnc_symbol
cat("Peaks with gene symbols:", sum(!is.na(gene_symbols) & gene_symbols != ""), "\n")
cat("Peaks without symbols:", sum(is.na(gene_symbols) | gene_symbols == ""), "\n\n")

cat("Top 10 MEF2C target genes:\n")
cat("--------------------------\n")
valid_genes <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]
top_genes <- head(sort(table(valid_genes), decreasing = TRUE), 10)
for (i in 1:length(top_genes)) {
    cat(sprintf("%-20s: %2d binding sites\n", names(top_genes)[i], top_genes[i]))
}

cat("\nSample of annotated peaks:\n")
cat("-------------------------\n")
has_symbols <- final_data[!is.na(final_data$hgnc_symbol) & final_data$hgnc_symbol != "", ]
print(head(has_symbols[, c("seqnames", "start", "end", "hgnc_symbol", "annotation")], 10))

cat("\nCheck complete.\n")
