\#!/usr/bin/env Rscript
# 016_pathway_analysis.R
# Functional Enrichment Analysis for MEF2C Target Genes

# Load required libraries
suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(DOSE)
  library(enrichplot)
  library(dplyr)
})

# Set up output directory
output_dir <- "../results/pathway_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Starting pathway analysis for MEF2C target genes\n")

# Read target genes
gene_data <- read.csv("../results/annotations/MEF2C_target_genes_FINAL.csv", stringsAsFactors = FALSE)
target_genes <- unique(na.omit(gene_data$SYMBOL))
cat("Number of unique target genes:", length(target_genes), "\n")

# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(target_genes, fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = "org.Hs.eg.db")
gene_list <- gene_entrez$ENTREZID
cat("Converted genes to Entrez IDs:", length(gene_list), "\n")

# Gene Ontology Enrichment
cat("Performing Gene Ontology (GO) enrichment analysis\n")
go_bp <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                  qvalueCutoff = 0.10, readable = TRUE)
go_mf <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                  ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                  qvalueCutoff = 0.10, readable = TRUE)
go_cc <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                  ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                  qvalueCutoff = 0.10, readable = TRUE)

# KEGG Pathway Enrichment
cat("Performing KEGG pathway enrichment analysis\n")
kegg_enrich <- enrichKEGG(gene = gene_list, organism = 'hsa', pvalueCutoff = 0.05, qvalueCutoff = 0.10)

# Disease Ontology Enrichment
cat("Performing Disease Ontology enrichment analysis\n")
do_enrich <- enrichDO(gene = gene_list, pvalueCutoff = 0.05, qvalueCutoff = 0.10, readable = TRUE)

# Save results
cat("Saving enrichment results\n")
if (!is.null(go_bp) && nrow(go_bp) > 0) write.csv(go_bp, file.path(output_dir, "GO_BP_Enrichment.csv"), row.names = FALSE)
if (!is.null(go_mf) && nrow(go_mf) > 0) write.csv(go_mf, file.path(output_dir, "GO_MF_Enrichment.csv"), row.names = FALSE)
if (!is.null(go_cc) && nrow(go_cc) > 0) write.csv(go_cc, file.path(output_dir, "GO_CC_Enrichment.csv"), row.names = FALSE)
if (!is.null(kegg_enrich) && nrow(kegg_enrich) > 0) write.csv(kegg_enrich, file.path(output_dir, "KEGG_Enrichment.csv"), row.names = FALSE)
if (!is.null(do_enrich) && nrow(do_enrich) > 0) write.csv(do_enrich, file.path(output_dir, "DO_Enrichment.csv"), row.names = FALSE)

# Visualization
cat("Generating plots\n")
create_barplot <- function(enrich_result, title, filename) {
  if (is.null(enrich_result) || nrow(enrich_result) == 0) return()
  p <- barplot(enrich_result, showCategory = 15, font.size = 10) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(output_dir, filename), plot = p, width = 12, height = 8)
}

create_barplot(go_bp, "GO Biological Process", "GO_BP_Barplot.png")
create_barplot(go_mf, "GO Molecular Function", "GO_MF_Barplot.png")
create_barplot(go_cc, "GO Cellular Component", "GO_CC_Barplot.png")
create_barplot(kegg_enrich, "KEGG Pathways", "KEGG_Barplot.png")
create_barplot(do_enrich, "Disease Ontology", "DO_Barplot.png")

# Dotplots for GO BP
if (!is.null(go_bp) && nrow(go_bp) > 0) {
  p <- dotplot(go_bp, showCategory = 15) + ggtitle("GO Biological Process - DotPlot")
  ggsave(file.path(output_dir, "GO_BP_Dotplot.png"), plot = p, width = 12, height = 8)
}

cat("Pathway analysis complete\n")
cat("Results saved in:", output_dir, "\n")

# Quick summary
cat("\nSummary:\n")
if (!is.null(go_bp) && nrow(go_bp) > 0) cat("GO BP terms:", nrow(go_bp), "Top 3:", paste(head(go_bp$Description, 3), collapse=", "), "\n")
if (!is.null(kegg_enrich) && nrow(kegg_enrich) > 0) cat("KEGG pathways:", nrow(kegg_enrich), "Top 3:", paste(head(kegg_enrich$Description, 3), collapse=", "), "\n")
if (!is.null(do_enrich) && nrow(do_enrich) > 0) cat("Disease Ontology terms:", nrow(do_enrich), "Top 3:", paste(head(do_enrich$Description, 3), collapse=", "), "\n")
