#!/usr/bin/env Rscript

# MEF2C NDD ENRICHMENT ANALYSIS - FINAL FIX

cat("MEF2C NDD INTEGRATION ANALYSIS\n")

suppressMessages({
  library(dplyr)
  library(readr) 
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})

# Setup results directory
results_dir <- "results/disease_integration"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Step 1: Load MEF2C target genes
cat("\nStep 1: Loading MEF2C target genes\n")

target_genes <- read_csv("results/annotations/MEF2C_target_genes_FINAL.csv", show_col_types = FALSE)
cat("Number of rows loaded:", nrow(target_genes), "\n")

target_ensembl <- unique(na.omit(target_genes$ENSEMBL))
cat("Unique MEF2C target genes (ENSEMBL):", length(target_ensembl), "\n")

cat("Converting ENSEMBL IDs to ENTREZ IDs\n")
target_entrez <- bitr(target_ensembl, 
                      fromType = "ENSEMBL", 
                      toType = "ENTREZID", 
                      OrgDb = "org.Hs.eg.db")
cat("Number of genes successfully converted to ENTREZ IDs:", nrow(target_entrez), "\n")
target_entrez_ids <- target_entrez$ENTREZID

# Step 2: Define curated NDD gene sets
cat("\nStep 2: Defining NDD gene sets\n")

ndd_genes_curated <- list(
  "ASD_SFARI_HighConf" = c("SHANK3", "NLGN3", "NLGN4X", "NRXN1", "CNTNAP2", "CHD8", 
                          "DYRK1A", "SCN2A", "ADNP", "ARID1B", "SYNGAP1", "GRIN2B",
                          "PTEN", "TSC1", "TSC2", "MECP2", "FMR1"),
  
  "Intellectual_Disability" = c("MECP2", "FMR1", "UBE3A", "CUL4B", "ARX", "SYNGAP1",
                               "STXBP1", "SCN1A", "SCN2A", "KCNQ2", "CDKL5",
                               "FOXG1", "TCF4", "EHMT1", "KANSL1"),
  
  "Epilepsy_Core" = c("SCN1A", "SCN2A", "SCN8A", "KCNQ2", "KCNQ3", "CDKL5",
                     "PCDH19", "STXBP1", "PRRT2", "SPTAN1", "ALG13", "GRIN2A",
                     "GRIN2B", "CACNA1A"),
  
  "Synaptopathy" = c("ANK3", "NLGN1", "SH3GL2", "DSCAM", "NRXN1", "NLGN3",
                    "SYNGAP1", "GRIN2B", "GRIN2A", "CACNA1C", "CACNB2",
                    "SHANK1", "SHANK2", "SHANK3")
)

cat("Defined", length(ndd_genes_curated), "NDD categories\n")

# Step 3: Enrichment analysis
cat("\nStep 3: Running enrichment analysis\n")

background_genes <- keys(org.Hs.eg.db, keytype = "ENTREZID")
cat("Background gene set size:", length(background_genes), "\n")

calculate_enrichment <- function(ndd_gene_symbols, category_name, targets, background) {
  
  ndd_entrez <- bitr(ndd_gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID
  if (length(ndd_entrez) == 0) return(NULL)
  
  overlap <- intersect(targets, ndd_entrez)
  overlap_count <- length(overlap)
  if (overlap_count == 0) return(NULL)
  
  overlap_symbols <- bitr(overlap, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")$SYMBOL
  
  m <- length(intersect(background, ndd_entrez))
  n <- length(background) - m
  k <- length(targets)
  
  enrichment_p <- phyper(
    q = overlap_count - 1,
    m = m,
    n = n,
    k = k,
    lower.tail = FALSE
  )
  
  expected <- (m / length(background)) * k
  fold_enrichment <- overlap_count / expected
  
  return(data.frame(
    Category = category_name,
    NDD_Gene_Set_Size = length(ndd_gene_symbols),
    NDD_Genes_In_Background = m,
    MEF2C_Target_Genes = k,
    Overlap_Count = overlap_count,
    Expected_Overlap = round(expected, 2),
    Fold_Enrichment = round(fold_enrichment, 2),
    P_Value = enrichment_p,
    Overlap_Genes = paste(overlap_symbols, collapse = ", "),
    stringsAsFactors = FALSE
  ))
}

enrichment_results <- list()
for (category in names(ndd_genes_curated)) {
  cat("Analyzing category:", category, "\n")
  result <- calculate_enrichment(ndd_genes_curated[[category]], category, target_entrez_ids, background_genes)
  if (!is.null(result)) {
    enrichment_results[[category]] <- result
  }
}

# Step 4: Compile results and adjust for multiple testing
cat("\nStep 4: Results\n")

if (length(enrichment_results) == 0) {
  cat("No significant overlaps found.\n")
} else {
  final_results <- bind_rows(enrichment_results) %>%
    mutate(
      FDR_Adjusted = p.adjust(P_Value, method = "BH"),
      Significance = case_when(
        FDR_Adjusted < 0.001 ~ "***",
        FDR_Adjusted < 0.01 ~ "**",
        FDR_Adjusted < 0.05 ~ "*", 
        FDR_Adjusted < 0.1 ~ ".",
        TRUE ~ "ns"
      )
    ) %>%
    arrange(P_Value)
  
  write_csv(final_results, file.path(results_dir, "NDD_Enrichment_Results.csv"))
  print(final_results[, c("Category", "Overlap_Count", "Fold_Enrichment", "P_Value", "FDR_Adjusted", "Significance")])
  
  plot_data <- final_results %>%
    mutate(
      Log10_PValue = -log10(P_Value),
      Category = factor(Category, levels = Category[order(P_Value)])
    )
  
  enrichment_plot <- ggplot(plot_data, aes(x = Fold_Enrichment, y = Category)) +
    geom_point(aes(size = Overlap_Count, color = Log10_PValue)) +
    scale_color_gradient(low = "blue", high = "red", name = "-log10(P-value)") +
    scale_size_continuous(name = "Overlap Count", range = c(3, 8)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    labs(
      title = "MEF2C Target Gene Enrichment in Neurodevelopmental Disorders",
      subtitle = paste("Analysis of", length(target_entrez_ids), "MEF2C target genes"),
      x = "Fold Enrichment",
      y = "NDD Category"
    ) +
    theme_minimal()
  
  ggsave(file.path(results_dir, "NDD_Enrichment_Plot.pdf"), enrichment_plot, width = 10, height = 6)
  cat("Enrichment plot saved\n")
}

# Step 5: High-confidence MEF2C-NDD targets
cat("\nStep 5: High-confidence targets\n")

high_confidence <- target_genes %>%
  group_by(ENSEMBL, SYMBOL) %>%
  summarise(Binding_Sites = n(), .groups = 'drop') %>%
  filter(!is.na(ENSEMBL), !is.na(SYMBOL)) %>%
  arrange(desc(Binding_Sites)) %>%
  mutate(
    In_NDD_Databases = SYMBOL %in% unlist(ndd_genes_curated),
    NDD_Categories = sapply(SYMBOL, function(gene) {
      cats <- names(ndd_genes_curated)[sapply(ndd_genes_curated, function(x) gene %in% x)]
      if (length(cats) > 0) paste(cats, collapse = "; ") else NA
    })
  )

write_csv(high_confidence, file.path(results_dir, "MEF2C_Targets_NDD_Annotation.csv"))

high_confidence_ndd <- high_confidence %>% 
  filter(In_NDD_Databases) %>%
  arrange(desc(Binding_Sites))

write_csv(high_confidence_ndd, file.path(results_dir, "High_Confidence_NDD_Targets.csv"))

cat("Number of high-confidence MEF2C-NDD targets:", nrow(high_confidence_ndd), "\n")
if (nrow(high_confidence_ndd) > 0) {
  print(high_confidence_ndd[, c("SYMBOL", "Binding_Sites", "NDD_Categories")])
}

cat("\nAnalysis complete. Results are saved in", results_dir, "\n")
