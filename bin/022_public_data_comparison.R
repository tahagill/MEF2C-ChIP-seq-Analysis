#!/usr/bin/env Rscript

library(GenomicRanges)
library(ggplot2)

cat("MEF2C Microglia Specificity Analysis (4 Cell Types)\n")

# Load microglial MEF2C peaks
cat("1. Loading microglial MEF2C peaks\n")
your_peaks <- read.table("results/differential_final_20k/WT_vs_KO_significant_peaks.csv", header = TRUE, sep = ",")
cat("Number of significant microglial peaks:", nrow(your_peaks), "\n")

# Convert to GRanges object
your_gr <- GRanges(
  seqnames = your_peaks$seqnames,
  ranges = IRanges(start = your_peaks$start, end = your_peaks$end),
  strand = "*"
)

# Simple BED reader (first 3 columns only)
read_bed_simple <- function(file_path) {
  if(!file.exists(file_path)) {
    cat("File not found:", file_path, "\n")
    return(NULL)
  }
  cat("Reading:", file_path, "\n")
  bed_data <- read.table(file_path, header = FALSE, sep = "\t")
  gr <- GRanges(
    seqnames = bed_data[,1],
    ranges = IRanges(start = bed_data[,2], end = bed_data[,3]),
    strand = "*"
  )
  cat("Number of peaks loaded:", length(gr), "\n")
  return(gr)
}

# Load ENCODE MEF2C data for multiple cell types
cat("2. Loading ENCODE MEF2C peaks\n")
encode_files <- c(
  cardio = "results/public_data/comparison/MEF2C_cardio.bed",
  lymph = "results/public_data/comparison/MEF2C_lymph.bed",
  myotube = "results/public_data/comparison/MEF2C_myotube.bed",
  monocyte = "results/public_data/comparison/MEF2C_monocyte.bed",
  hepg2 = "results/public_data/comparison/MEF2C_hepg2.bed"
)

encode_peaks <- list()
for (cell_type in names(encode_files)) {
  peaks <- read_bed_simple(encode_files[cell_type])
  if(!is.null(peaks)) {
    encode_peaks[[cell_type]] <- peaks
  }
}
cat("Successfully loaded", length(encode_peaks), "cell type datasets\n")

# Calculate overlaps (minimum 50bp)
cat("3. Calculating overlaps\n")
overlaps <- list()
for (cell_type in names(encode_peaks)) {
  overlaps[[cell_type]] <- findOverlaps(your_gr, encode_peaks[[cell_type]], minoverlap = 50)
  cat("Overlap with", cell_type, ":", length(overlaps[[cell_type]]), "peaks\n")
}

# Identify microglia-specific peaks
all_overlapped <- unique(unlist(lapply(overlaps, queryHits)))
microglia_specific_peaks <- your_gr[!seq_along(your_gr) %in% all_overlapped]

# Create comparison table
cell_type_names <- c(
  "cardio" = "Cardiomyocytes",
  "lymph" = "Lymphoblastoid", 
  "myotube" = "Skeletal Myotubes",
  "monocyte" = "Monocytes",
  "hepg2" = "HepG2 Liver"
)

comparison_results <- data.frame(
  Cell_Type = c("Microglia (Your Data)", cell_type_names[names(encode_peaks)]),
  Total_Peaks = c(length(your_gr), sapply(encode_peaks, length)),
  Overlap_With_Microglia = c(NA, sapply(overlaps, length)),
  Percent_Overlap = c(NA, round(sapply(overlaps, length) / length(your_gr) * 100, 2))
)

cat("\nRESULTS:\n")
print(comparison_results)
cat("Microglia-specific peaks:", length(microglia_specific_peaks), "/", length(your_gr), 
    "(", round(length(microglia_specific_peaks)/length(your_gr)*100,1), "%)\n")

# Save results
write.csv(comparison_results, "results/public_data/comparison/MEF2C_cell_type_comparison.csv", row.names = FALSE)

microglia_specific_df <- data.frame(
  seqnames = seqnames(microglia_specific_peaks),
  start = start(microglia_specific_peaks),
  end = end(microglia_specific_peaks)
)
write.csv(microglia_specific_df, "results/public_data/comparison/MEF2C_microglia_specific_peaks.csv", row.names = FALSE)

# Visualization
cat("4. Creating bar plot visualization\n")
biological_order <- c("Monocytes", "Lymphoblastoid", "Skeletal Myotubes", "Cardiomyocytes", "HepG2 Liver")
available_types <- comparison_results$Cell_Type[-1]
available_types <- available_types[order(match(available_types, biological_order))]

comparison_ordered <- comparison_results
comparison_ordered$Cell_Type <- factor(comparison_ordered$Cell_Type, 
                                       levels = c("Microglia (Your Data)", available_types))

p <- ggplot(comparison_results[-1,], aes(x = Cell_Type, y = Percent_Overlap, fill = Cell_Type)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = paste0(Percent_Overlap, "%")), vjust = -0.5, size = 4, fontface = "bold") +
  labs(title = "MEF2C Binding Overlap: Microglia vs Multiple Cell Types",
       subtitle = paste("Microglia-specific peaks:", length(microglia_specific_peaks), "/", length(your_gr), 
                        "(", round(length(microglia_specific_peaks)/length(your_gr)*100,1), "%)"),
       x = "Cell Type", 
       y = "% of Microglial Peaks Overlapping",
       caption = "Cell types ordered by biological relatedness to microglia") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(face = "bold", size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.caption = element_text(size = 9, color = "gray50")) +
  ylim(0, max(comparison_results$Percent_Overlap[-1], na.rm=TRUE) * 1.4)

ggsave("results/public_data/comparison/MEF2C_cell_type_comparison.pdf", p, width = 10, height = 7)

# Overlap summary
overlap_summary <- data.frame(
  Cell_Type = names(encode_peaks),
  Overlap_Count = sapply(overlaps, length),
  Overlap_Percent = round(sapply(overlaps, length) / length(your_gr) * 100, 2)
)
write.csv(overlap_summary, "results/public_data/comparison/MEF2C_overlap_summary.csv", row.names = FALSE)

# Biological interpretation
cat("\nBIOLOGICAL INTERPRETATION:\n")
microglia_specific_percent <- round(length(microglia_specific_peaks)/length(your_gr)*100,1)

if(microglia_specific_percent > 90) {
  cat("Extreme microglia specificity. MEF2C binding is largely unique to microglia.\n")
} else if(microglia_specific_percent > 70) {
  cat("Strong microglia specificity. Most MEF2C binding is cell-type specific.\n")
} else if(microglia_specific_percent > 40) {
  cat("Moderate microglia specificity. MEF2C has shared and microglia-specific functions.\n")
} else {
  cat("Conserved binding patterns. MEF2C binding is largely shared across cell types.\n")
}

# Gradient analysis
if(length(encode_peaks) >= 3) {
  cat("\nGRADIENT ANALYSIS:\n")
  overlap_rates <- sapply(overlaps, length) / length(your_gr) * 100
  cat("Overlap rates across cell types:\n")
  for(i in 1:length(overlap_rates)) {
    cat(" ", names(overlap_rates)[i], ":", round(overlap_rates[i], 2), "%\n")
  }
  
  if("monocyte" %in% names(overlap_rates) && "cardio" %in% names(overlap_rates)) {
    if(overlap_rates["monocyte"] > overlap_rates["cardio"]) {
      cat("Biological gradient observed: Higher overlap with monocytes (immune lineage)\n")
    }
  }
}

cat("\nAnalysis complete. Results saved to results/public_data/comparison/\n")
cat("Key files:\n")
cat("- MEF2C_cell_type_comparison.csv\n")
cat("- MEF2C_microglia_specific_peaks.csv\n")
cat("- MEF2C_overlap_summary.csv\n")
cat("- MEF2C_cell_type_comparison.pdf\n")
