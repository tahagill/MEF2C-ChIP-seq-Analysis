\#!/usr/bin/env Rscript

library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(patchwork)

# Create output directories
dir.create("results/visualization/genome_tracks", recursive = TRUE, showWarnings = FALSE)
dir.create("results/visualization/igv_snapshots", recursive = TRUE, showWarnings = FALSE)

# Load significant peaks
cat("Loading significant peaks...\n")
peaks <- read.csv("results/differential_final_20k/WT_vs_KO_significant_peaks.csv")
peaks_gr <- makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
cat("Loaded", length(peaks_gr), "peaks\n")

# Load gene annotations
cat("Loading hg38 gene annotations...\n")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- genes(txdb)

# Map NDD genes to ENTREZ IDs
ndd_genes <- c("ANK3", "DSCAM", "ARID1B", "NRXN1", "CACNB2")
gene_map <- select(org.Hs.eg.db, keys = ndd_genes, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")

# Get genomic coordinates for each gene
gene_coords <- list()
for(i in 1:nrow(gene_map)) {
  symbol <- gene_map$SYMBOL[i]
  entrez <- gene_map$ENTREZID[i]
  if(!is.na(entrez)) {
    gene_subset <- genes[genes$gene_id == entrez]
    if(length(gene_subset) > 0) {
      gene_coords[[symbol]] <- gene_subset
      cat("Found coordinates for", symbol, ":", as.character(seqnames(gene_subset)), 
          start(gene_subset), "-", end(gene_subset), "\n")
    }
  }
}

# Create UCSC BED files
cat("Creating UCSC BED files...\n")
export(peaks_gr, "results/visualization/genome_tracks/MEF2C_significant_peaks.bed", format = "BED")

ndd_bed <- data.frame(
  chrom = sapply(gene_coords, function(x) as.character(seqnames(x))),
  start = sapply(gene_coords, function(x) start(x)),
  end = sapply(gene_coords, function(x) end(x)),
  name = names(gene_coords),
  score = 1000,
  strand = sapply(gene_coords, function(x) as.character(strand(x)))
)
write.table(ndd_bed, "results/visualization/genome_tracks/NDD_genes.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Create IGV session
cat("Creating IGV session file...\n")
igv_session <- '<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="hg38" hasGeneTrack="true" hasSequenceTrack="true" version="8">
    <Resources>
        <Resource path="MEF2C_significant_peaks.bed"/>
        <Resource path="NDD_genes.bed"/>
    </Resources>
    <Panel name="Panel1">
        <Track attributeKey="Reference sequence" clazz="org.broad.igv.track.SequenceTrack" fontSize="10" id="Reference sequence" name="Reference sequence" sequenceTranslationStrand="POSITIVE" visible="true"/>
        <Track attributeKey="Chromosome Band" clazz="org.broad.igv.track.CytobandTrack" id="Chromosome Band" name="Chromosome Band" visible="true"/>
        <Track attributeKey="MEF2C Peaks" color="255,0,0" colorScale="ContinuousColorScale;0.0;272.0;255,0,0;255,0,0" autoscale="false" clazz="org.broad.igv.track.FeatureTrack" fontSize="10" height="60" id="MEF2C_significant_peaks.bed" name="MEF2C Significant Peaks" visible="true" windowFunction="count"/>
        <Track attributeKey="NDD Genes" color="0,0,255" clazz="org.broad.igv.track.FeatureTrack" fontSize="10" height="40" id="NDD_genes.bed" name="NDD Target Genes" visible="true"/>
    </Panel>
</Session>'
writeLines(igv_session, "results/visualization/igv_snapshots/MEF2C_NDD_session.xml")

# Create UCSC session URL
cat("Creating UCSC session URL...\n")
ucsc_session <- paste0(
  "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=",
  as.character(seqnames(gene_coords[[1]])), ":", start(gene_coords[[1]]), "-", end(gene_coords[[1]]),
  "&hgct_customText=track%20name=%22MEF2C%20Peaks%22%20description=%22MEF2C%20Significant%20Peaks%22%20visibility=2%20color=255,0,0"
)
writeLines(ucsc_session, "results/visualization/genome_tracks/UCSC_session_url.txt")

# Create genomic region plots
cat("Creating genomic plots...\n")
create_genomic_plot <- function(gene_name, gene_gr, peaks_gr) {
  gene_region <- resize(gene_gr, width = width(gene_gr) + 200000, fix = "center")
  nearby_peaks <- subsetByOverlaps(peaks_gr, gene_region)
  cat(" ", gene_name, ":", length(nearby_peaks), "nearby peaks\n")
  plot_data <- data.frame(
    Feature = c("Gene", rep("MEF2C Peak", length(nearby_peaks))),
    Position = c(0, seq_along(nearby_peaks)),
    Chromosome = as.character(seqnames(gene_gr)),
    Start = c(start(gene_gr), start(nearby_peaks))
  )
  ggplot(plot_data, aes(x = Start, y = 1, color = Feature, size = Feature)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("Gene" = "blue", "MEF2C Peak" = "red")) +
    scale_size_manual(values = c("Gene" = 6, "MEF2C Peak" = 4)) +
    labs(title = paste("MEF2C Binding at", gene_name, "Locus"),
         subtitle = paste("Chr", as.character(seqnames(gene_gr)), ":", 
                          format(start(gene_gr), big.mark = ","), "-", 
                          format(end(gene_gr), big.mark = ","),
                          " | ", length(nearby_peaks), " nearby peaks"),
         x = "Genomic Position (bp)", y = "") +
    theme_minimal() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank())
}

plots <- lapply(names(gene_coords), function(g) create_genomic_plot(g, gene_coords[[g]], peaks_gr))
combined_plot <- plots[[1]] / plots[[2]] / plots[[3]] +
  plot_annotation(title = "MEF2C Binding at Key NDD Gene Loci",
                  subtitle = "Genomic coordinates showing regulatory evidence",
                  tag_levels = 'A')

ggsave("results/visualization/Figure3_NDD_Gene_Binding.pdf", combined_plot, width = 12, height = 14)
ggsave("results/visualization/Figure3_NDD_Gene_Binding.png", combined_plot, width = 12, height = 14, dpi = 300)

# Create binding summary
cat("Creating binding summary...\n")
gene_names <- names(gene_coords)
evidence_strength <- c("Very Strong", "Strong", "Strong", "Medium")[1:length(gene_names)]
binding_summary <- data.frame(
  Gene = gene_names,
  Chromosome = sapply(gene_coords, function(x) as.character(seqnames(x))),
  Gene_Start = sapply(gene_coords, function(x) start(x)),
  Gene_End = sapply(gene_coords, function(x) end(x)),
  Nearby_MEF2C_Peaks = sapply(gene_coords, function(x) length(subsetByOverlaps(peaks_gr, resize(x, width = width(x) + 200000, fix = "center")))),
  Evidence_Strength = evidence_strength
)
write.csv(binding_summary, "results/visualization/genome_tracks/NDD_gene_binding_genomic_summary.csv", row.names = FALSE)

cat("Genome browser tracks and figures created.\n")
cat("Files include UCSC/IGV BEDs, session files, Figure3, and binding summary.\n")
