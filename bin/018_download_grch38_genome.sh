#!/bin/bash
# 018_download_grch38_genome.sh
# Download and prepare GRCh38 primary assembly (no ALT contigs)

echo "=== GRCh38 GENOME DOWNLOAD AND PREPARATION ==="

GENOME_DIR="resources/genomes/hg38"
mkdir -p "$GENOME_DIR"
cd "$GENOME_DIR" || exit

echo "Downloading GRCh38 primary assembly from Ensembl..."
wget -c http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

echo "Extracting genome FASTA..."
gunzip -f Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

echo "Renaming file to standard name..."
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38_noalt_as.fa

echo "Indexing genome for samtools/bedtools..."
samtools faidx GRCh38_noalt_as.fa

echo "Creating chromosome sizes file..."
cut -f1,2 GRCh38_noalt_as.fa.fai > chrom.sizes

echo "Download and preparation complete."
echo "Genome FASTA: $GENOME_DIR/GRCh38_noalt_as.fa"
echo "Chromosome sizes: $GENOME_DIR/chrom.sizes"
echo "File size: $(du -h GRCh38_noalt_as.fa | cut -f1)"
echo "Number of chromosomes: $(grep -c '^>' GRCh38_noalt_as.fa)"
