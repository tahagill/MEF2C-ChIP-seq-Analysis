#!/bin/bash

BASE_DIR=~/pheonix/MEF2C_ChIPseq
BOWTIE2_DIR="$BASE_DIR/resources/genomes/hg38/bowtie2"

mkdir -p "$BOWTIE2_DIR"
cd "$BOWTIE2_DIR" || exit

echo "Downloading pre-built Bowtie2 hg38 index from alternative source..."
wget -c https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip

echo "Unzipping..."
unzip -o GRCh38_noalt_as.zip

echo "Verifying files..."
ls -la *.bt2

echo "Pre-built index ready at: $BOWTIE2_DIR/"
