#!/bin/bash

BASE_DIR=~/pheonix/MEF2C_ChIPseq
FASTQ_DIR="$BASE_DIR/data/raw"
QC_DIR="$BASE_DIR/results/qc"

mkdir -p "$QC_DIR"

fastqc "$FASTQ_DIR"/*.fastq -o "$QC_DIR" -t 4

multiqc "$QC_DIR" -o "$QC_DIR"
