#!/bin/bash


BASE_DIR=~/pheonix/MEF2C_ChIPseq
QC_DIR="$BASE_DIR/results/qc"
MULTIQC_DIR="$QC_DIR/multiqc_trimming"

mkdir -p "$MULTIQC_DIR"

echo "Running MultiQC on trimmed data FastQC reports..."
multiqc "$QC_DIR/trimmed" -o "$MULTIQC_DIR" -n "multiqc_report_trimmed"

echo "MultiQC report generated: $MULTIQC_DIR/multiqc_report_trimmed.html"
