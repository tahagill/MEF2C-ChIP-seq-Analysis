#!/bin/bash

BASE_DIR=~/pheonix/MEF2C_ChIPseq
ALIGNED_DIR="$BASE_DIR/data/processed/aligned"
RESULTS_DIR="$BASE_DIR/results/alignment_qc"

mkdir -p "$RESULTS_DIR"

echo "ALIGNMENT VALIDATION REPORT" > "$RESULTS_DIR/alignment_summary.txt"
echo "Generated on $(date)" >> "$RESULTS_DIR/alignment_summary.txt"
echo "----------------------------------------" >> "$RESULTS_DIR/alignment_summary.txt"

for bam_file in "$ALIGNED_DIR"/*.bam; do
    sample=$(basename "$bam_file" .bam)
    echo "" >> "$RESULTS_DIR/alignment_summary.txt"
    echo "$sample" >> "$RESULTS_DIR/alignment_summary.txt"
    
    if [[ ! -f "$bam_file" ]]; then
        echo "BAM file missing." >> "$RESULTS_DIR/alignment_summary.txt"
        continue
    fi
    
    if [[ ! -f "${bam_file}.bai" ]]; then
        echo "BAM index missing." >> "$RESULTS_DIR/alignment_summary.txt"
    else
        echo "BAM index found." >> "$RESULTS_DIR/alignment_summary.txt"
    fi
    
    echo "Alignment stats:" >> "$RESULTS_DIR/alignment_summary.txt"
    grep -E "(overall alignment rate|reads paired|reads mapped|properly paired)" "$ALIGNED_DIR/${sample}_alignment_stats.txt" >> "$RESULTS_DIR/alignment_summary.txt"
    
    echo "BAM stats:" >> "$RESULTS_DIR/alignment_summary.txt"
    samtools flagstat "$bam_file" | head -5 >> "$RESULTS_DIR/alignment_summary.txt"
    
    echo "File size: $(ls -lh "$bam_file" | awk '{print $5}')" >> "$RESULTS_DIR/alignment_summary.txt"
done

echo "Running MultiQC..."
multiqc "$ALIGNED_DIR" -o "$RESULTS_DIR" -n "alignment_multiqc_report"

echo ""
echo "Validation complete."
echo "Summary: $RESULTS_DIR/alignment_summary.txt"
echo "MultiQC: $RESULTS_DIR/alignment_multiqc_report.html"
echo "Review for:"
echo "  - BAM files present"
echo "  - Alignment rates above 70%"
echo "  - Properly paired reads above 90%"
echo "  - No major issues detected"
