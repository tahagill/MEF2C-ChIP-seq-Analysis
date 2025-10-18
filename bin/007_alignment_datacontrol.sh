#!/bin/bash

BASE_DIR=~/pheonix/MEF2C_ChIPseq
ALIGNED_DIR="$BASE_DIR/data/processed/aligned"
RESULTS_DIR="$BASE_DIR/results/alignment_qc_new"

mkdir -p "$RESULTS_DIR"

NEW_SAMPLES=("SRR35220292" "SRR35220288" "SRR35220289" "SRR35220284")

echo "NEW SAMPLES ALIGNMENT VALIDATION REPORT" > "$RESULTS_DIR/alignment_summary_NEW.txt"
echo "Generated on $(date)" >> "$RESULTS_DIR/alignment_summary_NEW.txt"
echo "Samples: ${NEW_SAMPLES[*]}" >> "$RESULTS_DIR/alignment_summary_NEW.txt"
echo "----------------------------------------" >> "$RESULTS_DIR/alignment_summary_NEW.txt"

for SAMPLE in "${NEW_SAMPLES[@]}"; do
    bam_file="$ALIGNED_DIR/${SAMPLE}.bam"
    echo "" >> "$RESULTS_DIR/alignment_summary_NEW.txt"
    echo "$SAMPLE" >> "$RESULTS_DIR/alignment_summary_NEW.txt"
    
    if [[ ! -f "$bam_file" ]]; then
        echo "BAM file missing." >> "$RESULTS_DIR/alignment_summary_NEW.txt"
        continue
    fi
    
    if [[ ! -f "${bam_file}.bai" ]]; then
        echo "BAM index missing." >> "$RESULTS_DIR/alignment_summary_NEW.txt"
    else
        echo "BAM index OK." >> "$RESULTS_DIR/alignment_summary_NEW.txt"
    fi
    
    echo "Alignment stats:" >> "$RESULTS_DIR/alignment_summary_NEW.txt"
    grep -E "(overall alignment rate|reads paired|reads mapped|properly paired)" "$ALIGNED_DIR/${SAMPLE}_alignment_stats.txt" >> "$RESULTS_DIR/alignment_summary_NEW.txt"
    
    echo "BAM stats:" >> "$RESULTS_DIR/alignment_summary_NEW.txt"
    samtools flagstat "$bam_file" | head -5 >> "$RESULTS_DIR/alignment_summary_NEW.txt"
    
    echo "File size: $(ls -lh "$bam_file" | awk '{print $5}')" >> "$RESULTS_DIR/alignment_summary_NEW.txt"
done

echo "Generating MultiQC report for new samples..."
multiqc "$ALIGNED_DIR" -o "$RESULTS_DIR" -n "alignment_multiqc_report_NEW" \
    --filename "alignment_multiqc_report_NEW.html" \
    --ignore "*" \
    --include "${NEW_SAMPLES[@]}"

echo ""
echo "New samples validation complete."
echo "Summary: $RESULTS_DIR/alignment_summary_NEW.txt"
echo "MultiQC report: $RESULTS_DIR/alignment_multiqc_report_NEW.html"
echo "Check for:"
echo "  - All new BAM files present"
echo "  - Alignment rates above 70%"
echo "  - Properly paired reads above 90%"
echo "  - No major issues"
