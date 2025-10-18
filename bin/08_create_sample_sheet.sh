#!/bin/bash

BASE_DIR=~/pheonix/MEF2C_ChIPseq
ALIGNED_DIR="$BASE_DIR/data/processed/aligned"
SAMPLESHEET_DIR="$BASE_DIR/metadata"

mkdir -p "$SAMPLESHEET_DIR"

# Create sample sheet for diffBind
cat > "$SAMPLESHEET_DIR/sample_sheet_all.csv" << EOF
SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,Peaks,PeakCaller
SRR35220290,microglia,MEF2C,WT,1,$ALIGNED_DIR/SRR35220290.bam,,,
SRR35220291,microglia,MEF2C,WT,2,$ALIGNED_DIR/SRR35220291.bam,,,
SRR35220292,microglia,MEF2C,WT,3,$ALIGNED_DIR/SRR35220292.bam,,,
SRR35220287,microglia,MEF2C,HET,1,$ALIGNED_DIR/SRR35220287.bam,,,
SRR35220288,microglia,MEF2C,HET,2,$ALIGNED_DIR/SRR35220288.bam,,,
SRR35220289,microglia,MEF2C,HET,3,$ALIGNED_DIR/SRR35220289.bam,,,
SRR35220282,microglia,MEF2C,KO,1,$ALIGNED_DIR/SRR35220282.bam,,,
SRR35220283,microglia,MEF2C,KO,2,$ALIGNED_DIR/SRR35220283.bam,,,
SRR35220284,microglia,MEF2C,KO,3,$ALIGNED_DIR/SRR35220284.bam,,,
EOF

echo "Sample sheet created: $SAMPLESHEET_DIR/sample_sheet_all.csv"
echo "Total samples: 9 (3 WT + 3 HET + 3 KO)"
