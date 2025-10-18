#!/bin/bash
echo "ENCODE MEF2C Data Download (4 Verified Cell Types)"

mkdir -p results/public_data/comparison
cd results/public_data/comparison

echo "Downloading MEF2C datasets from 4 cell types..."



echo "1. Downloading MEF2C in cardiomyocytes..."
wget "https://www.encodeproject.org/files/ENCFF002CJQ/@@download/ENCFF002CJQ.bed.gz" -O MEF2C_cardio.bed.gz


echo "2. Downloading MEF2C in lymphoblastoid cells..."
wget "https://www.encodeproject.org/files/ENCFF002CJW/@@download/ENCFF002CJW.bed.gz" -O MEF2C_lymph.bed.gz


echo "3. Downloading MEF2C in skeletal muscle myotubes..."
wget "https://www.encodeproject.org/files/ENCFF002CJP/@@download/ENCFF002CJP.bed.gz" -O MEF2C_myotube.bed.gz


echo "4. Downloading MEF2C in monocytes..."
wget "https://www.encodeproject.org/files/ENCFF002CJR/@@download/ENCFF002CJR.bed.gz" -O MEF2C_monocyte.bed.gz


echo "5. Downloading MEF2C in HepG2 liver cells (backup)..."
wget "https://www.encodeproject.org/files/ENCFF002CJS/@@download/ENCFF002CJS.bed.gz" -O MEF2C_hepg2.bed.gz

echo "Checking downloads..."
success_count=0
for file in *.bed.gz; do
    if [ -s "$file" ]; then
        echo "✓ $file downloaded successfully"
        ((success_count++))
    else
        echo "✗ $file failed - removing"
        rm "$file" 2>/dev/null
    fi
done

echo "Unzipping successful downloads..."
gunzip -f *.bed.gz 2>/dev/null

echo "Final files downloaded: $success_count/5"
ls -lh *.bed 2>/dev/null || echo "No BED files downloaded"
