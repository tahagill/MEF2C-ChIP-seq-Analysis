\BASE_DIR=~/pheonix/MEF2C_ChIPseq
BIN_DIR="$BASE_DIR/bin"
LOG_DIR="$BASE_DIR/logs/auto_pipeline"

mkdir -p "$LOG_DIR"

SCRIPTS=("002_fastqc.sh" "003_trimgalore_rawdata.sh" "004_fastqc_trimmed.sh" "006_alignment_bwt.sh" "007_alignment_datacontrol.sh")

for script in "${SCRIPTS[@]}"; do
    if [[ -f "$BIN_DIR/$script" ]]; then
        chmod +x "$BIN_DIR/$script"
        cd "$BIN_DIR"

        echo "Running $script..."
        "./$script" 2>&1 | tee "$LOG_DIR/${script}.log"

        if [[ ${PIPESTATUS[0]} -eq 0 ]]; then
            echo "$script completed successfully."
        else
            echo "$script failed. Check log: $LOG_DIR/${script}.log"
            echo "Stopping pipeline."
            exit 1
        fi
    else
        echo "Could not find script: $BIN_DIR/$script"
        exit 1
    fi
done

echo ""
echo "All scripts finished successfully on $(date)"
echo "Logs are saved in $LOG_DIR"
echo "Next step: run MACS2 for peak calling."
