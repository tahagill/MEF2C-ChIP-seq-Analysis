#!/usr/bin/env Rscript
# DiffBind analysis with 20,000 peak set

cat("DiffBind analysis with 20,000 peaks\n")
cat("Started:", date(), "\n")

library(DiffBind)

base_dir <- "~/pheonix/MEF2C_ChIPseq"
sample_sheet <- file.path(base_dir, "metadata/sample_sheet_all.csv")
peaks_dir <- file.path(base_dir, "results/peaks")
output_dir <- file.path(base_dir, "results/differential_final_20k")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Create 20,000 peak set
cat("Creating 20,000 peak set...\n")
peaks_file <- file.path(peaks_dir, "consensus_peaks_small.bed")
all_peaks <- read.table(peaks_file, sep = "\t")
compromise_peaks <- all_peaks[1:20000, ]
write.table(compromise_peaks, 
            file.path(peaks_dir, "consensus_peaks_20k.bed"),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

cat("Using 20,000 most significant WT peaks\n")

# Read samples
samples <- read.csv(sample_sheet)
cat("Samples:", paste(samples$SampleID, collapse = ", "), "\n")

process_comparison <- function(comp_name, group1, group2) {
    cat("\n==================================================\n")
    cat("Processing:", comp_name, "\n")
    cat("==================================================\n")
    
    contrast_samples <- samples[samples$Condition %in% c(group1, group2), ]
    cat("Samples:", paste(contrast_samples$SampleID, collapse = ", "), "\n")
    
    gc(full = TRUE, verbose = FALSE)
    
    # Create DBA object
    cat("Step 1/4: Creating DBA object...\n")
    dba <- dba(sampleSheet = contrast_samples)
    
    # Count reads
    cat("Step 2/4: Counting reads in 20,000 peaks...\n")
    dba <- dba.count(dba, 
                     peaks = file.path(peaks_dir, "consensus_peaks_20k.bed"),
                     minOverlap = 2,
                     score = DBA_SCORE_READS,
                     bParallel = FALSE)
    
    cat("Peaks after counting:", nrow(dba$peaks), "\n")
    
    # Set up contrast
    cat("Step 3/4: Setting up contrast...\n")
    dba <- dba.contrast(dba, categories = DBA_CONDITION)
    
    # Analyze
    cat("Step 4/4: Running DESeq2...\n")
    dba <- dba.analyze(dba, method = DBA_DESEQ2, bParallel = FALSE)
    
    # Get results
    report <- dba.report(dba, method = DBA_DESEQ2)
    
    if (is.null(report) || length(report) == 0) {
        cat("No results generated for", comp_name, "\n")
        sig_count <- 0
        empty_df <- data.frame(Note = "No significant peaks found")
        write.csv(empty_df, file.path(output_dir, paste0(comp_name, "_full_report.csv")))
        write.csv(empty_df, file.path(output_dir, paste0(comp_name, "_significant_peaks.csv")))
    } else {
        report_df <- as.data.frame(report)
        if ("FDR" %in% colnames(report_df) && "Fold" %in% colnames(report_df)) {
            significant <- report_df[!is.na(report_df$FDR) & report_df$FDR < 0.05 & abs(report_df$Fold) > 1, ]
            sig_count <- nrow(significant)
        } else {
            cat("Statistical columns missing\n")
            sig_count <- 0
            significant <- data.frame(Note = "Statistical columns missing")
        }
        write.csv(report_df, file.path(output_dir, paste0(comp_name, "_full_report.csv")))
        write.csv(significant, file.path(output_dir, paste0(comp_name, "_significant_peaks.csv")))
    }
    
    cat(comp_name, "complete -", sig_count, "significant peaks\n")
    
    rm(dba, report, contrast_samples)
    if (exists("report_df")) rm(report_df)
    if (exists("significant")) rm(significant)
    gc(full = TRUE)
    Sys.sleep(5)
    
    return(sig_count)
}

cat("\nStarting analysis\n")
wt_ko_sig <- process_comparison("WT_vs_KO", "WT", "KO")
cat("\nMemory cleaned\n")
het_ko_sig <- process_comparison("HET_vs_KO", "HET", "KO")
cat("\nMemory cleaned\n")
wt_het_sig <- process_comparison("WT_vs_HET", "WT", "HET")

cat("\nAnalysis completed\n")
cat("WT vs KO:", wt_ko_sig, "\n")
cat("HET vs KO:", het_ko_sig, "\n")
cat("WT vs HET:", wt_het_sig, "\n")
cat("Finished:", date(), "\n")

writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
