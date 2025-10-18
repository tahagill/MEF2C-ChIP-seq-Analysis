#!/usr/bin/env Rscript

cat("Installing DiffBind with error handling\n")

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(Ncpus = 4)

is_installed <- function(pkg) {
    requireNamespace(pkg, quietly = TRUE)
}

# Install BiocManager if not present
if (!is_installed("BiocManager")) {
    cat("Installing BiocManager...\n")
    install.packages("BiocManager")
} else {
    cat("BiocManager already installed\n")
}

library(BiocManager)

# Install DiffBind if needed
if (is_installed("DiffBind")) {
    cat("DiffBind already installed\n")
    library(DiffBind)
    cat("Version:", as.character(packageVersion("DiffBind")), "\n")
} else {
    cat("Installing DiffBind and dependencies...\n")
    cat("This may take 15-30 minutes.\n\n")
    
    tryCatch({
        BiocManager::install("DiffBind", ask = FALSE, update = FALSE, force = TRUE)
        
        if (is_installed("DiffBind")) {
            library(DiffBind)
            cat("DiffBind installation successful\n")
            cat("Version:", as.character(packageVersion("DiffBind")), "\n")
        } else {
            cat("DiffBind installation failed\n")
            cat("Possible fixes:\n")
            cat("1. Install system dependencies:\n")
            cat("   sudo apt install libxml2-dev libcurl4-openssl-dev libssl-dev\n")
            cat("2. Run this script again\n")
            quit(status = 1)
        }
    }, error = function(e) {
        cat("Error during installation:\n")
        cat(conditionMessage(e), "\n")
        cat("Install system libraries first:\n")
        cat("sudo apt install libxml2-dev libcurl4-openssl-dev libssl-dev zlib1g-dev libbz2-dev liblzma-dev\n")
        quit(status = 1)
    })
}

cat("DiffBind installation complete\n")
