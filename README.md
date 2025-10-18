# MEF2C ChIP-seq Analysis in Human Microglia

Computational pipeline for analyzing MEF2C transcription factor binding patterns in iPSC-derived microglia. This repository contains all code and documentation from my undergraduate thesis investigating cell-type specific MEF2C regulatory programs and their relevance to neurodevelopmental disorders.

## Background

MEF2C haploinsufficiency causes severe neurodevelopmental delay, epilepsy, and autistic features. While MEF2C is well-characterized in neurons and cardiac tissue, its transcriptional targets in microglia—the brain's resident immune cells critical for synaptic pruning—remain largely unknown. This analysis identifies MEF2C genomic binding sites in human microglia and examines their functional relevance.

**Key Results:**
- 1,258 high-confidence MEF2C binding sites identified through differential analysis
- Enrichment for synaptic process genes (vesicle transport, synaptic vesicle cycle)
- Direct binding near neurodevelopmental disorder genes (ANK3, NRXN1, ARID1B, DSCAM)
- Substantial cell-type specificity with patterns following developmental lineage
- Network analysis identifies ANK3 as central hub in 110-gene synaptic module

## Repository Organization
```
MEF2C-ChIP-seq-Analysis/
│
├── scripts/                          # Sequential analysis pipeline (001-027)
│   ├── 001_download_samples.sh       # SRA data acquisition
│   ├── 002_fastqc.sh                 # Initial quality control
│   ├── 003_trimgalore_rawdata.sh     # Adapter trimming
│   ├── 004_fastqc_trimmed.sh         # Post-trim QC
│   ├── 006_alignment_bwt.sh          # Bowtie2 alignment to GRCh38
│   ├── 009_macs2_narrow_peak_calling_pooled.sh  # Peak calling
│   ├── 011_create_consensus_peaks.sh # Consensus peak generation
│   ├── 013_run_diffBind_small_peaks.R  # Differential binding (DiffBind)
│   ├── 014_peak_annotation.R         # ChIPseeker genomic annotation
│   ├── 015_peak_annotation_report_fixed.R  # Annotation summary
│   ├── 017_pathway_analysis.R        # GO/KEGG enrichment (clusterProfiler)
│   ├── 019_motif_analysis_meme.sh    # De novo motif discovery (MEME)
│   ├── 020_NDD_Enrichment.R          # Neurodevelopmental disorder gene overlap
│   ├── 021_download_encode_data.sh   # ENCODE comparison datasets
│   ├── 022_public_data_comparison.R  # Cross-cell-type analysis
│   ├── 023a_ppi_network_analysis.R   # STRING PPI network
│   ├── 023b_clean_network_visualisation.R  # Network visualization
│   ├── 025_genome_browser_tracks.R   # Track generation
│   ├── 026_integration_figures.R     # Multi-panel figures
│   └── 027_genome_browser_track.R    # Final browser tracks
│
├── environments/
│   ├── chipseq_diffbind.yml          # Main conda environment
│   └── meme_env.yml                  # Motif analysis environment
│
├── config/
│   └── sample_metadata.csv           # Sample annotations and SRA accessions
│
└── README.md                          # This file
```

## Data Source

**Public ChIP-seq Data:**
- **Repository:** NCBI Sequence Read Archive (SRA)
- **Accessions:** SRR35220282-SRR35220292
- **Samples:** 9 total (3 biological replicates × 3 genotypes)
  - Wild-type MEF2C+/+ (n=3)
  - Heterozygous MEF2C+/− (n=3)
  - Knockout MEF2C−/− (n=3)
- **Cell Type:** Human iPSC-derived microglia (isogenic lines)
- **Target:** MEF2C transcription factor

**Comparative Data:**
- ENCODE MEF2C ChIP-seq datasets from cardiomyocytes, skeletal myotubes, monocytes, lymphoblastoid cells, and HepG2 liver cells

## System Requirements

### Hardware
- **RAM:** 18GB minimum (differential binding analysis optimized for this constraint)
- **Storage:** ~100GB for raw data, alignments, and outputs
- **OS:** Ubuntu Linux (tested on 20.04 LTS)
- **CPU:** Multi-core recommended for alignment (8+ cores ideal)

### Software Dependencies

**Core Analysis Environment** (`chipseq_diffbind.yml`):
```yaml
- sra-tools=3.0.0          # Data download
- fastqc=0.11.9            # Quality control
- multiqc=1.11             # Aggregate QC reports
- trimgalore=0.6.7         # Adapter trimming
- bowtie2=2.4.5            # Read alignment
- samtools=1.15            # BAM manipulation
- macs2=2.2.7.1            # Peak calling
- bedtools=2.30.0          # Genomic interval operations
- R=4.2.0
  - DiffBind              # Differential binding analysis
  - ChIPseeker            # Peak annotation
  - clusterProfiler       # Pathway enrichment
  - org.Hs.eg.db          # Human genome annotation
  - STRINGdb              # Protein interaction networks
  - igraph                # Network analysis
  - ggplot2               # Visualization
  - ComplexHeatmap        # Heatmaps
```

**Motif Analysis Environment** (`meme_env.yml`):
```yaml
- meme=5.4.1              # Motif discovery
```

## Installation

### 1. Clone Repository
```bash
git clone https://github.com/tahagill/MEF2C-ChIP-seq-Analysis.git
cd MEF2C-ChIP-seq-Analysis
```

### 2. Install Conda/Mamba
If not already installed:
```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
```

### 3. Create Analysis Environments
```bash
# Main analysis environment
mamba env create -f environments/chipseq_diffbind.yml

# Motif analysis environment
mamba env create -f environments/meme_env.yml
```

### 4. Download Reference Genome
```bash
# Create reference directory
mkdir -p reference

# Download GRCh38 Bowtie2 index (no alt contigs)
cd reference
wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
unzip GRCh38_noalt_as.zip
cd ..
```

## Usage

### Automated Pipeline Execution
```bash
# Activate environment
conda activate chipseq_diffbind

# Execute full pipeline (creates all output directories automatically)
bash scripts/run_all.sh

# Expected runtime: 24-36 hours depending on hardware
```

### Manual Step-by-Step Execution

#### Phase 1: Data Acquisition and Quality Control
```bash
conda activate chipseq_diffbind

# Download raw data from SRA
bash scripts/001_download_samples.sh

# Initial quality assessment
bash scripts/002_fastqc.sh

# Adapter trimming and quality filtering
bash scripts/003_trimgalore_rawdata.sh

# Post-trimming QC
bash scripts/004_fastqc_trimmed.sh
```

#### Phase 2: Alignment
```bash
# Align to GRCh38 with Bowtie2
bash scripts/006_alignment_bwt.sh

# Output: Sorted, indexed BAM files in alignments/
```

#### Phase 3: Peak Calling
```bash
# Call peaks with MACS2 (pooled replicates per genotype)
bash scripts/009_macs2_narrow_peak_calling_pooled.sh

# Generate consensus peak set
bash scripts/011_create_consensus_peaks.sh
```

#### Phase 4: Differential Binding Analysis
```bash
# Identify genotype-dependent binding sites
# Note: Uses top 20,000 WT peaks due to RAM constraints
Rscript scripts/013_run_diffBind_small_peaks.R

# Output: 1,258 significant differential peaks (FDR < 0.05)
```

#### Phase 5: Genomic Annotation
```bash
# Annotate peaks to genomic features
Rscript scripts/014_peak_annotation.R

# Generate annotation report
Rscript scripts/015_peak_annotation_report_fixed.R

# Output: 1,148 genes with proximal MEF2C binding
```

#### Phase 6: Functional Enrichment
```bash
# GO and KEGG pathway analysis
Rscript scripts/017_pathway_analysis.R

# Output: Enriched synaptic processes, vesicle transport pathways
```

#### Phase 7: Motif Discovery
```bash
# Switch to motif environment
conda activate meme_env

# De novo motif discovery with MEME
bash scripts/019_motif_analysis_meme.sh

# Output: 10 significant motifs including canonical MEF2C sequence
conda deactivate
```

#### Phase 8: Disease Integration
```bash
conda activate chipseq_diffbind

# Test enrichment for NDD gene sets
Rscript scripts/020_NDD_Enrichment.R

# Download ENCODE data for comparison
bash scripts/021_download_encode_data.sh

# Cross-cell-type specificity analysis
Rscript scripts/022_public_data_comparison.R
```

#### Phase 9: Network Analysis
```bash
# Construct PPI network from STRING database
Rscript scripts/023a_ppi_network_analysis.R

# Network visualization and community detection
Rscript scripts/023b_clean_network_visualisation.R

# Output: 719-node network with 22 communities
```

#### Phase 10: Final Visualizations
```bash
# Generate genome browser tracks
Rscript scripts/025_genome_browser_tracks.R
Rscript scripts/027_genome_browser_track.R

# Create integrated multi-panel figures
Rscript scripts/026_integration_figures.R
```

## Key Outputs
```
results/
├── peaks/
│   ├── WT_pooled_peaks.narrowPeak              # 106,199 peaks
│   ├── HET_pooled_peaks.narrowPeak             # 73,018 peaks  
│   ├── KO_pooled_peaks.narrowPeak              # 14,848 peaks
│   └── differential_binding/
│       └── significant_peaks_WT_vs_KO.bed      # 1,258 differential peaks
│
├── annotation/
│   ├── annotated_peaks.csv                     # Genomic features
│   ├── peak_annotation_summary.txt             # Distribution stats
│   └── target_genes.txt                        # 1,148 genes
│
├── enrichment/
│   ├── GO_BP_enrichment.csv                    # Biological processes
│   ├── KEGG_pathway_enrichment.csv             # KEGG pathways
│   └── NDD_gene_overlap.csv                    # Disease gene enrichment
│
├── motifs/
│   └── meme_output/
│       ├── meme.html                           # Motif report
│       └── meme.txt                            # Motif sequences
│
├── networks/
│   ├── PPI_network.graphml                     # Full network
│   ├── hub_genes.csv                           # Top 78 hubs
│   └── communities.csv                         # 22 functional modules
│
├── cell_type_comparison/
│   └── overlap_statistics.csv                  # Cross-tissue analysis
│
└── figures/
    ├── Figure1_MEF2C_Regulatory_Network.png
    ├── Figure2_CellType_Specificity.pdf
    ├── Figure3_NDD_Gene_Binding.pdf
    ├── GO_enrichment_dotplot.pdf
    ├── PPI_network_visualization.pdf
    └── genome_browser_tracks.bed
```

## Analysis Parameters

### Peak Calling (MACS2)
- **Mode:** Narrow peaks (transcription factor binding)
- **q-value cutoff:** 0.05
- **Genome size:** hs (human)
- **Fragment length:** Auto-detected from data
- **Duplicate handling:** Keep all (0% duplication rate observed)

### Differential Binding (DiffBind)
- **FDR threshold:** 0.05 (Benjamini-Hochberg correction)
- **Comparison:** Wild-type vs Knockout
- **Peak subset:** Top 20,000 most significant WT peaks (memory optimization)
- **Normalization:** TMM (edgeR)
- **Minimum overlap:** 2 samples required for consensus peaks

### Functional Enrichment (clusterProfiler)
- **Background:** All human genes
- **p-value cutoff:** 0.05
- **q-value cutoff:** 0.05 (FDR correction)
- **Minimum gene set size:** 10 genes
- **Maximum gene set size:** 500 genes

### Network Analysis (STRINGdb)
- **Confidence score:** > 400 (medium-high confidence)
- **Species:** Homo sapiens (9606)
- **Network type:** Physical + functional interactions
- **Community detection:** Louvain algorithm

### Motif Discovery (MEME)
- **Mode:** ANR (any number of repetitions)
- **Number of motifs:** 10
- **Width range:** 6-20 bp
- **Search:** Both strands
- **E-value threshold:** Auto (reported up to 1e-3)

## Quality Control Metrics

### Sequencing Statistics
- **Total reads:** 313M raw → 311M trimmed (99.2% retention)
- **Median read length:** 101 bp → 96-98 bp post-trimming
- **GC content:** 40-43% (expected for human)
- **Adapter contamination:** < 0.1% after trimming

### Alignment Statistics
- **Alignment rates:** 91.7-97.9% across samples
- **PCR duplicates:** 0.0% (excellent library complexity)
- **Read distribution:** Balanced across genotypes
  - KO: 107.9M reads
  - HET: 108.5M reads  
  - WT: 94.5M reads

### Peak Quality
- **Peak score range:** 3.1 to 15,420 (MACS2)
- **Enrichment:** Up to 25-fold over background
- **Reproducibility:** High concordance between replicates
- **Motif validation:** Canonical MEF2C motif identified (E=2.4e-217)

## Computational Considerations

### Memory Optimization
Due to 18GB RAM constraints, differential binding analysis was performed on a subset of peaks:
- **Strategy:** Top 20,000 most significant wild-type peaks
- **Rationale:** Focuses on highest-confidence binding sites
- **Trade-off:** May miss weaker but biologically relevant sites
- **Validation:** Results consistent with full peak set trends

This conservative approach ensures:
1. Computationally feasible analysis on personal hardware
2. Focus on high-confidence MEF2C binding events
3. Robust statistical power for differential analysis
4. Reproducibility on similar systems

### Runtime Estimates
On Ubuntu 20.04 with 18GB RAM and 8 CPU cores:
- Data download: 2-3 hours
- QC and trimming: 1-2 hours
- Alignment: 4-6 hours
- Peak calling: 1-2 hours
- Differential analysis: 2-3 hours
- Annotation and enrichment: 1 hour
- Motif discovery: 2-3 hours
- Network analysis: 30 minutes
- Visualization: 1 hour
- **Total:** ~24-36 hours

## Interpretation Guidelines

**Important:** This analysis identifies MEF2C genomic binding sites through ChIP-seq but does not prove transcriptional regulation. Key considerations:

1. **ChIP-seq identifies binding, not regulation**
   - Binding sites are candidate regulatory elements
   - RNA-seq validation required to confirm gene expression changes
   - Not all bound sites are functionally active

2. **Cross-study comparisons have technical limitations**
   - Peak calling parameters vary between studies
   - Batch effects and experimental protocols differ
   - Biological patterns (e.g., lineage gradients) more informative than absolute overlap percentages

3. **Pathway enrichments indicate associations**
   - Enriched terms suggest functional relevance
   - Multiple testing correction applied (FDR < 0.05)
   - Convergent evidence strengthens confidence

4. **Network predictions are computational**
   - Based on literature and experimental databases
   - Hub genes warrant experimental validation
   - Interaction confidence scores provided

**Validation experiments needed:**
- RNA-seq from MEF2C-deficient microglia
- Enhancer reporter assays for bound regions
- ChIP-qPCR confirmation of select binding sites
- Functional assays (e.g., microglial synaptic pruning capacity)

## Troubleshooting

### Out of Memory Errors
```bash
# If DiffBind fails with memory error:
# Edit script 013: Reduce peak count from 20,000 to 15,000 or 10,000
# Line ~25: top_peaks <- ranked_peaks[1:15000,]
```

### MACS2 Installation Issues
```bash
# If pip install fails:
conda install -c bioconda macs2=2.2.7.1

# Or use mamba (faster):
mamba install -c bioconda macs2
```

### Bowtie2 Index Download Problems
```bash
# Alternative source:
wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/GRCh38_noalt_as.zip

# Or build from scratch (slower):
bowtie2-build GRCh38.fa GRCh38_noalt_as
```

### R Package Errors
```bash
# In R console:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DiffBind", "ChIPseeker", "clusterProfiler", 
                       "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene"))
```

### STRING Database Connection
```R
# If STRINGdb fails:
# Check internet connection
# Or download pre-built database:
# https://string-db.org/cgi/download
```

## Citation

If you use this analysis pipeline or adapt the approach, please cite:
```
Ahmad, T. (2024). ChIP-seq Analysis Defines Microglia-Specific MEF2C Binding 
Landscapes at Synaptic Genes Implicated in Neurodevelopmental Disorders. 
Undergraduate Thesis, Middle East Technical University, Department of Biological Sciences.
GitHub: https://github.com/tahagill/MEF2C-ChIP-seq-Analysis
```

## Reproducibility

This analysis is fully reproducible using the provided:
- Sequential numbered scripts (001-027)
- Conda environment specifications
- Public data accessions
- Documented parameters

Expected results should match within normal computational variance (< 5% for count-based metrics).

## License

**Code:** MIT License - free to use, modify, and distribute with attribution

**Data:** Original ChIP-seq data is publicly available through NCBI SRA. ENCODE data subject to ENCODE data use policies.

## Acknowledgments

- Original data generators for depositing ChIP-seq datasets to public repositories
- ENCODE Consortium for providing comparative MEF2C ChIP-seq data
- Bioconductor community for ChIP-seq analysis tools
- Open-source bioinformatics software developers

## Contact

**Author:** Taha Ahmad  
**Institution:** Middle East Technical University  
**Department:** Biological Sciences  
**Project Type:** Undergraduate Thesis (Independent Computational Research)  
**GitHub:** [@tahagill](https://github.com/tahagill)

---

**Note:** This represents independent undergraduate computational research performed entirely on personal computing resources. All data sources are publicly available. This work is not affiliated with or endorsed by any institution and was conducted as an independent academic project.

**Last Updated:** October 2024
