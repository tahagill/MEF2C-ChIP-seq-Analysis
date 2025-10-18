#!/usr/bin/env Rscript

library(igraph)
library(ggplot2)
library(dplyr)
library(ggraph)
library(tidygraph)
library(patchwork)

# Load network data
cat("Loading network data...\n")
hub_genes <- read.csv("results/network_analysis/PPI_hub_genes.csv")
communities <- read.csv("results/network_analysis/PPI_communities.csv")

# Load original target genes
cat("Loading target genes...\n")
target_genes <- read.csv("results/annotations/MEF2C_target_genes_FINAL.csv")
genes <- unique(target_genes$SYMBOL)

# Create focused subnetwork: hubs + ANK3 neighborhood + NDD genes
cat("Creating focused subnetworks...\n")
top_hubs <- head(hub_genes$Gene, 20)
ndd_genes <- c("ANK3", "DSCAM", "ARID1B", "CACNB2", "NRXN1", "NLGN1", "SH3GL2")

# Visualization 1: Hub Gene Bar Plot
cat("Creating hub gene bar plot...\n")
top_20_hubs <- head(hub_genes, 20)
p1 <- ggplot(top_20_hubs, aes(x=reorder(Gene, Degree), y=Degree, fill=Degree)) +
  geom_bar(stat="identity") +
  coord_flip() +
  labs(title="Top 20 Hub Genes in MEF2C Network",
       subtitle="Degree Centrality (Number of Connections)",
       x="Gene", y="Degree") +
  scale_fill_gradient(low="blue", high="red") +
  theme_minimal() +
  theme(legend.position="none")
ggsave("results/network_analysis/hub_gene_barplot.pdf", p1, width=10, height=8)

# Visualization 2: Community Size Plot
cat("Creating community size visualization...\n")
p2 <- ggplot(communities, aes(x=reorder(as.factor(Community), Size), y=Size, fill=Size)) +
  geom_bar(stat="identity") +
  labs(title="MEF2C Network Functional Modules",
       subtitle="Community Sizes from Louvain Clustering", 
       x="Community ID", y="Number of Genes") +
  scale_fill_gradient(low="lightblue", high="darkblue") +
  theme_minimal()
ggsave("results/network_analysis/community_sizes.pdf", p2, width=10, height=6)

# Visualization 3: ANK3-focused table
cat("Creating ANK3 network analysis...\n")
ank3_info <- hub_genes[hub_genes$Gene == "ANK3", ]
community_7_genes <- hub_genes[hub_genes$Community == 7, ] %>% arrange(desc(Degree))
write.csv(community_7_genes, "results/network_analysis/ANK3_community_genes.csv", row.names=FALSE)

# Visualization 4: Network Properties Plot
cat("Creating network properties visualization...\n")
network_stats <- data.frame(
  Metric = c("Nodes", "Edges", "Hub Genes", "Communities", "Network Density"),
  Value = c(719, 4340, 78, nrow(communities), 0.0168)
)
p3 <- ggplot(network_stats, aes(x=Metric, y=Value, fill=Metric)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=Value), vjust=-0.5, size=5) +
  labs(title="MEF2C PPI Network Statistics", y="Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave("results/network_analysis/network_statistics.pdf", p3, width=10, height=6)

# Visualization 5: Degree Distribution
cat("Creating degree distribution plot...\n")
p4 <- ggplot(hub_genes, aes(x=Degree)) +
  geom_histogram(binwidth=5, fill="steelblue", alpha=0.8) +
  geom_vline(xintercept=quantile(hub_genes$Degree, 0.9), linetype="dashed", color="red") +
  annotate("text", x=quantile(hub_genes$Degree, 0.9)+10, y=50, 
           label="Hub Threshold (90%)", color="red") +
  labs(title="Degree Distribution of MEF2C Target Genes",
       subtitle="Most genes have few connections, hubs have many",
       x="Degree (Number of Connections)", y="Number of Genes") +
  theme_minimal()
ggsave("results/network_analysis/degree_distribution.pdf", p4, width=10, height=6)

# Multi-panel summary figure
cat("Creating multi-panel summary figure...\n")
summary_plot <- (p1 | p3) / (p2 | p4) +
  plot_annotation(title="MEF2C Target Gene PPI Network Analysis",
                  subtitle="Comprehensive Network Characterization",
                  theme=theme(plot.title=element_text(size=16, face="bold")))
ggsave("results/network_analysis/network_analysis_summary.pdf", summary_plot, width=16, height=12)

# ANK3-specific analysis
cat("Creating ANK3-specific analysis...\n")
ank3_analysis <- list(
  ANK3_Network_Properties = ank3_info,
  ANK3_Community_Top_Genes = head(community_7_genes, 15),
  ANK3_Significance = "ANK3 is a major network hub (Degree=30) in the synaptic function community"
)
write.csv(ank3_analysis$ANK3_Community_Top_Genes, 
          "results/network_analysis/ANK3_community_top_genes.csv", row.names=FALSE)

cat("Network analysis complete. Generated clean, publication-ready visualizations:\n")
cat("- hub_gene_barplot.pdf\n")
cat("- community_sizes.pdf\n")
cat("- network_statistics.pdf\n")
cat("- degree_distribution.pdf\n")
cat("- network_analysis_summary.pdf\n")
cat("- ANK3_community_genes.csv\n")
