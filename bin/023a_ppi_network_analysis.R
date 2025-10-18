#!/usr/bin/env Rscript



library(STRINGdb)
library(igraph)
library(ggplot2)
library(dplyr)

cat("PHASE 2a: PPI Network Analysis of MEF2C Targets\n")

# Load MEF2C target genes
cat("1. Loading MEF2C target genes\n")
target_genes <- read.csv("results/annotations/MEF2C_target_genes_FINAL.csv")
genes <- unique(target_genes$SYMBOL)
cat("Total MEF2C target genes:", length(genes), "\n")

# Initialize STRING database
cat("2. Connecting to STRING database\n")
string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=400)
cat("STRING database initialized\n")

# Map genes to STRING IDs
cat("3. Mapping genes to STRING IDs\n")
mapped_genes <- string_db$map(data.frame(gene=genes), "gene", removeUnmappedRows=TRUE)
cat("Mapped genes:", nrow(mapped_genes), "\n")

# Retrieve protein-protein interactions
cat("4. Retrieving protein-protein interactions\n")
interactions <- string_db$get_interactions(mapped_genes$STRING_id)
cat("High-confidence interactions:", nrow(interactions), "\n")

# Build network
cat("5. Building PPI network\n")
network <- graph_from_data_frame(interactions[, c("from", "to")], directed=FALSE)

# Assign gene names to vertices
vertex_names <- setNames(mapped_genes$gene, mapped_genes$STRING_id)
V(network)$name <- vertex_names[V(network)$name]
V(network)$label <- V(network)$name

# Compute network metrics
cat("6. Calculating network metrics\n")
V(network)$degree <- degree(network)
V(network)$betweenness <- betweenness(network)
V(network)$eigen_centrality <- eigen_centrality(network)$vector

# Identify hub genes (top 10% by degree)
degree_threshold <- quantile(V(network)$degree, 0.9)
hub_genes <- V(network)$name[V(network)$degree >= degree_threshold]

cat("7. Network statistics\n")
cat("Nodes:", vcount(network), "\n")
cat("Edges:", ecount(network), "\n")
cat("Density:", round(graph.density(network), 4), "\n")
cat("Hub genes (top 10%):", length(hub_genes), "\n")
cat("Top hub genes:", paste(head(sort(hub_genes, decreasing=TRUE), 5), collapse=", "), "\n")

# Identify clusters/communities
cat("8. Identifying functional modules\n")
communities <- cluster_louvain(network)
V(network)$community <- communities$membership

# Create results directory
dir.create("results/network_analysis", showWarnings=FALSE)

# Save network results
cat("9. Saving results\n")

# Hub gene table
hub_gene_table <- data.frame(
  Gene = V(network)$name,
  Degree = V(network)$degree,
  Betweenness = V(network)$betweenness,
  Eigen_Centrality = V(network)$eigen_centrality,
  Community = V(network)$community
) %>% arrange(desc(Degree))

write.csv(hub_gene_table, "results/network_analysis/PPI_hub_genes.csv", row.names=FALSE)

# Community summary
community_summary <- hub_gene_table %>%
  group_by(Community) %>%
  summarise(
    Size = n(),
    Top_Genes = paste(head(Gene, 3), collapse=", "),
    Mean_Degree = mean(Degree)
  ) %>% arrange(desc(Size))

write.csv(community_summary, "results/network_analysis/PPI_communities.csv", row.names=FALSE)

# Basic network visualization
cat("10. Creating network visualization\n")
png("results/network_analysis/PPI_network_basic.png", width=1200, height=1000)
set.seed(123)
plot(communities, network,
     vertex.size = log(V(network)$degree + 1) * 2,
     vertex.color = V(network)$community,
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     main = "MEF2C Target Gene PPI Network",
     sub = paste("Nodes:", vcount(network), "Edges:", ecount(network)))
dev.off()

# ANK3-specific analysis
cat("11. ANK3-focused network analysis\n")
if("ANK3" %in% V(network)$name) {
  ank3_neighbors <- neighbors(network, V(network)[name == "ANK3"])
  cat("ANK3 degree:", degree(network, V(network)[name == "ANK3"]), "\n")
  cat("ANK3 neighbors:", length(ank3_neighbors), "\n")
  cat("Neighbor genes:", paste(ank3_neighbors$name, collapse=", "), "\n")
}

# Summary report
cat("12. Generating summary report\n")
sink("results/network_analysis/PPI_analysis_summary.txt")
cat("MEF2C TARGET GENE PPI NETWORK ANALYSIS\n")
cat("=======================================\n\n")
cat("Input genes:", length(genes), "\n")
cat("Genes mapped to STRING:", nrow(mapped_genes), "\n")
cat("High-confidence interactions:", nrow(interactions), "\n")
cat("Network nodes:", vcount(network), "\n")
cat("Network edges:", ecount(network), "\n")
cat("Network density:", round(graph.density(network), 4), "\n\n")

cat("TOP HUB GENES:\n")
print(head(hub_gene_table, 10))

cat("\nFUNCTIONAL MODULES (Communities):\n")
print(community_summary)

cat("\nANK3 NETWORK PROPERTIES:\n")
if("ANK3" %in% V(network)$name) {
  ank3_info <- hub_gene_table[hub_gene_table$Gene == "ANK3", ]
  cat("Degree:", ank3_info$Degree, "\n")
  cat("Betweenness:", round(ank3_info$Betweenness, 2), "\n")
  cat("Community:", ank3_info$Community, "\n")
}
sink()


cat("Results saved to: results/network_analysis/\n")
cat("Key files:\n")
cat("- PPI_hub_genes.csv\n")
cat("- PPI_communities.csv\n")
cat("- PPI_network_basic.png\n")
cat("- PPI_analysis_summary.txt\n")
