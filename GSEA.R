library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(fgsea)
library(msigdbr)
library(ggplot2)

res_gsea <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
ranked_genes <- res_gsea$log2FoldChange
names(ranked_genes) <- rownames(res_gsea)
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

#Load Hallmark Gene Sets
msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H" )
hallmark_gene_sets <-split(msig_hallmark$gene_symbol, msig_hallmark$gs_name)

common_genes <- intersect(names(ranked_genes), unlist(hallmark_gene_sets))
ranked_genes <- ranked_genes[common_genes]

set.seed(42)

fgsea_results <- fgsea(pathways = hallmark_gene_sets,
                       stats = ranked_genes,
                       minSize = 15,
                       maxSize = 500)

fgsea_results <- fgsea_results %>% arrange(padj)
head(fgsea_results)

top_pathways <- fgsea_results [fgsea_results$padj < 0.05, ]$pathway[1:10]
ggplot(fgsea_results %>% filter(pathway %in% top_pathways),
       aes(reorder(pathway, NES), NES, fill = padj < 0.05 )) +
  geom_col() +
  coord_flip() +
  labs(title = "Top Enriched Pathways", x = "Pathway", y = "Normalized Enrichement Score") +
  theme_minimal()

#Specific Pathway
head(top_pathways)

plotEnrichment(hallmark_gene_sets[["HALLMARK_KRAS_SIGNALING_UP"]], ranked_genes) +
  labs(title = "HALLMARK_KRAS_SIGNALING_UP")
####

ridgeplot(fgsea_results, showCategory = 10 ) + ggtitle("GSEA Ridge Plot")

####
library(pheatmap)
# Create a list to hold leading edge genes
leading_edge_genes <- list()

# Loop through each result in gsea_results
for (i in 1:nrow(fgsea_results@result)) {
  # Extract the core enrichment gene ids
  genes <- strsplit(fgsea_results@result$core_enrichment[i], "/")[[1]]
  
  # Store them in the list using the pathway ID as the name
  leading_edge_genes[[fgsea_results_df@result$ID[i]]] <- genes
}

# Remove any NULL or empty entries
leading_edge_genes <- leading_edge_genes[sapply(leading_edge_genes, function(x) length(x) > 0)]


print(leading_edge_genes)

pheatmap(do.call(rbind, lapply(leading_edge_genes, function(x) ranked_genes[x])))

# Prepare heatmap data
heatmap_data <- do.call(rbind, lapply(leading_edge_genes, function(x) {
  if (length(x) > 0 && all(x %in% names(ranked_genes))) {
    return(ranked_genes[x])
  } else {
    return(NULL)
  }
}))

# Clean heatmap data by removing NA, NaN, or Inf values
heatmap_data <- heatmap_data[complete.cases(heatmap_data), ]

# Generate the heatmap if data is valid
if (nrow(heatmap_data) > 0) {
  pheatmap(heatmap_data)
} else {
  warning("No valid data available for the heatmap after cleaning.")
}



