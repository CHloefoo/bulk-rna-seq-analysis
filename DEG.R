library(DESeq2)
library(tidyverse)
library(ggplot2)

gene_counts <- read.table('/SFS/home/fuyib/PID22957/gene.counts.tsv', row.names = "gene_id", header = TRUE)
gene_counts <- gene_counts[, sort(colnames(gene_counts))]
head(gene_counts)

colSums(gene_counts)

# Add Annotation
gene_annotation <- read_table('gene_annotation.tsv')

gene_counts$gene_id <- rownames(gene_counts)

merged_data <- merge(gene_annotation, gene_counts, by = "gene_id" )
head(merged_data)

merged_data$gene_name[is.na(merged_data$gene_name)] <- "Unknown_Gene"
merged_data$gene_name <- make.unique(as.character(merged_data$gene_name))

rownames(merged_data) <- merged_data$gene_name
merged_data <- merged_data[, -c(1,2)]

all(colnames(merged_data) == rownames(metadata))

# Load the metadata
metadata <- read_excel("metadata.xlsx")
metadata$ShortID <- substr(metadata$'Sample ID', nchar(metadata$'Sample ID') - 4 , nchar(metadata$'Sample ID'))
metadata$ShortID <- paste0("X", metadata$ShortID)

metadata$NumbericPart <- as.numeric(gsub("[^0-9]", "", metadata$ShortID))
metadata <- metadata[order(metadata$NumbericPart), ]

metadata <- metadata[metadata$ShortID != "X968D1", ]
metadata <- metadata[!duplicated(metadata$ShortID), ]

metadata <- metadata %>% column_to_rownames(var = "ShortID")

all(colnames(merged_data) == rownames(metadata))


metadata$ReplicateGroup <- as.factor(metadata$ReplicateGroup)

#DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = merged_data,
                              colData = metadata,
                              design = ~ ReplicateGroup)

dds <- dds[rowSums(counts(dds)) > 10, ]

dds <- DESeq(dds)

dds

#Normalize
res <- counts(dds, normalized = TRUE)

res <- results(dds, contrast = c("ReplicateGroup","Diseased","Normal"))

res <- res[order(res$padj), ]

head(res)


#
df <-  rownames_to_column(as.data.frame(res), var = "gene_id")
head(df)

# Volcano plot of results

df$Significance <- ifelse(df$padj < 0.05 & abs(df$log2FoldChange) > 1, "Significant", "Not Significant")
df$label <- ifelse(df$Significance == "Significant", df$gene_id, "")


ggplot(df, aes(x=log2FoldChange, y=-log10(padj), color = Significance, label = label)) +
  geom_point(alpha=0.6) +
  geom_hline(yintercept=-log10(0.05), linetype='dashed', color='red') +
  ggrepel::geom_text_repel(max.overlaps = 15 , size = 3 )+
  labs(title='Volcano Plot of Differentially Expressed Genes',
       x='Log2 Fold Change',
       y='-Log10 P-value') +
  theme_minimal()


#########

#Heatmap Top 50 GEDs
topGenes <- head(rownames(res[order(res$padj), ]), 50)

vsd <- vst(dds, blind = FALSE)

head(vsd)

df <- assay(vsd)[topGenes, ]

pheatmap(df, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         scale = "row")

###################
pcaData <- plotPCA(vsd, intgroup = "ReplicateGroup", returnData = TRUE)

pcaData$Sample <- rownames(pcaData)

ggplot(pcaData, aes(PC1, PC2, color = ReplicateGroup, label = Sample))+
  geom_point(size = 4) +
  geom_text(vjust = 1.5, size = 4) +
  theme_minimal() +
  labs(title = "PCA of Gene Expression Data", x = "PC1", y = "PC2")+
  theme(plot.title = element_text(hjust = 0.5))

