#Libraray packages
library(dplyr)
library(ggplot2)
library(stringr)
library(reshape2)
library(pheatmap)
library(DESeq2)
library(readr)
library(GenomicFeatures)

########################################
# Determine the number of detected genes in each sample
gene_counts <- read.table('pathway/gene.counts.tsv', row.names = 1, header = TRUE)
head(gene_counts)

# Filter out genes with all zero counts using dplyr
gene_counts <- gene_counts[rowSums(gene_counts) > 0, ]
dim(gene_counts)

detected_genes <- colSums(gene_counts > 0)
print(detected_genes)

#Histogram

# Sample data in a data frame
data <- data.frame(Sample_ID = c("","","",...),
                   Numbers = c(detected_genes_results))

ggplot(data, aes(x = Sample_ID, y = Numbers)) +
  geom_bar(stat = "identity", fill = "lightblue", color = "black") +
  geom_hline(yintercept = detected_genes, color = "red", linetype = "dashed", linewidth = 1.2) +  # Horizontal line at y = the number of detected genes
  geom_text(aes(label = Numbers), vjust = -0.5, color = "black") +  # Numbers above bars
  annotate("text", x = 1, y = detected_genes + 100, label = "of the number detected_genes", color = "red", vjust = -0.5) +  # Number above constant line
  labs(title = "Detected Genes per Sample ID", x = "Sample ID", y = "Number of Detected Genes") +
  theme_minimal()

################################################
#Expression Levels
#Checking the geneID in both files to make them fit in the same names
#Example: geneIDs start with  "X..." 

gene_annotation <- read_table('gene_annotation.tsv')

gene_counts$gene_id <- rownames(gene_counts)

merged_data <- merge(gene_annotation, gene_counts, by = "gene_id" )
head(merged_data)

merged_data$gene_name[is.na(merged_data$gene_name)] <- "Unknown_Gene"
merged_data$gene_name <- make.unique(as.character(merged_data$gene_name))

rownames(merged_data) <- merged_data$gene_name
merged_data <- merged_data[, -c(1,2)]
head(merged_data)

top_genes <- head(order(apply(merged_data, 1, var), decreasing = TRUE), 50)
head(top_genes)

long_counts<- reshape2::melt(merged_data, variable.name = "Sample", value.name = "Expression")
head(long_counts)

#Boxplot
ggplot(long_counts, aes(x = Sample, y = Expression)) +
  geom_boxplot(fill = "skyblue", alpha = 0.7) +
  labs(title = "Gene Expression Distribution", 
       x = "Sample ID", 
       y = "Expression Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Violin Plot
ggplot(long_counts, aes(x = Sample, y = Expression, fill = Sample)) +
  geom_violin(alpha = 0.7) +
  labs(title = "Gene Expression Distribution", 
       x = "Sample ID", 
       y = "Expression Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#Heatmap

pheatmap(merged_data[top_genes, ],
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE)


################################################
#Coverage Per Gene

txdb <- txdbmaker::makeTxDbFromGFF("/SFS/home/fuyib/genes.gtf", format ="gtf")
gene_lengths <- exonsBy(txdb, by = "gene")
gene_lengths <- sapply(gene_lengths, function(x) sum(width(x)))

gene_lengths_df <- data.frame(gene_id = names(gene_lengths), gene_length = gene_lengths)
gene_counts <- merge(gene_counts, gene_lengths_df, by = "gene_id", all.x = TRUE)

gene_counts$Coverage <- rowSums(gene_counts[, -c(1, ncol(gene_counts))]) / gene_counts$gene_length
head(gene_counts$gene_length)



#Histogram
ggplot(gene_counts, aes(x = Coverage)) + 
  geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7, color = "black") +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "Coverage Per Gene",
       x = "Coverage (Read Counts / Gene Length)",
       y = "Number of Genes")


#Density Plot
ggplot(gene_counts, aes(x = Coverage)) + 
  geom_density(fill = "blue", alpha = 0.5) +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "Coverage Per Gene",
       x = "Coverage (Read Counts / Gene Length)",
       y = "Density")

################################
#Read Distribution Across Genes
gene_sums <- rowSums(merged_data)
gene_data <- data.frame(Gene = rownames(merged_data), Reads = gene_sums)

summary(gene_data$Reads)
quantile(gene_data$Reads, probs = seq(0 ,1, 0.1))

#Histogram
ggplot(gene_data, aes(x = Reads)) +
  geom_histogram(binwidth = 0.2, fill = "blue", alpha = 0.7, color = "black") +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "Read Distribution Across Genes",
       x = "Total Reads per Gene",
       y = "Number of Genes")

#Density Plot
ggplot(gene_data, aes(x = Reads)) +
  geom_density(fill = "blue", alpha = 0.5) +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "Read Distribution Across Genes",
       x = "Total Reads per Gene",
       y = "Density")

#Boxplot
ggplot(gene_data, aes(y = Reads)) +
  geom_boxplot(fill = "skyblue", alpha = 0.7) +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Read Distribution Across Genes",
       y = "Total Reads per Gene (Log10 Scale")

