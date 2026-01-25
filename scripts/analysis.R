---
title: "RNA-seq Data Analysis: Sample Comparison Using PCA, Heatmap, and UPGMA Clustering"
output:
  pdf_document: default
  html_document: default
date: "2026-01-25"
---
# Introduction

In this study, raw RNA-seq count data obtained from the GEO database were analyzed to investigate gene expression patterns.

The analysis workflow includes:

Filtering low-expression genes,

Log transformation of count data,

Principal Component Analysis (PCA) to explore sample variability,

Heatmap visualization of highly variable genes,

Hierarchical clustering using the UPGMA method.

The main objective of this study is to demonstrate basic exploratory RNA-seq analysis steps using a manual and transparent workflow.

#  Dataset Description 
The following RNA-seq dataset was used in this study:

Source: GEO (Gene Expression Omnibus)

File: GSE164073_raw_counts_GRCh38.p13_NCBI.tsv.gz

Content: Gene-level raw count matrix

Dimension: 39,376 genes × 18 samples

The dataset consists of human (Homo sapiens) samples.

#  Data Preprocessing 

```{r}

# File name
file <- "GSE164073_raw_counts_GRCh38.p13_NCBI.tsv.gz"  #Downloaded from the GEO database


# Read the file
counts_tbl <- read.delim(gzfile(file), check.names = FALSE)

# Check dimensions
dim(counts_tbl)




```





#   Principal Component Analysis (PCA)

```{r}
# Extract gene IDs from the first column
gene_id <- counts_tbl[[1]]

# Convert count data (excluding first column) into a numeric matrix
counts <- as.matrix(counts_tbl[, -1])

# Assign gene IDs as row names of the matrix
rownames(counts) <- gene_id


# Check the dimensions of the count matrix (genes × samples)
dim(counts)

# Display the first 3 genes and first 3 samples
counts[1:3, 1:3]


# Filter out low-expression genes
# Keep genes with total counts greater than 10 across all samples
keep <- rowSums(counts) > 10

# Create a filtered count matrix
counts_f <- counts[keep, , drop = FALSE]


# Check dimensions after filtering
dim(counts_f)


# Apply log2 transformation to stabilize variance
# Adding +1 prevents log(0) errors
expr_log <- log2(counts_f + 1)


# Scale and transpose the data for PCA
# Samples become rows, genes become columns
expr_scaled <- scale(t(expr_log))


# Perform Principal Component Analysis (PCA)
pca <- prcomp(expr_scaled)


# Plot the first two principal components
plot(pca$x[,1], pca$x[,2],
     pch = 19,                         # Solid circle points
     xlab = "PC1",                     # X-axis label
     ylab = "PC2",                     # Y-axis label
     main = "PCA (log2(counts+1), scaled)")  # Plot title


# Add sample labels to each point
text(pca$x[,1], pca$x[,2],
     labels = colnames(counts_f),      # Sample names
     pos = 3,                          # Position above points
     cex = 0.6)                        # Text size



```


How to Read This PCA Plot

This PCA plot summarizes the global gene expression patterns of all samples in two main components (PC1 and PC2), which capture the largest sources of variation in the dataset.

Each dot represents one RNA-seq sample

The distance between dots reflects similarity in gene expression

Samples that are close together have similar transcriptomic profiles

Samples that are far apart show distinct expression patterns

This PCA plot shows how samples cluster based on global gene expression patterns, revealing underlying biological variation and sample heterogeneity.


## Enhanced Visualization of PCA Results

```{r}
# Convert PCA scores into a data frame
# Each row represents one sample
# PC1 and PC2 are the first two principal components
pca_df <- data.frame(
  sample = rownames(pca$x),   # Sample names
  PC1 = pca$x[,1],            # Scores for Principal Component 1
  PC2 = pca$x[,2]             # Scores for Principal Component 2
)


# Create an improved PCA scatter plot using base R
plot(pca_df$PC1, pca_df$PC2,
     pch = 19,                # Solid circle points
     cex = 1.3,               # Point size
     col = "steelblue",       # Point color
     xlab = "PC1",            # X-axis label
     ylab = "PC2",            # Y-axis label
     main = "PCA (log2(counts+1), scaled)")  # Plot title


# Add sample labels next to each point
text(pca_df$PC1, pca_df$PC2,
     labels = pca_df$sample,  # Sample names
     pos = 3,                 # Position above points
     cex = 0.7)               # Label size


# Add a background grid for better readability
grid()

```






#  Heatmap
```{r}
# Calculate variance for each gene across all samples
# This measures how much each gene's expression changes between samples
gene_var <- apply(expr_log, 1, var)


# Select the top 50 genes with the highest variance
# These genes show the strongest expression differences
top <- order(gene_var, decreasing = TRUE)[1:50]


# Create a heatmap of the selected genes
heatmap(
  
  # Extract expression values of the top variable genes
  # Transpose and scale them for visualization
  scale(t(expr_log[top, ])),
  
  scale = "none",                 # Disable additional scaling by heatmap()
  main = "Top 50 Variable Genes"  # Title of the heatmap
)


```

How to Read This Heatmap 

This heatmap shows the expression patterns of the top 50 most variable genes across all RNA-seq samples.

Rows represent samples

Columns represent genes

Colors indicate relative expression levels

Dark red: High expression

Light/yellow: Low expression

The dendrograms (tree structures) display hierarchical clustering based on gene expression similarity.

The heatmap reveals distinct expression patterns among samples based on the most variable genes. Samples cluster into several groups, suggesting shared transcriptional profiles. Certain gene clusters show coordinated upregulation and downregulation, indicating potential biological differences between sample groups.


## Heatmap Visualization of Highly Variable Genes

```{r}
# Install and load the pheatmap package for advanced heatmap visualization
install.packages("pheatmap")
library(pheatmap)


# Extract expression values of the top 50 most variable genes
# Rows = genes, Columns = samples
mat <- expr_log[top, ]


# Apply Z-score normalization to each gene
# This standardizes expression values across samples
mat_s <- t(scale(t(mat)))


# Generate a basic heatmap without sample annotation
pheatmap(mat_s)


# Extract sample names from the filtered count matrix
samples <- colnames(counts_f)


# Create temporary (dummy) group labels for visualization
# These labels are used only to add color annotations
group <- rep(c("Group1","Group2"), length.out = length(samples))
names(group) <- samples


# Convert group labels into a factor
cols <- as.factor(group)


# Create an annotation data frame for the heatmap
# This will add colored bars above the samples
ann <- data.frame(Group = cols)
rownames(ann) <- colnames(mat_s)


# Generate an annotated heatmap with sample group information
pheatmap(mat_s,
         annotation_col = ann,        # Add group annotation
         show_colnames = TRUE,        # Display sample names
         show_rownames = FALSE,       # Hide gene names for clarity
         main = "Top 50 Variable Genes")  # Heatmap title


```

Difference Between the Two Heatmaps

First heatmap:

Shows only gene expression values.

No information about sample groups.

Just helps you see general patterns.

Second heatmap (this one):

Includes group information (Group1 / Group2 at the top).

Has colored annotations.

Uses standardized values (z-score).

Looks more organized and scientific.

Why the Second One Is Better

This heatmap shows:

Which samples belong to which group.

How samples cluster together.

Whether groups are separated by gene expression.


#  UPGMA


```{r}
# Calculate the distance matrix between samples
# Each distance represents how different two samples are
# based on their gene expression profiles
dist_matrix <- dist(t(expr_log))


# Perform hierarchical clustering using the UPGMA method
# "average" refers to average linkage clustering
hc <- hclust(dist_matrix, method = "average")


# Plot the hierarchical clustering dendrogram
plot(hc,
     main = "UPGMA Clustering of RNA-seq Samples",
     xlab = "Samples",
     ylab = "Height")


```


How to Read This UPGMA Clustering Plot

This dendrogram shows the hierarchical clustering of RNA-seq samples based on their gene expression profiles using the UPGMA (average linkage) method.

Each label at the bottom represents one sample.

Samples that join together at lower heights are more similar.

Samples that merge at higher heights are more different.

The vertical height indicates the distance between clusters.

What This Analysis Shows

This clustering confirms the patterns observed in PCA and heatmap analyses. Samples that are close in PCA space and show similar heatmap profiles also tend to cluster together here, supporting the consistency of the results.



