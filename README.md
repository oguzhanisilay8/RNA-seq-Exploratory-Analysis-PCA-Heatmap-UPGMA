# RNA-seq Exploratory Analysis  

## PCA, Heatmap and UPGMA Clustering

This repository contains an exploratory analysis of RNA-seq gene expression data using basic statistical and visualization techniques.

The goal of this project is to investigate global gene expression patterns and sample relationships using dimensionality reduction and clustering methods.

The analysis demonstrates a simple and transparent RNA-seq exploratory workflow implemented in R.

---

## Project Overview

RNA-seq experiments measure expression levels of thousands of genes across multiple samples.  
Before performing differential expression analysis, it is important to explore the structure of the dataset.

This project focuses on exploratory analysis of RNA-seq count data using:

- gene filtering
- log transformation
- Principal Component Analysis (PCA)
- heatmap visualization
- hierarchical clustering (UPGMA)

These methods help reveal global transcriptional patterns and potential sample heterogeneity.

---

## Dataset

The dataset used in this study was obtained from the GEO (Gene Expression Omnibus) database.

**Dataset ID**

GSE164073

**File used**

GSE164073_raw_counts_GRCh38.p13_NCBI.tsv.gz

**Dataset characteristics**

- organism: *Homo sapiens*
- genes: 39,376
- samples: 18
- data type: raw RNA-seq count matrix

The dataset contains gene-level expression counts for multiple human samples.
---
## Analysis Workflow

The analysis pipeline used in this repository follows the steps below:

```
RNA-seq raw count matrix
        ↓
Filtering low-expression genes
        ↓
Log2 transformation of counts
        ↓
Principal Component Analysis (PCA)
        ↓
Heatmap of highly variable genes
        ↓
Hierarchical clustering (UPGMA)
```

This workflow is commonly used to explore the structure of RNA-seq datasets before downstream statistical analysis.

---

## Methods

### Gene Filtering

Genes with extremely low counts were removed to reduce noise in the dataset.

Filtering rule used:

```
rowSums(counts) > 10
```

Only genes with sufficient total expression across samples were retained.

---

### Log Transformation

Raw RNA-seq counts often show strong variance differences.  
To stabilize variance, the following transformation was applied:

```
log2(counts + 1)
```

Adding +1 prevents log transformation errors caused by zero counts.

---

### Principal Component Analysis (PCA)

PCA was used to reduce the dimensionality of the gene expression matrix and visualize relationships between samples.

In the PCA plot:

- each point represents one sample
- distances between points reflect similarity in gene expression profiles

Samples that cluster together show similar transcriptional patterns.

---

### Heatmap of Highly Variable Genes

To highlight genes contributing most to sample differences, the **top 50 most variable genes** were selected based on variance.

A heatmap was generated to visualize expression patterns across samples.

Heatmaps help reveal:

- sample clusters
- coordinated gene expression patterns
- potential biological differences between groups.

---

### Hierarchical Clustering (UPGMA)

Hierarchical clustering was performed using the **UPGMA (average linkage)** method.

This approach groups samples based on their global gene expression similarity.

The resulting dendrogram shows how samples cluster together and provides an additional view of dataset structure.

---

## Repository Structure

```
RNA-seq-Exploratory-Analysis-PCA-Heatmap-UPGMA

│
├── figures
│   ├── pca_plot.png
│   ├── heatmap.png
│   └── upgma_dendrogram.png
│
├── rnaseq_exploratory_analysis.Rmd
├── rnaseq_exploratory_analysis.html
│
├── README.md
└── LICENSE
```

---

## How to Run the Analysis

1. Download the RNA-seq dataset from GEO

2. Place the dataset file in the project directory

```
GSE164073_raw_counts_GRCh38.p13_NCBI.tsv.gz
```

3. Open the R Markdown file

```
rnaseq_exploratory_analysis.Rmd
```

4. Run the analysis in **RStudio**

5. Knit the document to generate the HTML report.

---

## Software

Analysis was performed using **R**.

Main functions used:

- `prcomp()` for PCA
- `heatmap()` and `pheatmap()` for heatmap visualization
- `hclust()` for hierarchical clustering

---

## Reproducibility

The R session information and package versions used in the analysis are provided at the end of the R Markdown report using:

```
sessionInfo()
```

---

## Author

Oğuzhan Işılay

Biotechnology student interested in:

- bioinformatics
- computational biology
- AI applications in biology

---

## License

This project is released under the MIT License.
