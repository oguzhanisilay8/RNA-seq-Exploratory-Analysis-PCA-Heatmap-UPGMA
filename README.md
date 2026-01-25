# RNA-seq-Exploratory-Analysis-PCA-Heatmap-UPGMA

This repository presents an exploratory RNA-seq data analysis workflow implemented in R.  
The project focuses on understanding transcriptomic similarities and differences among samples through dimensionality reduction, clustering, and visualization techniques.

---

##  Project Overview

In this project, RNA-seq count data were analyzed to explore gene expression patterns and sample relationships using:

- Principal Component Analysis (PCA)
- Heatmap visualization of highly variable genes
- UPGMA hierarchical clustering

The aim is to demonstrate how transcriptomic data can be processed and interpreted in an exploratory bioinformatics workflow.

---

## Dataset Information

- Source: Gene Expression Omnibus (GEO, NCBI)
- Accession: GSE164073
- Data type: Raw gene count matrix
- Organism: Homo sapiens

---

##  Analysis Workflow

The main steps of the analysis include:

1. Importing raw RNA-seq count data
2. Quality control and low-count gene filtering
3. Log2 normalization
4. Data scaling
5. Principal Component Analysis (PCA)
6. Identification of highly variable genes
7. Heatmap visualization
8. UPGMA hierarchical clustering

All analyses were conducted using R and open-source packages.

## Project Structure
```text
RNA-seq-Exploratory-Analysis-PCA-Heatmap-UPGMA/
├── scripts/
│ └── analysis.R # Main analysis script
├── results/ # Output figures (PCA, heatmap, dendrogram)
├── session_info.txt # R session information for reproducibility
├── README.md
├── LICENSE
├── CITATION.cff
├── CHANGELOG.md
└── CODE_OF_CONDUCT.md
```

## How to Run


## How to Run

1. Download the dataset from GEO (GSE164073).
2. Place the count file in the project directory (or update the file path in `scripts/analysis.R`).
3. Open R / RStudio and run:

```r
source("scripts/analysis.R")
```

## AI Assistance Disclosure

Artificial intelligence tools were partially used during the development of this project, mainly for:

Improving code explanations

Refining documentation

Supporting learning and debugging processes

All data analysis, interpretation, implementation, and final validation were performed by the author.
The project reflects the author's independent understanding and applied bioinformatics skills.

## Author
Oğuzhan Işılay
Biotechnology & Bioinformatics Student
Mersin University

LinkedIn: https://www.linkedin.com/in/o%C4%9Fuzhan-i%C5%9F%C4%B1lay-3a793a348/
