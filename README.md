# Lung Cancer Differential Gene Expression Analysis

## Overview
This project analyzes gene expression data to identify differentially expressed genes in lung cancer samples compared to normal tissue.

## Dataset
GEO accession: GSE19188  
Downloaded from NCBI GEO.
The raw GSE expression data is **too large to upload** here. You can download it directly from NCBI GEO:

[GSE19188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19188)
Place the downloaded file in a folder called `data/` before running the scripts.

## Methods
- Data preprocessing in R
- Differential expression analysis using limma
- Threshold: |log2FC| > 1
- Adjusted p-value (FDR < 0.05)

## Results
- Total significant genes
- Top 10 genes identified

## Key Findings
Brief biological interpretation of major genes.

## Reproducibility
All scripts used are included in the `scripts/` folder.




