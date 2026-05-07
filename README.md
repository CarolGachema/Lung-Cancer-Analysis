# Lung Cancer Differential Gene Expression Analysis


## Overview

This project performs a comprehensive **differential gene expression (DEG) analysis** on publicly available lung cancer microarray data, comparing gene expression profiles between **tumour tissue** and **healthy lung tissue**. The goal is to identify genes that are significantly upregulated or downregulated in lung cancer, providing molecular insights into tumourigenesis and potential therapeutic targets.

Lung cancer remains the leading cause of cancer-related mortality worldwide, accounting for approximately 1.8 million deaths annually. Understanding the transcriptomic landscape of lung tumours is critical for developing targeted diagnostics and therapies.


## Dataset

| Parameter | Details |
|-----------|---------|
| **GEO Accession** | [GSE19188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19188) |
| **Platform** | Affymetrix Human Genome U133 Plus 2.0 Array (GPL570) |
| **Comparison** | Tumour vs Healthy lung tissue |
| **Source** | NCBI Gene Expression Omnibus (GEO) |


## Methods

```
Raw GEO Data → Expression Matrix Extraction → Group Classification
→ Design Matrix → Contrast Definition → limma Linear Model
→ Empirical Bayes Statistics → FDR Correction → DEG Filtering
→ Gene Annotation → Visualisation → Results Export
```

### Statistical Framework
- **Package:** `limma` (Linear Models for Microarray Analysis)
- **Moderation:** Empirical Bayes (eBayes)
- **Multiple Testing Correction:** Benjamini-Hochberg FDR
- **Significance Thresholds:** adj.P.Val < 0.05 AND |log2FC| > 1

## Key Results

| Metric | Value |
|--------|-------|
| **Total Significant DEGs** | 3,672 |
| **Upregulated in Tumour** | 1,129 |
| **Downregulated in Tumour** | 2,543 |

### Top Differentially Expressed Genes

| Gene | Expression | Biological Relevance |
|------|-----------|---------------------|
| **TNXB** | Differentially expressed | Extracellular matrix glycoprotein; influences tumour invasion and tissue remodelling |
| **AGER** | Differentially expressed | Receptor for advanced glycation end products (RAGE); linked to lung cancer progression and inflammation |
| **ADAMTS8** | Differentially expressed | Putative tumour suppressor; implicated in ECM remodelling in cancer |
| **ADH1B** | Differentially expressed | Alcohol dehydrogenase; reflects altered metabolic landscape of tumour cells |
| **VEPH1** | Differentially expressed | Emerging evidence links it to tumour regulatory pathways |
| **LRRC36** | Differentially expressed | Leucine-rich repeat protein involved in protein-protein interactions and signalling |

> Multiple probes mapping to the same gene (e.g., AGER, ADAMTS8) indicate consistent and robust expression changes across independent probe sets — strengthening confidence in these findings.


## Visualisations

### Volcano Plot
The volcano plot below visualises effect size (log2 fold change) against statistical significance (−log10 adjusted p-value). Genes in the upper corners represent the most biologically and statistically significant DEGs. The top 10 genes are annotated directly on the plot.

 `figures/volcano_plot.png`

### Heatmap
A heatmap of the top 10 most significant DEGs across all samples, with row scaling applied to normalise expression. Tumour samples are labelled in red, healthy samples in blue.

 `figures/heatmap_top10_genes.png`


## Repository Structure

```
Lung-Cancer-DEG-Analysis/
│
├── data/
│   └── GSE19188_series_matrix.txt.gz    # Raw GEO series matrix (download from NCBI)
│
├── scripts/
│   └── lung_cancer_analysis.R           # Full annotated R analysis pipeline
│
├── results/
│   ├── all_genes_results.csv            # All genes with full statistics
│   ├── significant_genes.csv            # Filtered significant DEGs
│   ├── upregulated_genes.csv            # Genes upregulated in tumour tissue
│   └── downregulated_genes.csv          # Genes downregulated in tumour tissue
│
└── figures/
    ├── volcano_plot.png                 # Volcano plot with top 10 annotated genes
    └── heatmap_top10_genes.png          # Heatmap of top 10 DEGs
```


## How to Reproduce This Analysis

### Prerequisites
Install the following R packages before running the analysis:

```r
install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "hgu133plus2.db", "AnnotationDbi"))
```

### Steps
1. Clone this repository
2. Download the GSE19188 dataset from [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19188) and place it in the `data/` folder
3. Open `scripts/lung_cancer_analysis.R` in RStudio
4. Run the script end to end
5. Results will be saved to `results/` and `figures/`

---

## Biological Interpretation

The predominance of **downregulated genes** (2,543 vs 1,129 upregulated) is consistent with the widespread transcriptional silencing and tumour suppressor loss documented in lung cancer biology. Key findings include:

- **AGER downregulation** reflects disruption of normal pulmonary homeostasis and innate immune signalling in tumour tissue
- **ADAMTS8 downregulation** supports its role as a tumour suppressor restraining ECM remodelling in healthy lung tissue
- **ADH1B dysregulation** reflects the altered metabolic reprogramming characteristic of tumour cells (Warburg effect)
- Multiple probe sets mapping to the same genes strengthen confidence in these findings as genuine biological signals


## Limitations

- Exploratory analysis — findings require functional validation
- No clinical metadata stratification (e.g., cancer subtype, stage, treatment history)
- Microarray technology offers less precise quantification than RNA-sequencing
- Future analyses should incorporate pathway enrichment (GO, KEGG) for deeper mechanistic insight


## Author

**Caroline Gachema**


*Data sourced from NCBI GEO — publicly available for research use.*

