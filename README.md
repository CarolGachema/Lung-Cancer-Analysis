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
- Top 10 genes identified (in results folder)

## Key Findings
# Top 10 Differentially Expressed Genes in Lung Cancer

The top 10 significantly differentially expressed genes are summarized in **Table 1** (see `Figures/` folder). Below are brief descriptions with links to [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/) pages for each gene.

## Top Differentially Expressed Genes in Lung Cancer

The following genes ranked highest based on differential expression analysis.  
Some genes appeared multiple times due to multiple probes mapping to the same gene, indicating consistent expression changes.

| Gene Symbol | Expression* | Function / Relevance | NCBI Gene Link |
|-------------|------------|----------------------|----------------|
| **TNXB** | Differentially expressed | Encodes Tenascin-X, an extracellular matrix glycoprotein involved in matrix organization and cell adhesion. Alterations may influence tumor invasion and tissue remodeling. | https://www.ncbi.nlm.nih.gov/gene/?term=TNXB |
| **TNXA** | Differentially expressed (multiple probes) | Pseudogene related to Tenascin family members. While not protein-coding, its expression may reflect genomic or regulatory alterations in tumor tissue. | https://www.ncbi.nlm.nih.gov/gene/?term=TNXA |
| **AGER** | Differentially expressed (multiple probes) | Encodes the receptor for advanced glycation end products (RAGE). Highly expressed in lung tissue and involved in inflammation and immune signaling; dysregulation is linked to lung cancer progression. | https://www.ncbi.nlm.nih.gov/gene/?term=AGER |
| **ADH1B** | Differentially expressed | Alcohol dehydrogenase involved in ethanol metabolism and oxidative stress pathways. Altered metabolic gene expression is commonly observed in tumor cells. | https://www.ncbi.nlm.nih.gov/gene/?term=ADH1B |
| **ADAMTS8** | Differentially expressed (multiple probes) | Member of the ADAMTS metalloproteinase family. Often reported as a tumor suppressor gene and implicated in extracellular matrix remodeling in cancers. | https://www.ncbi.nlm.nih.gov/gene/?term=ADAMTS8 |
| **VEPH1** | Differentially expressed | Encodes Ventricular Zone Expressed PH Domain Containing 1. Suggested roles in cell signaling and development; emerging evidence links it to tumor regulatory pathways. | https://www.ncbi.nlm.nih.gov/gene/?term=VEPH1 |
| **LRRC36** | Differentially expressed | Leucine-rich repeat-containing protein; limited characterization, but leucine-rich repeat proteins are often involved in protein–protein interactions and signaling pathways. | https://www.ncbi.nlm.nih.gov/gene/?term=LRRC36 |

## Reproducibility
All scripts used are included in the `scripts/` folder.




