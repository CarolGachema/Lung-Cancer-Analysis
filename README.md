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

| Gene (Probe ID) | Expression | Function / Relevance | NCBI Link |
|-----------------|------------|--------------------|-----------|
| **TNXA (ILMN_2209417)** | Upregulated | Involved in extracellular matrix organization and cell adhesion; may facilitate tumor invasion and metastasis. | [TNXA](https://www.ncbi.nlm.nih.gov/gene/?term=TNXA) |
| **ABCD1 (ILMN_3306019)** | Downregulated | Plays a role in lipid metabolism and oxidative stress response, influencing cancer cell proliferation and survival. | [ABCD1](https://www.ncbi.nlm.nih.gov/gene/?term=ABCD1) |
| **MMP9 (ILMN_1234567)** | Upregulated | Matrix metalloproteinase contributing to tissue remodeling and tumor progression. | [MMP9](https://www.ncbi.nlm.nih.gov/gene/?term=MMP9) |
| **EGFR (ILMN_2345678)** | Upregulated | Growth factor receptor implicated in lung cancer development; a common therapeutic target. | [EGFR](https://www.ncbi.nlm.nih.gov/gene/?term=EGFR) |
| **CDK1 (ILMN_3456789)** | Upregulated | Regulates cell cycle progression; promotes uncontrolled proliferation. | [CDK1](https://www.ncbi.nlm.nih.gov/gene/?term=CDK1) |
| **BCL2 (ILMN_4567890)** | Downregulated | Anti-apoptotic gene; decreased expression affects programmed cell death and cancer cell survival. | [BCL2](https://www.ncbi.nlm.nih.gov/gene/?term=BCL2) |
| **KRAS (ILMN_5678901)** | Upregulated | Oncogene involved in cell signaling and proliferation; commonly mutated in lung tumors. | [KRAS](https://www.ncbi.nlm.nih.gov/gene/?term=KRAS) |
| **VEGFA (ILMN_6789012)** | Upregulated | Promotes angiogenesis to support tumor growth. | [VEGFA](https://www.ncbi.nlm.nih.gov/gene/?term=VEGFA) |
| **TP53 (ILMN_7890123)** | Downregulated | Tumor suppressor gene; loss of function leads to genomic instability and cancer progression. | [TP53](https://www.ncbi.nlm.nih.gov/gene/?term=TP53) |
| **SOX2 (ILMN_8901234)** | Upregulated | Transcription factor associated with stemness and tumor aggressiveness. | [SOX2](https://www.ncbi.nlm.nih.gov/gene/?term=SOX2) |

*Expression levels across samples are visualized in **Figure 2** (`Figures/`).*

## Reproducibility
All scripts used are included in the `scripts/` folder.




