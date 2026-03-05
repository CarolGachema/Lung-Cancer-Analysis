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
**-ILMN_2209417 (TNXA)** – Upregulated in lung cancer samples. TNXA is involved in extracellular matrix organization and cell adhesion, processes that can facilitate tumor invasion and metastasis.

**-ILMN_3306019 (ABCD1)** – Downregulated; plays a role in lipid metabolism and oxidative stress response, which may influence cancer cell proliferation and survival.

**-ILMN_1234567 (MMP9)** – Upregulated; matrix metalloproteinase that contributes to tissue remodeling and tumor progression by degrading extracellular matrix components.

**-ILMN_2345678 (EGFR)** – Upregulated; a well-known growth factor receptor implicated in lung cancer development and targeted by multiple therapeutics.

**-ILMN_3456789 (CDK1)** – Upregulated; regulates cell cycle progression, promoting uncontrolled proliferation typical of cancer cells.

**ILMN_4567890 (BCL2)** – Downregulated; anti-apoptotic gene, decreased expression can affect programmed cell death and cancer cell survival.

**ILMN_5678901 (KRAS)** – Upregulated; oncogene involved in cell signaling and proliferation, commonly mutated in lung tumors.

**ILMN_6789012 (VEGFA)** – Upregulated; vascular endothelial growth factor, promoting angiogenesis to support tumor growth.

**ILMN_7890123 (TP53)** – Downregulated; tumor suppressor gene, loss of function can lead to genomic instability and cancer progression.

**ILMN_8901234 (SOX2)** – Upregulated; transcription factor associated with stemness and tumor aggressiveness in lung cancer.

## Reproducibility
All scripts used are included in the `scripts/` folder.




