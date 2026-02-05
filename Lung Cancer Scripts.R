############################################
#Lung Cancer Differential Expression Analysis
#Dataset: GSE19188
#Method: Limma
#############################################

#Loading required libraries
library(GEOquery)
library(limma)

setwd("C:/Users/Salome/Documents/MY PROJECTS/Lung Cancer Analysis")

getwd()

#Loading GEO Data (Already downloaded from NCBI GEO)
gse <- getGEO(
  filename = "data/GSE19188_series_matrix.txt.gz",
  getGPL = FALSE
)

#Extract Expression Data and Metadata
expression_data <- exprs(gse)
metadata <- pData(gse)

#Define Experimental Groups
group <- ifelse(
  grepl("healthy", metadata$characteristics_ch1, ignore.case = TRUE),
  "Healthy",
  "Tumor"
)

#Convert to factor
group <- factor(group, levels = c("Healthy", "Tumor"))

#Checking sample counts
table(group)

#Design Matrix
design <- model.matrix(~0 + group)
colnames(design) <- c("Healthy", "Tumor")
design

#Define Contrast (Tumor vs Healthy)
contrast <- makeContrasts(
  Tumor - Healthy,
  levels = design
)

#Fit linear model and apply statistics
fit <- lmFit(expression_data, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

#Extract Differential Expression Results
results <- topTable(
  fit2,
  adjust.method = "fdr",
  number = Inf
)

head(results)

#Volcano plot overview
#Shows effect size (logFC) vs significance (-log10 p-value)
plot(
  results$logFC,
  -log10(results$adj.P.Val),
  pch = 20,
  col = "grey",
  main = "Lung Cancer: Tumor vs Healthy",
  xlab = "Log2 Fold Change",
  ylab = "-log10 Adjusted P-value"
)

#Threshold lines
abline(v = c(-1, 1), col = "blue", lty = 2)
abline(h = -log10(0.05), col = "red", lty = 2)
   
#Identify significant genes
sig_genes <- results[
  results$adj.P.Val < 0.05 & abs(results$logFC) > 1,
]

nrow(sig_genes) #How many significant genes

#Annotation of genes
annotation(gse)

#Load annotation database
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
library(AnnotationDbi)

#Map probe IDs to gene symbols
gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = rownames(sig_genes),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

#Add gene symbols to results table
sig_genes$GeneSymbol <- gene_symbols

head(sig_genes[, c("GeneSymbol", "logFC", "adj.P.Val")])

#Up and downregulated genes.
up_genes <- sig_genes[sig_genes$logFC > 1, ]
down_genes <- sig_genes[sig_genes$logFC < -1, ]

nrow(up_genes)
nrow(down_genes)


#label top genes on volcano plot
top10 <- sig_genes[order(sig_genes$adj.P.Val), ][1:10, ]
points(
  top10$logFC,
  -log10(top10$adj.P.Val),
  col = "red",
  pch = 20
)

text(
  top10$logFC,
  -log10(top10$adj.P.Val),
  labels = top10$GeneSymbol,
  pos = 3,
  cex = 0.8
)

#Heatmap
top10_genes <- rownames(sig_genes[order(sig_genes$adj.P.Val), ][1:10, ])
heatmap_data <- expression_data[top20_genes, ]
heatmap(
  heatmap_data,
  scale = "row",
  ColSideColors = ifelse(group == "Tumor", "red", "blue"),
  main = "Top 10 Differentially Expressed Genes"
)
#Saving heatmap
png("figures/heatmap_top10_genes.png", width = 900, height = 700)
dev.off()

# Save results 
write.csv(results, "results/all_genes_results.csv")
write.csv(sig_genes, "results/significant_genes.csv")
write.csv(up_genes, "results/upregulated_genes.csv")
write.csv(down_genes, "results/downregulated_genes.csv")
png("figures/volcano_plot.png", width = 800, height = 600)
dev.off()

load("results/lung_cancer_workspace.RData")



