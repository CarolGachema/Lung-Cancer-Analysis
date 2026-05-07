############################################
# Lung Cancer Differential Gene Expression Analysis
# Dataset: GSE19188 
# Comparison: Tumor vs Healthy Lung Tissue
# Method: Limma (Linear Models for Microarray Analysis)
#############################################



#Load Required Libraries
library(GEOquery)  # For downloading and loading GEO datasets from NCBI
library(limma)     # For differential expression analysis using linear models


# Set and Verify Working Directory
setwd("C:/Users/Caroline/Documents/MY PROJECTS/Lung Cancer Analysis")


# Load GEO Dataset.

# Load the pre downloaded series matrix file from NCBI GEO
gse <- getGEO(
  filename = "data/GSE19188_series_matrix.txt.gz",
  getGPL = FALSE
)


# Extract Expression Data and Sample Metadata
expression_data <- exprs(gse) # Gene expression matrix (rows = probes, cols = samples)
metadata <- pData(gse)        # Sample metadata containing group labels and clinical info


# Define Experimental Groups

# Classify each sample as either "Healthy" or "Tumor"
# based on the characteristics_ch1 field in the sample metadata
# Samples containing "healthy" are control tissue; all others are tumour samples
group <- ifelse(
  grepl("healthy", metadata$characteristics_ch1, ignore.case = TRUE),
  "Healthy",
  "Tumor"
)


# Convert group to a factor with Healthy as the reference level
# This ensures the contrast is correctly interpreted as Tumor vs Healthy
group <- factor(group, levels = c("Healthy", "Tumor"))


# Verify sample counts per group to confirm correct classification
table(group)


# Build the Design Matrix

# Create a no-intercept design matrix for two-group comparison
# Each column represents one experimental group
design <- model.matrix(~0 + group)
colnames(design) <- c("Healthy", "Tumor")
design


# Define the Contrast of Interest

# Define Tumor minus Healthy as our contrast
# This identifies genes differentially expressed in tumour
# tissue compared to healthy lung tissue
contrast <- makeContrasts(
  Tumor - Healthy,
  levels = design
)


# Fit Linear Model and Apply Empirical Bayes Statistics
fit <- lmFit(expression_data, design)   # Fit linear model to the expression data
fit2 <- contrasts.fit(fit, contrast)    # Apply the Tumor vs Healthy contrast
fit2 <- eBayes(fit2)                    # Apply empirical Bayes moderation for
                                        # more robust and stable variance estimates

# Extract Differential Expression Results

# Retrieve statistics for all genes with FDR-adjusted p-values
# number = Inf ensures all genes are returned, not just the top subset
results <- topTable(
  fit2,
  adjust.method = "fdr",
  number = Inf
)


# Preview the top differentially expressed genes
head(results)


# Generate Volcano Plot

# Volcano plot visualises both effect size (logFC on x-axis) and
# statistical significance (-log10 adjusted p-value on y-axis)
# Genes in the upper left/right corners are the most biologically significant
plot(
  results$logFC,
  -log10(results$adj.P.Val),
  pch = 20,
  col = "grey",
  main = "Lung Cancer: Tumor vs Healthy",
  xlab = "Log2 Fold Change",
  ylab = "-log10 Adjusted P-value"
)


# Add threshold lines to the volcano plot:
# Vertical blue dashed lines = fold change thresholds (logFC = ±1, i.e. 2-fold change)
# Horizontal red dashed line = significance threshold (adj.P.Val = 0.05)
abline(v = c(-1, 1), col = "blue", lty = 2)
abline(h = -log10(0.05), col = "red", lty = 2)
 
  
# Identify significant genes

# Apply dual significance thresholds to identify biologically meaningful DEGs:
# FDR-adjusted p-value < 0.05 (statistically significant after correction)
# |log2 fold change| > 1 (at least 2-fold difference in expression)
sig_genes <- results[
  results$adj.P.Val < 0.05 & abs(results$logFC) > 1,
]


nrow(sig_genes) # Report total number of significant DEGs identified


#Annotation of genes

# Check GEO annotation associated with the dataset
annotation(gse)


# Install and load the Affymetrix HG-U133 Plus 2.0 annotation package
# This package maps Affymetrix probe IDs to HGNC gene symbols and other identifiers
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
library(AnnotationDbi)


# Map probe IDs to gene symbols

# Translate Affymetrix probe IDs into human-readable HGNC gene symbols
# Keys are the probe IDs from the significant genes results table
# multiVals = "first" retains only the first gene symbol when multiple
# genes map to a single probe
gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = rownames(sig_genes),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)


# Add gene symbols as a new column in the significant genes table
sig_genes$GeneSymbol <- gene_symbols


# Preview the top annotated significant genes with key statistics
head(sig_genes[, c("GeneSymbol", "logFC", "adj.P.Val")])


# Upregulated and downregulated genes.

# Upregulated genes: higher expression in Tumor compared to Healthy
# logFC > 1 means at least 2-fold higher expression in tumour tissue
up_genes <- sig_genes[sig_genes$logFC > 1, ]


# Downregulated genes: lower expression in Tumor compared to Healthy
# logFC < -1 means at least 2-fold lower expression in tumour tissue
down_genes <- sig_genes[sig_genes$logFC < -1, ]


# Report counts of upregulated and downregulated genes
nrow(up_genes)
nrow(down_genes)


# Annotate Top Genes on Volcano Plot

# Highlight and label the top 10 most significant DEGs on the volcano plot
# Ordered by adjusted p-value (most significant first)
top10 <- sig_genes[order(sig_genes$adj.P.Val), ][1:10, ]
points(
  top10$logFC,
  -log10(top10$adj.P.Val),
  col = "red",
  pch = 20
)


# Add gene symbol labels above each highlighted point
text(
  top10$logFC,
  -log10(top10$adj.P.Val),
  labels = top10$GeneSymbol,
  pos = 3,
  cex = 0.8
)


png("figures/volcanoplot_top_10_genes.png", width = 800, height = 600)
dev.off


# Generate Heatmap of Top 10 DEGs

# Extract expression values for the top 10 most significant genes
# across all samples for visualisation in a heatmap
top10_genes <- rownames(sig_genes[order(sig_genes$adj.P.Val), ][1:10, ])
heatmap_data <- expression_data[top10_genes, ]


# Generate heatmap with:
# Row scaling (scale = "row") to normalise expression across samples
# Red column labels for Tumor samples, blue for Healthy samples
heatmap(
  heatmap_data,
  scale = "row",
  ColSideColors = ifelse(group == "Tumor", "red", "blue"),
  main = "Top 10 Differentially Expressed Genes"
)


# Save heatmap as PNG to figures folder
png("figures/heatmap_top10_genes.png", width = 900, height = 700)
dev.off()


# Save All Results to CSV
write.csv(results, "results/all_genes_results.csv")
write.csv(sig_genes, "results/significant_genes.csv")
write.csv(up_genes, "results/upregulated_genes.csv")
write.csv(down_genes, "results/downregulated_genes.csv")


# Save volcano plot as PNG to figures folder
png("figures/volcano_plot.png", width = 800, height = 600)
dev.off()


# Save and Reload Workspace
# Useful for resuming analysis without re-running the full pipeline
load("results/lung_cancer_workspace.RData")



