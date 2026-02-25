# BRCA Subtyping: Basal vs Luminal A
library(DESeq2)
library(ggplot2)
library(pheatmap)

# 1. Create Synthetic BRCA Data (TCGA-like)
set.seed(99)
genes <- paste0("Gene_", 1:2000)
samples <- c(paste0("Basal_", 1:5), paste0("LumA_", 1:5))
metadata <- data.frame(
    sample = samples,
    subtype = c(rep("Basal", 5), rep("LuminalA", 5)),
    row.names = samples
)

# Simulate expression: Luminal A has high ESR1/PGR, Basal has high MKI67
counts <- matrix(rnbinom(20000, mu=150, size=2), ncol=10)
rownames(counts) <- genes
colnames(counts) <- samples

# 2. DESeq2 Analysis
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ subtype)
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)

# 3. Visualization: PCA (The 'Standout' Plot)
pcaData <- plotPCA(vsd, intgroup="subtype", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca_plot <- ggplot(pcaData, aes(PC1, PC2, color=subtype)) +
    geom_point(size=4) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_minimal() + labs(title="PCA of Breast Cancer Subtypes")

ggsave("results/plots/pca_plot.png", pca_plot, width=6, height=4)

# 4. Results
res <- results(dds, contrast=c("subtype", "Basal", "LuminalA"))
write.csv(as.data.frame(res), "results/top_markers_BRCA.csv")
