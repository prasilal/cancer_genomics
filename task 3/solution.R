library(DESeq2)
library(pheatmap)


# Load count matrix
x = read.table("sample.counts", row.names=1, header=T, sep=",")
s = read.table("sample.info", header=T, row.names=1, colClasses=c("character", "factor"))

# Create DESeq2 object
dds = DESeqDataSetFromMatrix(countData = x, colData = s, design = ~ condition)

# Run a differential expression analysis (Tumour vs. Normal) using a log-fold change threshold of 1
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","Tumour","Normal"), lfcThreshold = 1)

# Generate an MA-plot
plotMA(res)

# Plot the normalized counts for the GJB2 gene
plotCounts(dds, gene="GJB2", intgroup="condition", normalized = TRUE)

# Generate a PCA plot of the samples using the transformed count data
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="condition")

rld <- rlog(dds) 
plotPCA(rld, intgroup="condition")

# Visualize the differential gene expression results as a heatmap
# Take the top 20 genes according to the adjusted p-value
top_genes <- head(order(res$padj), 20)  
pheatmap(assay(vsd)[top_genes,])
pheatmap(assay(rld)[top_genes,])

# Export the significant results (padj < 0.01) to a CSV file
s_res <- subset(res, padj < 0.01)
write.csv(as.data.frame(s_res), file="significant_results.csv")

