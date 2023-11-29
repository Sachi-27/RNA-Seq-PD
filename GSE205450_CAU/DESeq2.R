library("DESeq2")

countData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
colData <- read.csv("metadata.csv", header=TRUE, sep=",", row.names=1)

print(all(rownames(colData) %in% colnames(countData)))
countData <- countData[, rownames(colData)]
print(all(rownames(colData) == colnames(countData)))

colData$Condition <- as.factor(colData$Condition)
design_formula <- ~ factor(Condition)
countData[countData == 0] <- 1
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = design_formula)

# cat("Number of genes", nrow(dds), "\n")
# gene_count_thresh <- 10 # Removing genes with low counts
# dds1 <- dds[ rowSums(counts(dds)) >= gene_count_thresh, ]
# cat("Number of genes, after thresholding count with", gene_count_thresh, "is", nrow(dds1), "\n")

dds <- DESeq(dds)
res <- results(dds, tidy = TRUE)
res <- res[order(res$padj),]
write.csv(res, file = "results.csv", row.names = FALSE) 
