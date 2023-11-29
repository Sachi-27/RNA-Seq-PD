pdf("plots.pdf")

#### PLOTCOUNTS
top_genes <- head(res[, 1], 6)

par(mfrow = c(2, 3))
for (gene_id in top_genes) {
  plotCounts(dds, gene = gene_id, intgroup = "Condition")
}

#### VOLCANO PLOT
par(mfrow=c(1,1))

# Basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))

padjthr <- 0.01

# Colored points for upregulated genes
with(subset(res, padj < padjthr & log2FoldChange > 0), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, padj < padjthr & log2FoldChange > 1), points(log2FoldChange, -log10(pvalue), pch=20, col="red2"))
with(subset(res, padj < padjthr & log2FoldChange > 2), points(log2FoldChange, -log10(pvalue), pch=20, col="red4"))

# Colored points for downregulated genes
with(subset(res, padj < padjthr & log2FoldChange < 0), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj < padjthr & log2FoldChange < -1), points(log2FoldChange, -log10(pvalue), pch=20, col="blue2"))
with(subset(res, padj < padjthr & log2FoldChange < -2), points(log2FoldChange, -log10(pvalue), pch=20, col="blue4"))

# Counts
up2 <- sum(res$padj < padjthr & (res$log2FoldChange > 2))
up1 <- sum(res$padj < padjthr & (res$log2FoldChange > 1)) - up2
up <- sum(res$padj < padjthr & (res$log2FoldChange > 0)) - up1 - up2
up <- ifelse(is.na(up), 0, up)
down2 <- sum(res$padj < padjthr & (res$log2FoldChange < -2))
down1 <- sum(res$padj < padjthr & (res$log2FoldChange < -1)) - down2 
down <- sum(res$padj < padjthr & (res$log2FoldChange < 0)) - down1 - down2
down <- ifelse(is.na(down), 0, down)
cat(up, up1, up2, down, down1, down2, "\n")

# Add legend
legend("top", legend=c(up, up1, up2, down, down1, down2), 
       col=c("red", "red2", "red4", "blue", "blue2", "blue4"), pch=20, bty="n")

abline(v = c(-1, 1, -2, 2), col="gray", lty=2)

# Add gene names for points with y-axis > 30
#selected <- -log10(res$padj) > 30
#text(res$log2FoldChange[selected], -log10(res$padj)[selected], lab=rownames(res)[selected], cex=0.5)


library("affy")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("ggfortify")

# Remove Genes with low counts
dds1 <- dds[ rowSums(counts(dds)) >= 10, ]
nrow(dds1)
PCAdata <- prcomp(t(assay(dds1)))
autoplot(PCAdata, data = colData,colour = "Condition", label = FALSE, main="PCA Unnormalised")

clusters2 <- hclust(dist( t( assay(dds1) ) ), method ="ward.D")
plot(clusters2, main = "Dendrogram Unnormalised",labels = FALSE)

plotDensity(assay(dds1), col=1:24,lwd=2,lty=1,xlab("Density"),ylab("Counts"), main = "Density plot Unnormalised")

vst = vst(dds1, blind=FALSE)
PCAdata <- prcomp(t(assay(vst)))
autoplot(PCAdata, data = colData,colour = "Condition",label = FALSE, main="PCA Normalised")

clusters2 <- hclust(dist( t( assay(vst) ) ),method ="ward.D")
plot(clusters2, main = "Dendrogram Normalised", label = colData$Condition)

plot(clusters2, main = "Dendrogram Normalised", labels = FALSE)

plotDensity(assay(vst), lwd=2,lty=1,xlab("Density"),ylab("Counts"), main = "Density plot Normalised")

sampleDists <- dist(t(assay(vst)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$SampleID, vst$Condition, sep="-")
colnames(sampleDistMatrix) <- paste(dds1$SampleID, dds1$Condition, sep="-")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Specify the order of rows (samples) in the heatmap
order_rows <- order(vst$Condition, decreasing = FALSE)

pheatmap(sampleDistMatrix[order_rows, order_rows],
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors, main = "Sample Distance Matrix ")
		 
# MA PLOTS	 
plotMA(dds)

# CUSTOM MA PLOT
library(ggplot2)
library(scales) # needed for oob parameter
library(viridis)

# Coerce to a data frame
deseq2ResDF <- as.data.frame(res)

# Examine this data frame
head(deseq2ResDF)

# Set a boolean column for significance
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .1, "Significant", NA)

## Dispersion plot
plotDispEsts( dds, ylim = c(1e-6, 1e1) )

# Histogram of p-value
hist( res$pvalue, breaks=20, col="grey" )

# Close the PDF device
dev.off()

##############
pdf("plots2.pdf")
# Plot the results similar to DEseq2
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", size=2) + labs(x="mean of normalized counts", y="log fold change") + scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()

# Let's add some more detail
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=padj)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") + labs(x="mean of normalized counts", y="log fold change") + scale_colour_viridis(direction=-1, trans='sqrt') + theme_bw() + geom_density_2d(colour="black", size=2)

dev.off()
