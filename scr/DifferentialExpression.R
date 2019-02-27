rm(list = ls())

library("DESeq2")
library("gplots")
library("RColorBrewer")

sampleFiles <- c("ref1.count","ref2.count","ref3.count","fh1.count","fh2.count","fh3.count")
sampleTable <- data.frame(sampleName= c("ref1","ref2","ref3","fh1","fh2","fh3"),
                          fileName=c("ref1.count","ref2.count","ref3.count","fh1.count","fh2.count","fh3.count"),
                          condition= c("poor","poor","poor","rich","rich","rich"))
summary(sampleTable)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = "~/Analysis_Transcriptomic/results/4_Read_counts/", design = ~condition)
ddsHTSeq
colData(ddsHTSeq)$condition <-factor(colData(ddsHTSeq)$condition, levels=c('poor','rich'))
dds <- DESeq(ddsHTSeq)
normalize.counts <- as.data.frame(counts( dds, normalized=TRUE))
head(normalize.counts)
write.table(normalize.counts, file="~/Analysis_Transcriptomic/results/5_differential_expression/normalized.count", sep="\t", quote = FALSE)
summary(normalize.counts)
res <-results(dds)
res <-res[order(res$padj),]
head(res)
res <- results(dds)
#Expression versus significance
pdf("~/Analysis_Transcriptomic/results/6_visualize/fig1.pdf")
plotMA(dds, ylim=c(-10,10), main="DESeq2, Expression versus significance")
dev.off()

#Expression heatmap for the most expressed genes

rld <- rlogTransformation(dds, blind = TRUE)
vsd <-varianceStabilizingTransformation(dds, blind = TRUE)
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing = TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

pdf("~/Analysis_Transcriptomic/results/6_visualize/fig2.pdf")
heatmap.2(counts(dds, normalized= TRUE)[select,], col= hmcol, Rowv = FALSE,
          Colv = FALSE, scale = "none", dendrogram = "none",
          trace = "none", margin = c(10,6), main="                   Expression heatmap for
          the most expressed genes")
dev.off()

#Heatmap of similarity between replicates

distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition, sampleFiles, sep=" : "))
hc <- hclust(distsRL)
pdf("~/Analysis_Transcriptomic/results/6_visualize/fig3.pdf")
heatmap.2(mat, Rowv = as.dendrogram(hc), symm=TRUE, trace="none",
          col = rev(hmcol), margin = c(13,13), main = "                Heatmap of similarity
          between replicates")
dev.off()

#PCA plot
png("~/Analysis_Transcriptomic/results/6_visualize/fig4.png")
print(plotPCA(rld, intgroup=c("condition")))
dev.off()

#Writting all calculations to a file
write.table(res, file = "~/Analysis_Transcriptomic/results/5_differential_expression/diffExp.tab", sep = "\t", quote = FALSE)

#Extracting genes with p-adjust value below 0.05
resSig <- subset(res, padj < 0.05)
write.table(resSig, file = "~/Analysis_Transcriptomic/results/5_differential_expression/diffExp.0.05.tab", sep="\t", quote=FALSE)

resSigs <- subset(res, padj < 0.01)
write.table(resSigs, file = "~/Analysis_Transcriptomic/results/5_differential_expression/diffExp.0.01.tab", sep="\t", quote=FALSE)
