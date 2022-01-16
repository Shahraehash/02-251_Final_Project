if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library("DESeq2")
require(devtools)

#Cancer comparison file
data_Cancer <- read.csv("CancerComparison.csv")
data_Cancer <- data_Cancer[1:60483,1:7]

#File Comparison File
data_File <- read.csv("DataFileComparison.csv")
data_File <- data_File[1:60483,1:7]

#Sequencing Comparison File
data_Sequencing <- read.csv("SequencingComparison.csv")
data_Sequencing <- data_Sequencing[1:60483, 1:7]

#metadata for the different analysis
metadata <- read.csv("CancerMetadata.csv")
metadata_Cancer <- metadata[1:6,1:4]
metadata_File <- metadata[7:12,1:4]
metadata_Sequencing <- metadata[13:18,1:4]

#Differential Gene Expression Analysis Cancer
dds_Cancer <- DESeqDataSetFromMatrix(countData = data_Cancer, colData = metadata_Cancer, design =~ CancerType, tidy = TRUE)
dds_Cancer <- DESeq(dds_Cancer)
res_Cancer <- results(dds_Cancer)
res_Cancer <- res_Cancer[order(res_Cancer$padj),]
summary(res_Cancer)
head(res_Cancer)
#Plotting of the DE Cancer
par(mfrow=c(1,1))
with(res_Cancer, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE Genes between Cancer Types", xlim=c(-15,15)))
with(subset(res_Cancer, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res_Cancer, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#Differential Gene Expression Analysis File
dds_File <- DESeqDataSetFromMatrix(countData = data_File, colData = metadata_File, design =~ FileType, tidy = TRUE)
dds_File <- DESeq(dds_File)
res_File <- results(dds_File)
res_File <- res_File[order(res_File$padj),]
summary(res_File)
head(res_File)
#Plotting of the DE File
par(mfrow=c(1,1))
with(res_File, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE Genes between File Types", xlim=c(-15,15)))
with(subset(res_File, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res_File, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#Differential Gene Expression Analysis Sequencing
dds_Sequencing <- DESeqDataSetFromMatrix(countData = data_Sequencing, colData = metadata_Sequencing, design =~ SequencingType, tidy = TRUE)
dds_Sequencing <- DESeq(dds_Sequencing)
res_Sequencing <- results(dds_Sequencing)
res_Sequencing <- res_Sequencing[order(res_Sequencing$padj),]
summary(res_Sequencing)
head(res_Sequencing)
#Plotting of the DE File
par(mfrow=c(1,1))
with(res_Sequencing, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE Genes between Sequencing Types", xlim=c(-15,15)))
with(subset(res_Sequencing, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res_Sequencing, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))



