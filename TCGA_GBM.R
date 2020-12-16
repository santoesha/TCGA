#setwd("~/project/IntellanceII/TCGA")

#Load libraries
library(tidyverse)
library(DESeq2)

#Load TCGA data
counts <- read.delim("tcga_GBM.txt",row.names = c(1),header=TRUE,sep=" ")
counts <- counts[rowSums(counts) >= 3 * ncol(counts),]

# ---- Normalization and VST ----
condition <- as.factor(paste0('r',sample(c(1,2),ncol(counts),T)))
dds <- DESeqDataSetFromMatrix(countData=as.matrix(counts), colData=data.frame(row.names=colnames(counts), condition), design=~condition)
vst <- varianceStabilizingTransformation(DESeq(dds),blind=TRUE)
plotPCA(vst) #no outliers observed
counts <- assay(vst)

#write.csv(counts, "vst_normalized_tcga_gbm_counts.csv")

plot(counts["EGFR",],counts["FGF3",], xlab="counts EGFR", ylab="counts FGFR3",main="TCGA GBM counts")
