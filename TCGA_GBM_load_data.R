#setwd("~/project/IntellanceII/TCGA")

# ---- Load packages ----
library(dplyr)
library(readr)
library(tidyverse)
library(biomaRt)

# ---- Load TCGA GBM data ---- 
tcga_data <- list.files(path="./TCGA",pattern = "*.counts.gz",full.names = TRUE)
tcga_data <- sapply(tcga_data, read.delim,header=FALSE,simplify=FALSE,USE.NAMES = TRUE)
tcga_data <- bind_rows(tcga_data, .id="file")
colnames(tcga_data) <- c("file","ensembl gene id","number of reads")
tcga_data <- spread(tcga_data,"file", "number of reads")
tcga_data <- remove_rownames(tcga_data)
tcga_data <- column_to_rownames(tcga_data, var = "ensembl gene id")
tcga_data <- tcga_data[-c(1:5),]
rownames(tcga_data) <- gsub("\\..*","",rownames(tcga_data))
colnames(tcga_data) <- gsub("./TCGA/", "",colnames(tcga_data))

# ---- Metadata ----
characteristics <- read_tsv(file='./TCGA/gdc_sample_sheet.2020-07-07.tsv')
exposure <- read_tsv(file='./TCGA/exposure.tsv')
clinical <- read_tsv(file='./TCGA/clinical.tsv')

#Omit Solid Tissue Normal 
solid_tissue_normal <- characteristics[characteristics$`Sample Type` == "Solid Tissue Normal",c(2)] %>% pull('File Name')
tcga_data <- tcga_data[,!(names(tcga_data) %in% solid_tissue_normal)]

#Omit technical replicate TCGA-0156 and TCGA-0211
tcga_data <- tcga_data[,!(names(tcga_data) == "18fbf794-a85b-4b0d-9a5c-c1c0e0ede3d8.htseq.counts.gz")] #0156
tcga_data <- tcga_data[,!(names(tcga_data) == "f67cf68d-56b3-451d-ba21-566664213691.htseq.counts.gz")] #0211

#Match column ids 
common_samples <- intersect(names(tcga_data),characteristics$`File Name`)
names(tcga_data)[names(tcga_data) %in% common_samples] = characteristics[characteristics$`File Name` %in% common_samples,7] %>% pull('Sample ID')

# ---- Convert ensemblID to geneID ----
ensembl_ids = as.character(rownames(tcga_data))
mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

dat = getBM(
  values = ensembl_ids,
  filters = c("ensembl_gene_id"),
  attributes = c("ensembl_gene_id", "external_gene_name", "description"),
  mart = mart
)
dat$description <- NULL

tcga_data$ensembl_gene_id <- rownames(tcga_data)
tcga_data <- merge(dat,tcga_data,by="ensembl_gene_id")
duplicates <- tcga_data$external_gene_name[which(duplicated(tcga_data$external_gene_name))]
tcga_data <- tcga_data[-which(duplicated(tcga_data$external_gene_name)),]
rownames(tcga_data) <- tcga_data$external_gene_name

tcga_data$ensembl_gene_id <- NULL
tcga_data$external_gene_name <- NULL 

#write.table(tcga_data,file="tcga_GBM.txt")
