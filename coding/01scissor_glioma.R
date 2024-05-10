setwd("E:/scRNA-seq/GBM")

library(Scissor)
library(stringi)
library(rtracklayer)

gtf <- rtracklayer::import(gzfile("data/download/gencode.v45.chr_patch_hapl_scaff.annotation.gtf.gz"))
gtf <- as.data.frame(gtf)
gtf <- gtf[gtf$type == "gene", ]
gtf <- gtf[, c("gene_id", "gene_name")]
gtf$gene_id <-  stri_sub(gtf$gene_id,1,15)
saveRDS(gtf,file = "data/gtf.rds")

# load data
if(T){


load("data/seurat v4/rdata/tcells_sct_clean.rdata")
#tcells_sub <- subset(tcells, subset = orig.ident == "GSM5518615")
tcells_sub <- tcells
GBM_sc_dataset <- GetAssayData(object = tcells_sub, assay = "RNA", slot = "counts") # get count matrix data


# survival n=602
GBM_bulk_survival <- read.table(file = "data/download/TCGA_GBM_clinicial.txt",row.names=1,header = T,sep = "\t")


# FPKM to TPM
GBM_bulk_dataset_log <- read.table(file = "data/download/TCGA_GBM_buk_fpkm.gz",row.names=1,header = T)
GBM_bulk_dataset <- as.matrix(2^GBM_bulk_dataset_log - 1)

# RNA-seq n=172
if(T){
num_genes <- 60483
num_samples <- 173

fpkm_matrix <- matrix(rexp(num_genes * num_samples, rate = 0.1), nrow = num_genes)
colnames(fpkm_matrix) <- paste0("Sample_", 1:num_samples)
rownames(fpkm_matrix) <- paste0("Gene_", 1:num_genes)

sum_fpkm_per_sample <- colSums(fpkm_matrix)
scaling_factors <- sum_fpkm_per_sample / 1e6
tpm_matrix <- t(t(fpkm_matrix) / scaling_factors)
} # step by step
colnames(tpm_matrix) <- colnames(GBM_bulk_dataset)
rownames(tpm_matrix) <- rownames(GBM_bulk_dataset)
GBM_bulk_dataset <- as.matrix(tpm_matrix)

}

# colnames convert "." to "-"
colnames_dataset <- colnames(GBM_bulk_dataset)
colnames_new_dataset <- gsub("\\.", "-", colnames_dataset)
colnames(GBM_bulk_dataset) <- colnames_new_dataset

#  convert gene ID
tmp <- GBM_bulk_dataset
rownames(tmp) <- stri_sub(rownames(tmp),1,15)
# calculate median for each gene_id and rank genes by median
ids <- data.frame(gene_id=rownames(tmp),
                  median=apply(tmp,1,median)) #计算基因表达中位数，用于之后排序

ids <- merge(ids, gtf, by = "gene_id", all = F)
ids <- ids[order(ids$gene_name, ids$median, decreasing = T),]
table(duplicated(ids$gene_name))
ids <- ids[!duplicated(ids$gene_name),]

# 转化geneid为symbol 
tmp <- tmp[ids$gene_id,] #取出表达矩阵中ids有的行  
rownames(tmp) <- ids[match(rownames(tmp),ids$gene_id),"gene_name"] 
GBM_bulk_dataset <- tmp
colnames(GBM_bulk_dataset) <- stri_sub(colnames(GBM_bulk_dataset),1,15)


# find common 
colnames_survival <- rownames(GBM_bulk_survival)
colnames_new_dataset <- colnames(GBM_bulk_dataset)
common <- intersect(colnames_new_dataset,colnames_survival)
GBM_bulk_dataset <- GBM_bulk_dataset[,common]
GBM_bulk_survival <- GBM_bulk_survival[common,]


GBM_sc_dataset <- Seurat_preprocessing(GBM_sc_dataset,verbose = F)


# phenotype
all(colnames(GBM_bulk_dataset) == GBM_bulk_survival$TCGA_patient_barcode)
GBM_phenotype <- GBM_bulk_survival[,2:3]
colnames(GBM_phenotype) <- c("status","time")


# Scissor
#######################
infos1 <- Scissor(GBM_bulk_dataset,GBM_sc_dataset,GBM_phenotype,alpha = 0.05,
                    family = "cox",Save_file = "Scissor_GBM_survival.RData")
# visualize

# define Scissor+,Scissor-
Scissor_select <- rep("Background cells",ncol(GBM_sc_dataset))
names(Scissor_select) <- colnames(GBM_sc_dataset)

Scissor_select[infos1$Scissor_pos] <- "Scissor+ cell"
Scissor_select[infos1$Scissor_neg] <- "Scissor- cell"

# metatda,add Scissor
GBM_sc_dataset <- AddMetaData(GBM_sc_dataset,metadata = Scissor_select,col.name = "scissor")
Dim_plot1 <- DimPlot(GBM_sc_dataset,reduction = 'umap',group.by = 'scissor',cols = c('grey','royalblue','indianred1'),pt.size = 2)
Dim_plot1
