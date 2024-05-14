setwd("E:/scRNA-seq/GBM")
library(Scissor)
library(stringi)
library(rtracklayer)

gtf <- rtracklayer::import(gzfile("data/download/gencode.v45.chr_patch_hapl_scaff.annotation.gtf.gz"))
gtf <- as.data.frame(gtf)
gtf <- gtf[gtf$type == "gene", ]
gtf <- gtf[, c("gene_id", "gene_name")]
gtf$gene_id <-  stri_sub(gtf$gene_id,1,15)
saveRDS(gtf,file = "gtf.rds")
