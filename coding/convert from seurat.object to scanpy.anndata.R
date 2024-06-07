library(Seurat)

sc.tcells <- readRDS("./tmp/sc.tcells.after_umap.rpca~without QC.rds")


##########!SCP##########
library(SCP)
data("pancreas1k")
adata <- srt_to_adata(pancreas1k)
adata$write_h5ad("pancreas1k.h5ad")

##########!zellkonverter##########
sce_obj <- as.SingleCellExperiment(sce_obj, assay = c("RNA"))
library(zellkonverter)
writeH5AD(sce_obj, "sce_obj.h5ad", X_name = 'counts')

