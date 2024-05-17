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


##########!reticulate##########
require(Seurat)
require(reticulate)

seu <- readRDS('your_path_seurat_object_rds')

# load python anndata package
anndata <- reticulate::import('anndata')
# create anndata object
adata <- anndata$AnnData(X = seu@assays$RNA@layers$counts, obs = data.frame(row.names = rownames(seu)), var = seu@meta.data )
adata$write("your_path_scanpy_obj_h5ad")

# Of note that adata require an inversion in python scanpy
import scanpy as sc
adata = sc.read_h5ad("your_path_scanpy_obj_h5ad")
adata = adata.T
##########SeuratDisk##########
library(SeuratDisk)
SaveH5Seurat(sc.tcells, filename = "./tmp/pbmc3k.h5Seurat")
Convert("./tmp/pbmc3k.h5Seurat", dest = "h5ad")
