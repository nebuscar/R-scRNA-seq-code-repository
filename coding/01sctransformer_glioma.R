library(Seurat)

glioma_data <- Read10X("data/download/sc_dataset/raw")
glioma <- CreateSeuratObject(counts = glioma_data, project = "glioma", min.cells = 3, min.features = 200)
#rm(glioma_data)
metadata <- read.csv(file("data/download/sc_dataset/Meta_GBM.txt.gz"), header = T)
tcells.id <- metadata[metadata$Assignment %in% "TCells",][, 1]
glioma@meta.data$barcode <- colnames(glioma)
tcells <- subset(glioma, barcode %in% tcells.id)
#tcells <- subset(glioma, subset = orig.ident == "GSM5518615")
rm(glioma)
saveRDS(tcells,file = "data/seurat/tcells_v5.rds")
saveRDS(glioma, file = "data/seurat/glioma_v5.rds")

# tcells <- PercentageFeatureSet(tcells, pattern = "^MT-", col.name = "percent.mt")
# tcells <- SCTransform(tcells, vars.to.regress = "percent.mt", verbose = F)
# 
# tcells <- RunPCA(tcells, verbose = FALSE)
# tcells <- RunUMAP(tcells, dims = 1:30, verbose = FALSE)
# 
# tcells <- FindNeighbors(tcells, dims = 1:30, verbose = FALSE)
# tcells <- FindClusters(tcells, verbose = FALSE)
# DimPlot(tcells, label = TRUE)
# 
