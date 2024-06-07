workdir = "F:/Projects/Skull_scRNAseq"  # where you place the data and results
setwd(workdir)
library(dplyr)
library(Seurat)

file.list <- list.files("rawData/", pattern = "matrix$")

sk.rawData <- list()
for (file in file.list){
  sk.rawData[[file]] <- Read10X(paste("rawData/", file, sep = ""))
  colnames(sk.rawData[[file]]) <- paste(strsplit(file, "-")[[1]][1], colnames(sk.rawData[[file]]), sep = "_")
}

sk.scData <- CreateSeuratObject(sk.rawData, project = "skull", min.cells = 3, min.features = 200)
#sk.scData@meta.data$orig.ident <- vapply(strsplit(colnames(sk.scData), "_"), `[`, 1, FUN.VALUE=character(1))

sk.scData[["percent.mt"]] <- PercentageFeatureSet(sk.scData, pattern = "^MT-")
VlnPlot(sk.scData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, 
        raster = F, pt.size = 0)
table(sk.scData@meta.data$orig.ident)

# QC
sk.scData <- subset(sk.scData, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
VlnPlot(sk.scData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, 
        raster = F, pt.size = 0)
table(sk.scData@meta.data$orig.ident)

# non-integration
# run standard analysis workflow
sk.scData <- NormalizeData(sk.scData)
sk.scData <- FindVariableFeatures(sk.scData)
sk.scData <- ScaleData(sk.scData)
sk.scData <- RunPCA(sk.scData)
sk.scData <- FindNeighbors(sk.scData, dims = 1:30, reduction = "pca")
sk.scData <- FindClusters(sk.scData, resolution = 0.1, cluster.name = "unintegrated_clusters")
sk.scData <- RunUMAP(sk.scData, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
pdf("./tmp/dimplot_unintegrated.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sk.scData, reduction = "umap.unintegrated", 
        group.by = c("orig.ident", "seurat_clusters"),
        raster = FALSE)
dev.off()

################################################################################
# Integration
sk.scData <- IntegrateLayers(object = sk.scData, method = CCAIntegration, 
                             orig.reduction = "pca", new.reduction = "integrated.cca",
                             verbose = FALSE)
# re-join layers after integration
sk.scData[["RNA"]] <- JoinLayers(sk.scData[["RNA"]])
sk.scData <- FindNeighbors(sk.scData, reduction = "integrated.cca", dims = 1:30)
sk.scData <- FindClusters(sk.scData, resolution = 0.1)
sk.scData <- RunUMAP(sk.scData, dims = 1:30, reduction = "integrated.cca")

pdf("./tmp/dimplot_integrated.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sk.scData, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"),
        raster = FALSE, label = TRUE, label.size = 8)+ NoLegend()
dev.off()

saveRDS(sk.scData, file = "./tmp/sk.scData.rds")
################################################################################

sk.scData <- readRDS("./tmp/sk.scData.rds")
sk.scData <- FindClusters(sk.scData, resolution = 0.5)
sk.scData <- RunUMAP(sk.scData, dims = 1:30, reduction = "integrated.cca")

pdf("./result/dimplot_integrated_05.pdf", width = 6, height = 6, useDingbats = F)
DimPlot(sk.scData, reduction = "umap", group.by = c("seurat_clusters"),
        raster = FALSE, label = TRUE, label.size = 4)+ NoLegend()
dev.off()
