workdir = "your_work_directory"  # where you place the data and results
setwd(workdir)
library(Seurat)


##########Create seurat object##########
sc.raw <- Read10X("your_data_path")
sc.raw <- readRDS("your_data_path")
sc.data <- CreateSeuratObject(sc.raw, min.cells = 3, min.features = 200)

##########Quality control##########
sc.data[["percent.mt"]] <- PercentageFeatureSet(sc.data, pattern = "^MT-")
sc.data <- subset(sc.data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

##########Non-Intergration##########
# Standard preprocessing workflow
sc.data <- NormalizeData(sc.data)
sc.data <- FindVariableFeatures(sc.data)
sc.data <- ScaleData(sc.data)
sc.data <- RunPCA(sc.data)
sc.data <- FindNeighbors(sc.data, reduction = "pca", dims = 1:30) # reduction: choose others if IntegrateLayerss
sc.data <- FindClusters(sc.data, resolution = 0.8, cluster.name = "res = 0.1")
sc.data <- FindClusters(sc.data, resolution = 0.8, cluster.name = "res = 0.8")
sc.data <- RunUMAP(sc.data, reduction = "pca", dims = 1:30, reduction.name = "umap")
pdf("./tmp/dimplot_unintegrated.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.data, reduction = "umap.unintegrated", 
        group.by = c("orig.ident", "seurat_clusters"),
        raster = FALSE)
dev.off()

##########CCA-Intergration##########
sc.data[["RNA"]] <- split(sc.data[["RNA"]], f = sc.data$orig.ident) # split layers if multi datasets
sc.data <- IntegrateLayers(sc.data,
                           method = CCAIntegration,
                           orig.reduction = "pca",
                           new.reduction = "integrated.rpca",
                           verbose = F,
                           k.weight = 32)

sc.data[["RNA"]] <- JoinLayers(sc.data[["RNA"]])
sc.data <- FindNeighbors(sc.data, reduction = "integrated.cca", dims = 1:30)
sc.data <- FindClusters(sc.data, resolution = 0.1)
sc.data <- RunUMAP(sc.data, dims = 1:30, reduction = "integrated.cca")

pdf("./tmp/dimplot_integrated.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.data, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"),
        raster = FALSE, label = TRUE, label.size = 8)+ NoLegend()
dev.off()

saveRDS(sc.data, file = "./tmp/sc.data.rds")