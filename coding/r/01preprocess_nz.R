library(Seurat)


##########Create seurat object##########
sc.raw <- Read10X("your_data_path")
sc.raw <- readRDS("your_data_path")
sc.data <- CreateSeuratObject(sc.raw, min.cells = 3, min.features = 200)

##########Standard preprocessing workflow##########
sc.data[["percent.mt"]] <- PercentageFeatureSet(sc.data, pattern = "^MT-")
# sc.data[["RNA"]] <- split(sc.data[["RNA"]], f = sc.data$orig.ident) split layers if multi datasets
sc.data <- NormalizeData(sc.data)
sc.data <- FindVariableFeatures(sc.data)
sc.data <- ScaleData(sc.data)
sc.data <- RunPCA(sc.data)

##########CCA.IntegrateLayers##########
# sc.data <- IntegrateLayers(object = sc.data, 
#                           method = RPCAIntegration,
#                           orig.reduction = "pca", 
#                           new.reduction = "integrated.rpca",
#                           verbose = F, 
#                           k.weight = 32)

sc.data <- FindNeighbors(sc.data, reduction = "pca", dims = 1:30) # reduction: choose others if IntegrateLayerss
sc.data <- FindClusters(sc.data, resolution = 0.8, cluster.name = "res = 0.1")
sc.data <- FindClusters(sc.data, resolution = 0.8, cluster.name = "res = 0.8")
sc.data <- RunUMAP(sc.data, reduction = "pca", dims = 1:30, reduction.name = "umap")

