library(Seurat)
library(ggplot2)

##########1.Load data##########
data_location <- "./data/raw/scRNAseq.wsy/"
samples <- c("Ctrl1_202105", "Ctrl2_202105", "Pos1_202104", "Pos2_202105")

# Create Seurat.object indivdually
seurat_objects <- lapply(samples, function(sample) {
  counts <- Read10X(file.path(data_location, sample), gene.column = 1)
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample, min.cells = 3, min.features = 200)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  return(seurat_obj)
})

##########2.Remove BatchEffect:IntegrateData##########
features <- SelectIntegrationFeatures(seurat_objects)

for (i in seq_along(seurat_objects)) {
  seurat_objects[[i]] <- ScaleData(seurat_objects[[i]], features = features)
  seurat_objects[[i]] <- RunPCA(seurat_objects[[i]], features = features)
}

anchors <- FindIntegrationAnchors(object.list = seurat_objects, dims = 1:30)
sc.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

saveRDS(sc.integrated, "tmp/competition/sc.integrated~after IntegrateAnchors.rds")

# Set default Assay:integrated
DefaultAssay(sc.integrated) <- "integrated"


# rename Idents
# rename_vector <- c("Pos1_202104" = "STM_1", "Pos2_202105" = "STM_2", "Ctrl1_202105" = "Control_1", "Ctrl2_202105" = "Control_2")
# sc.integrated <- RenameIdents(object = sc.integrated, rename_vector)


##########3.Standard work flow##########
# ScaleData and RunPCA
sc.integrated <- ScaleData(sc.integrated, verbose = FALSE)
sc.integrated <- RunPCA(sc.integrated, npcs = 30, verbose = FALSE)

# UMAP
sc.integrated <- RunUMAP(sc.integrated, reduction = "pca", dims = 1:30, reduction.name = "umap.integrated")
sc.integrated <- FindNeighbors(sc.integrated, dims = 1:30)
sc.integrated <- FindClusters(sc.integrated, resolution = 1, cluster.name = "integrated_anchors_res.1")

# seurat clusters
pdf("./results/competition/dimplot_integrated_res=1~after IntegrateAnchors.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.integrated, reduction = "umap.integrated", group.by = c("orig.ident", "seurat_clusters"),
        raster = FALSE, label = TRUE, label.size = 8)+ NoLegend()
dev.off()

saveRDS(sc.integrated, file = "./tmp/competition/sc.integrated.after_umap.cca.rds")
