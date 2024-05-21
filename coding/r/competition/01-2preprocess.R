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

##########Remove BatchEffect##########
anchors <- FindIntegrationAnchors(object.list = seurat_objects, dims = 1:30)
sc.combined <- IntegrateData(anchorset = anchors, dims = 1:30)

saveRDS(sc.combined, "tmp/competition/sc.combined~after IntegrateAnchors.rds")

# Set default Assay:integrated
DefaultAssay(sc.combined) <- "integrated"


# rename Idents
# rename_vector <- c("Pos1" = "STM_1", "Pos2" = "STM_2", "Ctrl1" = "Control_1", "Ctrl2" = "Control_2")
# sc.combined <- RenameIdents(object = sc.combined, rename_vector)

##########
# ScaleData and RunPCA
sc.combined <- ScaleData(sc.combined, verbose = FALSE)
sc.combined <- RunPCA(sc.combined, npcs = 30, verbose = FALSE)

# UMAP
sc.combined <- RunUMAP(sc.combined, reduction = "pca", dims = 1:30)

# add metadata:percent.mt
sc.combined[["percent.mt"]] <- PercentageFeatureSet(sc.combined, pattern = "^MT-")
vln_plot <- VlnPlot(sc.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, 
                    raster = F, pt.size = 0)
ggsave(filename = "./results/competition/violin_plot~after IntegrateAnchors .png", plot = vln_plot, width = 10, height = 6)
table(sc.combined@meta.data$orig.ident)

# seurat clusters
pdf("./results/competition/dimplot_integrated_after IntegrateAnchors.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.combined, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"),
        raster = FALSE, label = TRUE, label.size = 8)+ NoLegend()
dev.off()

saveRDS(sc.combined, file = "./tmp/competition/sc.combined.after_umap.cca~without QC.rds")
