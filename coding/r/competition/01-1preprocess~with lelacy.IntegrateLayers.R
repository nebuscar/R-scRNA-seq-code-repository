library(Seurat)
library(ggplot2)

##########1.Load data##########
counts.ctrl1<- Read10X("./data/raw/scRNAseq.wsy/Ctrl1_202105/", gene.column = 1)
counts.ctrl2 <- Read10X("./data/raw/scRNAseq.wsy/Ctrl2_202105/", gene.column = 1)
counts.pos1 <- Read10X("./data/raw/scRNAseq.wsy/Pos1_202104/", gene.column = 1)
counts.pos2 <- Read10X("./data/raw/scRNAseq.wsy/Pos2_202105/", gene.column = 1)

# Create Seurat.object indivdually
sc.ctrl1 <- CreateSeuratObject(counts = counts.ctrl1, project = "Ctrl1", min.cells = 3, min.features = 200)
sc.ctrl2 <- CreateSeuratObject(counts = counts.ctrl2, project = "Ctrl2", min.cells = 3, min.features = 200)
sc.pos1 <- CreateSeuratObject(counts = counts.pos1, project = "Pos1", min.cells = 3, min.features = 200)
sc.pos2 <- CreateSeuratObject(counts = counts.pos2, project = "Pos2", min.cells = 3, min.features = 200)

# merge
sc.combined <- merge(sc.ctrl1, y = c(sc.ctrl2, sc.pos1, sc.pos2), add.cell.ids = c("Ctrl1", "Ctrl2", "Pos1", "Pos2"))

# rename Idents
rename_vector <- c("Pos1" = "STM_1", "Pos2" = "STM_2", "Ctrl1" = "Control_1", "Ctrl2" = "Control_2")
sc.combined <- RenameIdents(object = sc.combined, rename_vector)

rm(counts.ctrl1, counts.ctrl2, counts.pos1, counts.pos2,
   sc.ctrl1, sc.ctrl2, sc.pos1, sc.pos2)

# add metadata:percent.mt
sc.combined[["percent.mt"]] <- PercentageFeatureSet(sc.combined, pattern = "^MT-")
vln_plot <- VlnPlot(sc.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, 
        raster = F, pt.size = 0)
ggsave(filename = "./results/competition/violin_plot~before QC.png", plot = vln_plot, width = 10, height = 6)
table(sc.combined@meta.data$orig.ident)

saveRDS(sc.combined, "tmp/competition/sc.combined~with QC.rds")

##########2.Quality control##########
# Done

##########3-1Non-Integration##########
sc.combined_unintegrated <- NormalizeData(sc.combined)
sc.combined_unintegrated <- FindVariableFeatures(sc.combined_unintegrated)
sc.combined_unintegrated <- ScaleData(sc.combined_unintegrated)
sc.combined_unintegrated <- RunPCA(sc.combined_unintegrated)
sc.combined_unintegrated <- FindNeighbors(sc.combined_unintegrated, reduction = "pca", dims = 1:30)
sc.combined_unintegrated <- FindClusters(sc.combined_unintegrated, resolution = 0.1, cluster.name = "unintegrated_clusters_res=0.1")
sc.combined_unintegrated <- FindClusters(sc.combined_unintegrated, resolution = 0.8, cluster.name = "unintegrated_clusters_res=0.8")
sc.combined_unintegrated <- RunUMAP(sc.combined_unintegrated, reduction = "pca", dims = 1:30, reduction.name = "umap.unintegrated")
pdf("./results/competition/dimplot_unintegrated.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.combined_unintegrated, reduction = "umap.unintegrated", 
        group.by = c("orig.ident", "seurat_clusters"),
        raster = FALSE)
dev.off()

pdf("./results/competition/dimplot_unintegrated_res=0.1.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.combined_unintegrated, reduction = "umap.unintegrated", 
        group.by = c("orig.ident", "unintegrated_clusters_res=0.1"),
        raster = FALSE)
dev.off()

pdf("./results/competition/dimplot_unintegrated_res=0.8.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.combined_unintegrated, reduction = "umap.unintegrated", 
        group.by = c("orig.ident", "unintegrated_clusters_res=0.8"),
        raster = FALSE)
dev.off()

##########3-2CCA-Integration##########
# sc.combined[["RNA"]] <- split(sc.combined[["RNA"]], f = sc.combined$orig.ident)
sc.combined.cca_integration <- NormalizeData(sc.combined)
sc.combined.cca_integration <- FindVariableFeatures(sc.combined.cca_integration)
sc.combined.cca_integration <- ScaleData(sc.combined.cca_integration)
sc.combined.cca_integration <- RunPCA(sc.combined.cca_integration)
sc.combined.cca_integration <- IntegrateLayers(sc.combined.cca_integration,
                             method = CCAIntegration,
                             orig.reduction = "pca",
                             new.reduction = "integrated.cca",
                             verbose = F)
sc.combined.cca_integration[["RNA"]] <- JoinLayers(sc.combined.cca_integration[["RNA"]])
sc.combined.cca_integration <- FindNeighbors(sc.combined.cca_integration, reduction = "integrated.cca", dims = 1:30)
sc.combined.cca_integration <- FindClusters(sc.combined.cca_integration, resolution = 0.1, cluster.name = "integrated_cca_clusters_res.0.1")
sc.combined.cca_integration <- FindClusters(sc.combined.cca_integration, resolution = 0.8,cluster.name = "integrated_cca_clusters_res.0.8")
sc.combined.cca_integration <- FindClusters(sc.combined.cca_integration, resolution = 1,cluster.name = "integrated_cca_clusters_res.1")
sc.combined.cca_integration <- RunUMAP(sc.combined.cca_integration, dims = 1:30, reduction = "integrated.cca", reduction.name = "umap.cca")

# seurat clusters
pdf("./results/competition/dimplot_integrated_cca.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.combined.cca_integration, reduction = "umap.cca", group.by = c("orig.ident", "seurat_clusters"),
        raster = FALSE, label = TRUE, label.size = 8)+ NoLegend()
dev.off()

# res=0.1
pdf("./results/competition/dimplot_integrated_cca_res=0.1.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.combined.cca_integration, reduction = "umap.cca", group.by = c("orig.ident", "integrated_cca_clusters_res.0.1"),
        raster = FALSE, label = TRUE, label.size = 8)+ NoLegend()
dev.off()

# res=0.8
pdf("./results/competition/dimplot_integrated_cca_res=0.8.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.combined.cca_integration, reduction = "umap.cca", group.by = c("orig.ident", "integrated_cca_clusters_res.0.8"),
        raster = FALSE, label = TRUE, label.size = 8)+ NoLegend()
dev.off()

# res=1
pdf("./results/competition/dimplot_integrated_cca_res=1.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.combined.cca_integration, reduction = "umap.cca", group.by = c("orig.ident", "integrated_cca_clusters_res.1"),
        raster = FALSE, label = TRUE, label.size = 8)+ NoLegend()
dev.off()

saveRDS(sc.combined.cca_integration, file = "./tmp/competition/sc.combined.after_umap.cca~without QC.rds")

##########3-3RPCA-Integration##########
# sc.tcells[["RNA"]] <- split(sc.tcells[["RNA"]], f = sc.tcells$orig.ident)
sc.combined.rpca_integration <- NormalizeData(sc.combined)
sc.combined.rpca_integration <- FindVariableFeatures(sc.combined.rpca_integration)
sc.combined.rpca_integration <- ScaleData(sc.combined.rpca_integration)
sc.combined.rpca_integration <- RunPCA(sc.combined.rpca_integration)
sc.combined.rpca_integration <- IntegrateLayers(sc.combined.rpca_integration,
                             method = RPCAIntegration,
                             orig.reduction = "pca",
                             new.reduction = "integrated.rpca",
                             verbose = F,
                             )
sc.combined.rpca_integration[["RNA"]] <- JoinLayers(sc.combined.rpca_integration[["RNA"]])
sc.combined.rpca_integration <- FindNeighbors(sc.combined.rpca_integration, reduction = "integrated.rpca", dims = 1:30)
sc.combined.rpca_integration <- FindClusters(sc.combined.rpca_integration, resolution = 0.1, cluster.name = "integrated_rpca_clusters_res.0.1")
sc.combined.rpca_integration <- FindClusters(sc.combined.rpca_integration, resolution = 0.8,cluster.name = "integrated_rpca_clusters_res.0.8")
sc.combined.rpca_integration <- FindClusters(sc.combined.rpca_integration, resolution = 1,cluster.name = "integrated_rpca_clusters_res.1")
sc.combined.rpca_integration <- RunUMAP(sc.combined.rpca_integration, dims = 1:30, reduction = "integrated.rpca", reduction.name = "umap.rpca")

pdf("./results/competition/dimplot_integrated_rpca.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.combined.rpca_integration, reduction = "umap.rpca", group.by = c("orig.ident", "seurat_clusters"),
        raster = FALSE, label = TRUE, label.size = 8)+ NoLegend()
dev.off()

pdf("./results/competition/dimplot_integrated_rpca_res=0.1.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.combined.rpca_integration, reduction = "umap.rpca", group.by = c("orig.ident", "integrated_rpca_clusters_res.0.1"),
        raster = FALSE, label = TRUE, label.size = 8)+ NoLegend()
dev.off()

pdf("./results/competition/dimplot_integrated_rpca_res=0.8.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.combined.rpca_integration, reduction = "umap.rpca", group.by = c("orig.ident", "integrated_rpca_clusters_res.0.8"),
        raster = FALSE, label = TRUE, label.size = 8)+ NoLegend()
dev.off()

pdf("./results/competition/dimplot_integrated_rpca_res=1.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.combined.rpca_integration, reduction = "umap.rpca", group.by = c("orig.ident", "integrated_rpca_clusters_res.1"),
        raster = FALSE, label = TRUE, label.size = 8)+ NoLegend()
dev.off()

saveRDS(sc.combined.rpca_integration, file = "./tmp/competition/sc.combined.after_umap.rpca~without QC.rds")
