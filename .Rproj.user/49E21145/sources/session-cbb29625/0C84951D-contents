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
# sc.combined <- RenameIdents(object = sc.combined, `Pos1` = "STM_1")

rm(counts.ctrl1, counts.ctrl2, counts.pos1, counts.pos2,
   sc.ctrl1, sc.ctrl2, sc.pos1, sc.pos2)

# add metadata:percent.mt
sc.combined[["percent.mt"]] <- PercentageFeatureSet(sc.combined, pattern = "^MT-")
vln_plot <- VlnPlot(sc.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, 
        raster = F, pt.size = 0)
ggsave(filename = "./results/competition/violin_plot~before QC.png", plot = vln_plot, width = 10, height = 6)
table(sc.combined@meta.data$orig.ident)

saveRDS(sc.combined, "tmp/competition/sc.combined.rds")

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

##########3-2CCA-Intergration##########
sc.combined[["RNA"]] <- split(sc.combined[["RNA"]], f = sc.combined$orig.ident)
sc.combined <- NormalizeData(sc.combined)
sc.combined <- FindVariableFeatures(sc.combined)
sc.combined <- ScaleData(sc.combined)
sc.combined <- RunPCA(sc.combined)
sc.combined <- IntegrateLayers(sc.combined,
                             method = CCAIntegration,
                             orig.reduction = "pca",
                             new.reduction = "integrated.cca",
                             verbose = F,
                             k.weight = 32,
)
sc.tcells[["RNA"]] <- JoinLayers(sc.tcells[["RNA"]])
sc.tcells <- FindNeighbors(sc.tcells, reduction = "integrated.cca", dims = 1:30)
sc.tcells <- FindClusters(sc.tcells, resolution = 0.1, cluster.name = "integrated_cca_clusters_res.0.1")
sc.tcells <- FindClusters(sc.tcells, resolution = 0.8,cluster.name = "integrated_cca_clusters_res.0.8")
sc.tcells <- RunUMAP(sc.tcells, dims = 1:30, reduction = "integrated.cca", reduction.name = "umap.cca")

pdf("./results/dimplot_integrated_cca.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.tcells, reduction = "umap.cca", group.by = c("orig.ident", "seurat_clusters"),
        raster = FALSE, label = TRUE, label.size = 8)+ NoLegend()
dev.off()

pdf("./results/dimplot_integrated_cca_res=0.1.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.tcells, reduction = "umap.cca", group.by = c("orig.ident", "integrated_cca_clusters_res.0.1"),
        raster = FALSE, label = TRUE, label.size = 8)+ NoLegend()
dev.off()

pdf("./results/dimplot_integrated_caa_res=0.8.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.tcells, reduction = "umap.cca", group.by = c("orig.ident", "integrated_cca_clusters_res.0.8"),
        raster = FALSE, label = TRUE, label.size = 8)+ NoLegend()
dev.off()

saveRDS(sc.tcells, file = "./tmp/sc.tcells.after_umap.cca~without QC.rds")

##########3-3RPCA-Intergration##########
sc.tcells[["RNA"]] <- split(sc.tcells[["RNA"]], f = sc.tcells$orig.ident)
sc.tcells <- NormalizeData(sc.tcells)
sc.tcells <- FindVariableFeatures(sc.tcells)
sc.tcells <- ScaleData(sc.tcells)
sc.tcells <- RunPCA(sc.tcells)
sc.tcells <- IntegrateLayers(sc.tcells,
                             method = RPCAIntegration,
                             orig.reduction = "pca",
                             new.reduction = "integrated.rpca",
                             verbose = F,
                             k.weight = 32,
)
sc.tcells[["RNA"]] <- JoinLayers(sc.tcells[["RNA"]])
sc.tcells <- FindNeighbors(sc.tcells, reduction = "integrated.rpca", dims = 1:30)
sc.tcells <- FindClusters(sc.tcells, resolution = 0.1, cluster.name = "integrated_rpca_clusters_res.0.1")
sc.tcells <- FindClusters(sc.tcells, resolution = 0.8,cluster.name = "integrated_rpca_clusters_res.0.8")
sc.tcells <- RunUMAP(sc.tcells, dims = 1:30, reduction = "integrated.rpca", reduction.name = "umap.rpca")

pdf("./results/dimplot_integrated_rpca.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.tcells, reduction = "umap.rpca", group.by = c("orig.ident", "seurat_clusters"),
        raster = FALSE, label = TRUE, label.size = 8)+ NoLegend()
dev.off()

pdf("./results/dimplot_integrated_rpca_res=0.1.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.tcells, reduction = "umap.rpca", group.by = c("orig.ident", "integrated_rpca_clusters_res.0.1"),
        raster = FALSE, label = TRUE, label.size = 8)+ NoLegend()
dev.off()

pdf("./results/dimplot_integrated_rpca_res=0.8.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.tcells, reduction = "umap.rpca", group.by = c("orig.ident", "integrated_rpca_clusters_res.0.8"),
        raster = FALSE, label = TRUE, label.size = 8)+ NoLegend()
dev.off()

saveRDS(sc.tcells, file = "./tmp/sc.tcells.after_umap.rpca~without QC.rds")


