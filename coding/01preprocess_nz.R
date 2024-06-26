workdir = "your_work_directory"  # where you place the data and results
setwd(workdir)
library(Seurat)

##########1.Create seurat object##########
sc.raw.counts <- Read10X("data/raw/sc.data.SingleCellPortal/")
sc.raw <- CreateSeuratObject(counts = sc.raw.counts, project = "glioma", min.cells = 3, min.features = 200)
saveRDS(sc.raw, "./tmp/sc.raw.glioma.seuratobject.rds")

# subset T cells
raw.meta <- read.csv(gzfile("data/raw/sc.data.SingleCellPortal/Meta_GBM.txt.gz"))[-1, ]
tcells.id <- raw.meta[raw.meta$Assignment %in% "TCells", ][, 1]
sc.raw@meta.data$barcode <- colnames(sc.raw)
# sc.tcells <- subset(sc.raw, barcode %in% tcells.id) # 尽量避免使用subset进行连续子集操作
sc.tcells <- sc.raw[, colnames(sc.raw) %in% tcells.id]

# add metadata:percent.mt
sc.tcells[["percent.mt"]] <- PercentageFeatureSet(sc.tcells, pattern = "^MT-")
VlnPlot(sc.tcells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, 
        raster = F, pt.size = 0)
table(sc.tcells@meta.data$orig.ident)

# subset datasets:cells > 30
orig.ident.meta <- table(sc.tcells@meta.data$orig.ident)>30
sc.tcells <- subset(sc.tcells, subset = orig.ident %in% names(orig.ident.meta[orig.ident.meta == T]))
# count_dataset <- table(sc.tcells@meta.data$orig.ident)
# count_few <- count_dataset[count_dataset < 30]
# sc.tcells <- subset(sc.tcells, subset = orig.ident %in% names(count_few), invert = T)

##########2.Quality control##########
sc.tcells <- subset(sc.tcells, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
vln_plot <- VlnPlot(sc.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, 
                    raster = F, pt.size = 0)
ggsave(filename = "./results/violin_plot~before QC.png", plot = vln_plot, width = 10, height = 6)
table(sc.tcells@meta.data$orig.ident)

##########3-1Non-Integration##########
sc.tcells <- NormalizeData(sc.tcells)
sc.tcells <- FindVariableFeatures(sc.tcells)
sc.tcells <- ScaleData(sc.tcells)
sc.tcells <- RunPCA(sc.tcells)
sc.tcells <- FindNeighbors(sc.tcells, reduction = "pca", dims = 1:30)
sc.tcells <- FindClusters(sc.tcells, resolution = 0.1, cluster.name = "unintegrated_clusters_res=0.1")
sc.tcells <- FindClusters(sc.tcells, resolution = 0.8, cluster.name = "unintegrated_clusters_res=0.8")
sc.tcells <- RunUMAP(sc.tcells, reduction = "pca", dims = 1:30, reduction.name = "umap.unintegrated")
pdf("./results/dimplot_unintegrated.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.tcells, reduction = "umap.unintegrated", 
        group.by = c("orig.ident", "seurat_clusters"),
        raster = FALSE)
dev.off()

pdf("./results/dimplot_unintegrated_res=0.1.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.tcells, reduction = "umap.unintegrated", 
        group.by = c("orig.ident", "unintegrated_clusters_res=0.1"),
        raster = FALSE)
dev.off()

pdf("./results/dimplot_unintegrated_res=0.8.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.tcells, reduction = "umap.unintegrated", 
        group.by = c("orig.ident", "unintegrated_clusters_res=0.8"),
        raster = FALSE)
dev.off()

##########3-2CCA-Integration##########
sc.tcells[["RNA"]] <- split(sc.tcells[["RNA"]], f = sc.tcells$orig.ident)
sc.tcells <- NormalizeData(sc.tcells)
sc.tcells <- FindVariableFeatures(sc.tcells)
sc.tcells <- ScaleData(sc.tcells)
sc.tcells <- RunPCA(sc.tcells)
sc.tcells <- IntegrateLayers(sc.tcells,
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

pdf("./results/dimplot_integrated_cca_res=0.8.pdf", width = 15, height = 6, useDingbats = F)
DimPlot(sc.tcells, reduction = "umap.cca", group.by = c("orig.ident", "integrated_cca_clusters_res.0.8"),
        raster = FALSE, label = TRUE, label.size = 8)+ NoLegend()
dev.off()

saveRDS(sc.tcells, file = "./tmp/sc.tcells.after_umap.cca~without QC.rds")

##########3-3RPCA-Integration##########
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
