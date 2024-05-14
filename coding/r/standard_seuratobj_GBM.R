library(Seurat)

###########################################
tcells <- readRDS("data/seurat/tcells_v5.rds")
count_dataset <- table(tcells@meta.data$orig.ident)
count_few <- count_dataset[count_dataset < 30]
tcells <- subset(tcells, subset = orig.ident %in% names(count_few), invert = T)
#tcells[["RNA5"]] <- as(object = tcells[["RNA"]], Class = "Assay5")
#DefaultAssay(tcells)='RNA5'
#tcells[['RNA']]=NULL  
tcells[["RNA"]] <- split(tcells[["RNA"]], f = tcells$orig.ident)
tcells <- NormalizeData(tcells)
tcells <- FindVariableFeatures(tcells)
tcells <- ScaleData(tcells)
tcells <- RunPCA(tcells)

##################CCA.integrateLayers##################
tcells <- IntegrateLayers(object = tcells, 
                          method = RPCAIntegration,
                          orig.reduction = "pca", 
                          new.reduction = "integrated.rpca",
                          verbose = F, 
                          k.weight = 32)
saveRDS(tcells,file = "data/seurat/integrate_rpca.rds")

tcells <- readRDS("data/seurat/integrate_rpca.rds")
tcells <- FindNeighbors(tcells, reduction = "integrated.rpca", dims = 1:30)
tcells <- FindClusters(tcells, resolution = 0.8, cluster.name = "rpca_clusters")
tcells <- RunUMAP(tcells, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

p1 <- DimPlot(tcells, reduction = "umap.rpca", group.by = "orig.ident", label.size = 2)
p2 <- DimPlot(tcells, reduction = "umap.rpca", group.by = "rpca_clusters", label.size = 2) 
p1 + p2

##################add meta.scissor##################
scissor <- readRDS("data/downstream/sc_scissor.rds")
scissor_sub <- scissor@meta.data
scissor_sub <- scissor_sub[match(colnames(tcells), rownames(scissor_sub)),]
tcells@meta.data$scissor <- scissor_sub$scissor

plot3 <- DimPlot(tcells,group.by = "scissor",cols = c('grey','royalblue','indianred1'), size = 1)
plot1+plot2+plot3


saveRDS(tcells,file = "data/seurat/tcells_rm_batch-1.rds")

