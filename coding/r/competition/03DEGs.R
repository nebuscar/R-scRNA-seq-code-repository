library(Seurat)
library(ggplot2)

##########1.Load data##########
seurat <- readRDS("./tmp/competition/sc.combined.after_umap.rpca~without QC.rds")
levels(seurat)

##########
de.markers <- FindAllMarkers(seurat)
