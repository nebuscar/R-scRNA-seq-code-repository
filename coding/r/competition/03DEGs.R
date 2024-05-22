library(Seurat)
library(ggplot2)
library(dplyr)

##########1.Load data##########
seurat <- readRDS("./tmp/competition/sc.combined.after_umap.rpca~without QC.rds")
levels(seurat)

##########find markers, report only the positive##########
de.markers <- FindAllMarkers(seurat, only.pos = T)
de.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)
de.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
write.csv(de.markers, "./results/competition/res01_findallmarkers.csv", 
          row.names = TRUE, quote = FALSE)
write.csv(top10, "./results/competition/res01_top10gene.csv", 
          row.names = TRUE, quote = FALSE)
