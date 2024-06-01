library(Seurat)
library(ggplot2)
library(dplyr)

# markers
markers <- readxl::read_excel("./BMMC marker.xls")
markers.list <- sapply(markers$Marker, function(x) unlist(strsplit(x, ","), use.names=FALSE))
markers.list <- unlist(markers.list, use.names=FALSE)
markers.list <- unique(sort(markers.list))

##HSC
sk.tcell <- readRDS("./tmp/wsy/20240531/sk.tcell.rds")
sk.tcell <- FindClusters(sk.tcell, resolution = 0.5)
sk.tcell <- RunUMAP(sk.tcell, dims = 1:30, reduction = "integrated.cca")

nk <- subset(sk.tcell, idents = c("4"))
nk <- FindClusters(nk, resolution = 0.8)
nk <- RunUMAP(nk, dims = 1:30, reduction = "integrated.cca")
#rm(bm)

nk.markers <- FindAllMarkers(nk)
nk.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

# Markers & ident
# geneset <- c("TBX21", "EOMES", "CD5", "CD200R1", "IL17RB", "IL1R1", "RORC", "NRP1", "LIF")

geneset <- c("NCAM1", "FCGR3A",
             "IL32", 
             "IFI44L","CX3CR1",
             "BEND7", "ZNF90", "NACA2",
             "CCL3", "NFKBIA",
             "SPC25", "MKI67",
             "HSPA1A", "DNAJB1", 
             "CREM", "NR4A3",
             "GCSAM", "KLRC2",
             "CCL4", "CD160")
# tIdent <- c("CD56dimCD16hi-NFKBIA", "CD56dimCD16hi-CX3CR1",
#             "CD56dimCD16hi-IL32", "CD56dimCD16lo", "CD56brightCD16lo-CCL3") 
tIdent <- c("CD56dimCD16lo","CD56dimCD16hi−IL32", "CD56dimCD16hi-CX3CR1", 
            "Unknow", "CD56brightCD16lo-CCL3", "CD56dimCD16hi-NFKBIA")
names(tIdent) <- levels(nk)
nk <- RenameIdents(nk, tIdent)

# Visualization
p1 <- DimPlot(nk, reduction = "umap", label = TRUE, label.size = 5)+
  NoLegend()
p2 <- DotPlot(nk, features = geneset)+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 8),
        axis.text.y = element_text(angle = 0, hjust=1, size = 8),
        axis.title = element_blank())+
  #RotatedAxis()+
  #coord_flip()

pdf("./results/wsy/20240531/nk_08.pdf", width = 14, height = 8, useDingbats = F)
p1+p2
dev.off()

nk@meta.data$seurat_clusters <- Idents(nk)
sk.pct <- table(nk@meta.data[, c("seurat_clusters", "orig.ident")])
sk.pct <- as.data.frame(sk.pct)

sk.row.pct <- reshape(sk.pct, idvar = "orig.ident", timevar = "seurat_clusters", direction = "wide")
sk.row.pct <- sk.row.pct[, -ncol(sk.row.pct)]
row.names(sk.row.pct) <- sk.row.pct[,1]
sk.row.pct <- sk.row.pct[, -1]
sk.row.pct <- sk.row.pct/rowSums(sk.row.pct)*100
sk.row.pct <- round(sk.row.pct, digits = 2)

##########Group##########
nk_origIdent <- nk@meta.data$orig.ident
levels(nk_origIdent) <- c("GBM1", "GBM2", "PD", "PD", "PD", "PD", "ST1", "TBI", "TBI")
nk@meta.data$group <- nk_origIdent
nk_pd_tbi <- subset(nk, subset = group %in% c("PD", "TBI"))


marker_list <- list()
for (cluster in levels(nk_pd_tbi)){
  print(cluster)
  tmp <- subset(nk_pd_tbi, idents = cluster)
  Idents(tmp) <- tmp@meta.data$group
  markers <- FindMarkers(tmp, ident.1 = "PD", ident.2 = "TBI")
  marker_list[[cluster]] <- markers %>%
    dplyr::filter(avg_log2FC > 1) %>%
    filter(p_val_adj < 0.05) %>%
    ungroup()
  marker_list[[cluster]]$cluster <- cluster
}
markers <- do.call(rbind.data.frame, marker_list)
write.csv(markers, "./results/wsy/20240531/markers_nk_res8_PB_TBI.csv")
