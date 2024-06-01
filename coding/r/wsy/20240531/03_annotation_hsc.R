library(Seurat)
library(ggplot2)
library(dplyr)

# markers
markers <- readxl::read_excel("./BMMC marker.xls")
markers.list <- sapply(markers$Marker, function(x) unlist(strsplit(x, ","), use.names=FALSE))
markers.list <- unlist(markers.list, use.names=FALSE)
markers.list <- unique(sort(markers.list))

##HSC
bm <- readRDS("./tmp/wsy/20240531/sk.bone.rds")
bm <- FindClusters(bm, resolution = 0.5)
bm <- RunUMAP(bm, dims = 1:30, reduction = "integrated.cca")

hsc <- subset(bm, idents = c("1", "9", "10"))
hsc <- FindClusters(hsc, resolution = 1.3)
hsc <- RunUMAP(hsc, dims = 1:30, reduction = "integrated.cca")
#rm(bm)

hsc.markers <- FindAllMarkers(hsc)
hsc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

# Markers & ident
geneset <- c("MEIS1", "MSI2", "CD34", "HOPX", "EGR1", "PF4", "PPBP", "ITGA2B", "GATA1",
             "GYPA", "PECAM1", "CD79A", "TFRC", "BLVRB", "SPTA1", "LYZ", "MPO", "ELANE",
             "AZU1", "CSF3R", "AVP", "HLF", "SLAMF1", "FLT3")
# tIdent <- c("HSC", "CD34−_HSC", "EPC", "NPC", "pro_B", "pre_pro_B", "CLP", "MLP",
#             "myeloid−biased_HSC", "CMP", "MEP", "GMP")
tIdent <- c("ST-HSC", "pro_B", "LT-HSC", "CLP", "NPC", "MLP", "EPC", "MEP", "MPP3", "pre_pro_B", "CMP", "MPP4")
names(tIdent) <- levels(hsc)
hsc <- RenameIdents(hsc, tIdent)

# Visualization
p1 <- DimPlot(hsc, reduction = "umap", label = TRUE, label.size = 6)+
  NoLegend()
p2 <- DotPlot(hsc, features = geneset)+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 8),
        axis.text.y = element_text(angle = 0, hjust=1, size = 8),
        axis.title = element_blank())+
  #RotatedAxis()+
  #coord_flip()

pdf("./results/wsy/20240531/hsc_13.pdf", width = 14, height = 8, useDingbats = F)
p1+p2
dev.off()

# PCT
hsc@meta.data$seurat_clusters <- Idents(hsc)
sk.pct <- table(hsc@meta.data[, c("seurat_clusters", "orig.ident")])
sk.pct <- as.data.frame(sk.pct)
sk.pct <- sk.pct[-2,]

sk.row.pct <- reshape(sk.pct, idvar = "orig.ident", timevar = "seurat_clusters", direction = "wide")
sk.row.pct <- sk.row.pct[, -ncol(sk.row.pct)]
row.names(sk.row.pct) <- sk.row.pct[,1]
sk.row.pct <- sk.row.pct[, -1]
sk.row.pct <- sk.row.pct/rowSums(sk.row.pct)*100
sk.row.pct <- round(sk.row.pct, digits = 2)

##########Group##########
hsc_origIdent <- hsc@meta.data$orig.ident
levels(hsc_origIdent) <- c("GBM1", "GBM2", "PD", "PD", "PD", "PD", "ST1", "TBI", "TBI")
hsc@meta.data$group <- hsc_origIdent
hsc_pd_tbi <- subset(hsc, subset = group %in% c("PD", "TBI"))


marker_list <- list()
for (cluster in levels(hsc_pd_tbi)){
  print(cluster)
  tmp <- subset(hsc_pd_tbi, idents = cluster)
  Idents(tmp) <- tmp@meta.data$group
  markers <- FindMarkers(tmp, ident.1 = "PD", ident.2 = "TBI")
  marker_list[[cluster]] <- markers %>%
    dplyr::filter(avg_log2FC > 1) %>%
    filter(p_val_adj < 0.05) %>%
    ungroup()
  marker_list[[cluster]]$cluster <- cluster
}
markers <- do.call(rbind.data.frame, marker_list)
write.csv(markers, "./results/wsy/20240531/markers_hsc_res13_PB_TBI.csv")
