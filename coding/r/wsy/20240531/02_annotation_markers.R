library(Seurat)
library(ggplot2)
library(dplyr)

# markers
markers <- readxl::read_excel("./BMMC marker.xls")
markers.list <- sapply(markers$Marker, function(x) unlist(strsplit(x, ","), use.names=FALSE))
markers.list <- unlist(markers.list, use.names=FALSE)
markers.list <- unique(sort(markers.list))

DotPlot(sk.scData, group.by = "RNA_snn_res.0.1",
        features = markers.list)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
DimPlot(sk.scData, group.by = "RNA_snn_res.0.1",
        reduction = "umap", raster = FALSE, 
        label = TRUE, label.size = 4)+ NoLegend()
# further reduction
sk.scData <- readRDS("./tmp/wsy/20240531/sk.scData.rds")
sk.tcell <- subset(sk.scData, subset = RNA_snn_res.0.1 %in% c("0", "11"))
sk.bone <- subset(sk.scData, subset = RNA_snn_res.0.1 %in% c("6", "5", "7", "10", "9", "12"))
sk.mono <- subset(sk.scData, subset = RNA_snn_res.0.1 %in% c("1", "2", "3", "4"))
saveRDS(sk.tcell, "./tmp/wsy/20240531/sk.tcell.rds")
saveRDS(sk.bone, "./tmp/wsy/20240531/sk.bone.rds")
saveRDS(sk.mono, "./tmp/wsy/20240531/sk.mono.rds")

##########T cell##########
#future::multisession(workers = 8)
sk.tcell <- readRDS("./tmp/wsy/20240531/sk.tcell.rds")
sk.tcell <- FindClusters(sk.tcell, resolution = 0.5)
sk.tcell <- RunUMAP(sk.tcell, dims = 1:30, reduction = "integrated.cca")

geneset <- c("CD2", "CD3D", "TRBC2", "TRBC1",
             "CD4", "CCR7", "LEF1", "TCF7", "SELL",
             "CD8A", "CD8B", "GZMA", "LAG3", "TIGIT",
             "FCGR3A", "HOPX", "KLRD1", "KLRC1", "NKG7", 
             "LYZ", "MPO", "MSI2")
# tIdent <- c("CD4+ naive T", "CD4+ effector/memory", "CD8+ exhuasted T", 
#             "CD8+ effector T", "NK", "NK-T", "Others") 
tIdent <- c("CD4+ naive T", "CD8+ exhuasted T", 
            "CD4+ effector/memory","CD8+ effector T",
            "NK-1", "NK-T",  "Others", "CD4+ naive T-2", "NK-2") 
names(tIdent) <- levels(sk.tcell)
sk.tcell <- RenameIdents(sk.tcell, tIdent)

p1 <- DimPlot(sk.tcell, reduction = "umap", label = TRUE, label.size = 5)
  # NoLegend()
p2 <- DotPlot(sk.tcell,
        features = geneset)+
  theme(axis.text.x = element_text(angle = 90, hjust=1, size = 8),
        axis.text.y = element_text(angle = 0, hjust=1, size = 8),
        axis.title = element_blank())+
  scale_y_discrete(limits=rev)
  # NoLegend()

pdf("./results/wsy/20240531/tcell_05.pdf", width = 14, height = 8, useDingbats = F)
p1+p2
dev.off()

sk.tcell@meta.data$seurat_clusters <- Idents(sk.tcell)
sk.pct <- table(sk.tcell@mewta.data[, c("seurat_clusters", "orig.ident")])
sk.pct <- as.data.frame(sk.pct)

sk.row.pct <- reshape(sk.pct, idvar = "orig.ident", timevar = "seurat_clusters", direction = "wide")
row.names(sk.row.pct) <- sk.row.pct[,1]
sk.row.pct <- sk.row.pct[, -1]
sk.row.pct <- sk.row.pct/rowSums(sk.row.pct)*100
sk.row.pct <- round(sk.row.pct, digits = 2)

pdf("./results/wsy/20240531/pct_tcell.pdf", width = 8, height = 8, useDingbats = F)
ggplot(sk.pct, aes(fill = orig.ident, y = Freq, x = seurat_clusters)) + 
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 9),
        axis.text.y = element_text(angle = 0, hjust=1, size = 9),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
        )+
  ggsci::scale_fill_npg(name = "")+
  labs(y = "Cell count", x = "")
dev.off()

##########myeloid##########
myeloid <- readRDS("./tmp/wsy/20240531/sk.mono.rds")
myeloid <- FindClusters(myeloid, resolution = 0.5)
myeloid <- RunUMAP(myeloid, dims = 1:30, reduction = "integrated.cca")

# Markers & ident
geneset <- c("SELL", "CXCR2", 
             "DEFA4", "ELANE", "AZU1", "CD177", "EGR1", "MPO",
             "CD74", "CD1C", "FCGR3A", "CD14", "FCER1A", "LILRB4", "NKG7", "VCAM1", "CD163", "ITGAM",
             "TRBC2", "CD3D")
# tIdent <- c("NDN_4", "NDN_3", "LDN_3", "NDN_1", "LDN_1", "LDN_2", "Monocyte",
#             "NDN_2", "LDN_4", "Doublet", "NDN_5", "LDN_5", "DC", "Doublet_2", "macrophage") 
tIdent <- c("NDN_4", "NDN_3", "LDN_1", "LDN_2", "LDN_3", "NDN_1", "NDN_2", 
            "Monocyte", "LDN_4", "other", "NDN_5","DC", "Doublet") 
names(tIdent) <- levels(myeloid)
myeloid <- RenameIdents(myeloid, tIdent)

# Visualization
p1 <- DimPlot(myeloid, reduction = "umap", label = TRUE, label.size = 5)+
  NoLegend()
p2 <- DotPlot(myeloid, features = geneset)+
  theme(axis.text.x = element_text(angle = 90, hjust=1, size = 8),
        axis.text.y = element_text(angle = 0, hjust=1, size = 8),
        axis.title = element_blank())+
  #RotatedAxis()+
  #coord_flip()+
  # scale_y_discrete(limits = c("NDN_1", "NDN_2", "NDN_3", "NDN_4", "NDN_5",
  #                             "LDN_1", "LDN_2", "LDN_3", "LDN_4", "LDN_5",
  #                             "Monocyte", "Doublet", "macrophage", "DC", "Doublet_2"))+
  # NoLegend()

pdf("./results/wsy/20240531/myeloid_05.pdf", width = 14, height = 8, useDingbats = F)
p1+p2
dev.off()

myeloid@meta.data$seurat_clusters <- Idents(myeloid)
sk.pct <- table(myeloid@meta.data[, c("seurat_clusters", "orig.ident")])
sk.pct <- as.data.frame(sk.pct)

sk.row.pct <- reshape(sk.pct, idvar = "orig.ident", timevar = "seurat_clusters", direction = "wide")
row.names(sk.row.pct) <- sk.row.pct[,1]
sk.row.pct <- sk.row.pct[, -1]
sk.row.pct <- sk.row.pct/rowSums(sk.row.pct)*100
sk.row.pct <- round(sk.row.pct, digits = 2)

pdf("./result/subset_pct/pct_nm.pdf", width = 8, height = 8, useDingbats = F)
ggplot(sk.pct, aes(fill = orig.ident, y = Freq, x = seurat_clusters)) + 
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 9),
        axis.text.y = element_text(angle = 0, hjust=1, size = 9),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
  )+
  ggsci::scale_fill_npg(name = "")+
  labs(y = "Cell count", x = "")
dev.off()


##########BM##########
bm <- readRDS("./tmp/wsy/20240531/sk.bone.rds")
bm <- FindClusters(bm, resolution = 0.5)
bm <- RunUMAP(bm, dims = 1:30, reduction = "integrated.cca")

# Markers & ident
geneset <- c("CD34", "FLT3", "CD79A", "CD19", "SOX4", "CD24", "PAX5", "EBF1", "MS4A1", "IGHM", "IGHD", "IGHG1", "MZB1", "NKG7", "VPREB1", 
             "CSF3R", "EGR1", "HOPX", "MEIS1", "MSI2", "SELL", "SPTA1", "TFRC",
             "BLVRB", "GATA1", "GYPA",
             "CLEC4C", "FCER1A", "IL3RA", "MS4A2")
tIdent <- c("mature B", "HSC", "ETD", "CEP_2", "plasma", "CEP_1", "pDC",
            "immature B", "red_cell", "EEP", "pro_B", "basophil", "NKB", "trans B_2",
            "trans B_1") 
names(tIdent) <- levels(bm)
bm <- RenameIdents(bm, tIdent)

# Visualization
p1 <- DimPlot(bm, reduction = "umap", label = TRUE, label.size = 4)+
  NoLegend()
  # scico::scale_color_scico_d(palette = "managua")
p2 <- DotPlot(bm, features = geneset)+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 8),
        axis.text.y = element_text(angle = 0, hjust=1, size = 8),
        axis.title = element_blank())+
  #RotatedAxis()+
  #coord_flip()+
  NoLegend()+
  # scale_y_discrete(limits = c("pro_B", "immature B", "trans B_1", "trans B_2", "mature B", "plasma", "NKB",
  #                             "HSC", "EEP", "CEP_1", "CEP_2", "ETD", "red_cell",
  #                             "basophil", "pDC"))

pdf("./results/wsy/20240531/bm_05.pdf", width = 14, height = 8, useDingbats = F)
p1+p2
dev.off()

bm@meta.data$seurat_clusters <- Idents(bm)
sk.pct <- table(bm@meta.data[, c("seurat_clusters", "orig.ident")])
sk.pct <- as.data.frame(sk.pct)

sk.row.pct <- reshape(sk.pct, idvar = "orig.ident", timevar = "seurat_clusters", direction = "wide")
row.names(sk.row.pct) <- sk.row.pct[,1]
sk.row.pct <- sk.row.pct[, -1]
sk.row.pct <- sk.row.pct/rowSums(sk.row.pct)*100
sk.row.pct <- round(sk.row.pct, digits = 2)

pdf("./results/wsy/20240531/pct_bone.pdf", width = 8, height = 8, useDingbats = F)
ggplot(sk.pct, aes(fill = orig.ident, y = Freq, x = seurat_clusters)) + 
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 9),
        axis.text.y = element_text(angle = 0, hjust=1, size = 9),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
  )+
  ggsci::scale_fill_npg(name = "")+
  labs(y = "Cell count", x = "")
dev.off()
