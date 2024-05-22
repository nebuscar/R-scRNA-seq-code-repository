library(dplyr)
library(ggplot2)
library(Seurat)
devtools::install_github("chuiqin/irGSEA")

sc.nash <- readRDS("./data/sc.sub.rds")
scScore <- data.frame(ucell = sc.nash[["UCell"]]$data, 
                      singscore = sc.nash[["singscore"]]$data, 
                      cluster = sc.nash@meta.data$seurat_clusters)
rn <- do.call(rbind, strsplit(rownames(scScore), '_'))
scScore <- cbind(scScore, rn)
colnames(scScore) <- c("ucell", "singscore", "cluster", "type", "time", "indiv", "capture", "barcode")
scScore$group <- paste(scScore$type, scScore$time, scScore$indiv)

# 
scScore.mean <- scScore %>%
  group_by(type, time, cluster, indiv) %>%
  mutate(Ucell = as.numeric(ucell),
         singscore = as.numeric(singscore)) %>%
  summarise(mean_ucell = mean(ucell),
            mean_singscore = mean(singscore),
            sd_ucell = sd(ucell),
            sd_singscore = sd(singscore))

pdf("./results/wsy/voilin_score_group.pdf", width = 12, height = 6, useDingbats = F)
scScore %>%
  filter(cluster == 0) %>%
ggplot(aes(x = group, y = as.numeric(ucell), color = group)) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.1)+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 8),
        axis.text.y = element_text(angle = 0, hjust=1, size = 8),
        axis.title = element_blank(),
        legend.position="none",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()

scScore.stat <- scScore %>%
  filter(time %in% c("30weeks", "15weeks"), cluster == 0, type == "Hepatocyte") %>%
  mutate(ucell = as.numeric(ucell),
         singscore = as.numeric(singscore)) %>%
  summarise(p = wilcox.test(ucell~time)$p.value)


##########
idents.0.4 <- subset(sc.nash, subset = seurat_clusters %in% c("0", "4"))
idents.0.4 <- FindClusters(idents.0.4, resolution = 0.5, cluster.name = "res.0.5")
idents.0.4 <- RunUMAP(idents.0.4, dims = 1:30, reduction = "integrated.cca")
pdf("./results/wsy/dimplot_res=0.5.pdf", width = 15, height = 10, useDingbats = F)
p1 <- DimPlot(idents.0.4, reduction = "umap", group.by = "seurat_clusters",
        raster = FALSE, label = TRUE, label.size = 6)+ NoLegend()
p2 <- DimPlot(idents.0.4, reduction = "umap", group.by = "orig.ident",
              raster = FALSE, label = TRUE, label.size = 6)+ NoLegend()
p1 + p2
dev.off()

scScore_sub <- scScore[match(rownames(sc.nash@meta.data), rownames(scScore)), ]
sc.nash@meta.data$scScore_sub <- scScore_sub
tmp <- sc.nash@reductions$umap@cell.embeddings
tmp <- tmp[match(rownames(scScore_sub), rownames(tmp)), ]
scScore_sub <- cbind(scScore_sub,tmp)

p3 <- ggplot(scScore_sub,aes(x = umap_1, y = umap_2, color = ucell)) + 
  geom_point()+
  theme(axis.text.x = element_text(angle = 0, hjust=1, size = 16),
        axis.text.y = element_text(angle = 0, hjust=1, size = 16),
        axis.title = element_blank(),
        legend.position="none",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+scale_colour_gradient2(
          midpoint=0.15
          #limits=range(scScore$ucell)
        )
pdf("./results/wsy/signal_umap.pdf", width = 24, height = 12, useDingbats = F)
p1+p3

dev.off()

##########
scScore_sub$cluster <- sc.nash@meta.data$seurat_clusters
pdf("./results/wsy/voilin_score_sub_cluster_singscore.pdf", width = 12, height = 6, useDingbats = F)
scScore_sub %>%
  ggplot(aes(x = cluster, y = as.numeric(ucell), color = cluster)) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.1)+
  theme(axis.text.x = element_text(angle = 0, hjust=1, size = 8),
        axis.text.y = element_text(angle = 0, hjust=1, size = 8),
        axis.title = element_blank(),
        legend.position="none",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()

scScore_sub$group <- paste(scScore_sub$type, scScore_sub$time, scScore_sub$indiv)

scScore_pct <- as.data.frame(table(scScore_sub[c("group", "cluster")]))
pct <- reshape(scScore_pct, idvar = "group", timevar = "cluster", direction = "wide")

# pct <- pct[, -ncol(pct)]
row.names(pct) <- pct[,1]
pct <- pct[, -1]
pct <- pct/rowSums(pct)*100
pct <- round(pct, digits = 2)
write.table(pct, "./results/wsy/cluster_pct.csv", sep = ",", quote = F, col.names = NA)
