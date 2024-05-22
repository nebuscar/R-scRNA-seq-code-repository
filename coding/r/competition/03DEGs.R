library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)

##########1.Load data##########
seurat <- readRDS("./tmp/competition/sc.combined.after_umap.rpca~without QC.rds")
levels(seurat)

##########find all markers, report only the positive##########
de.markers <- FindAllMarkers(seurat, only.pos = T)
de.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
write.csv(de.markers, "./results/competition/res01_findallmarkers.csv", 
          row.names = TRUE, quote = FALSE)
write.csv(top10, "./results/competition/res01_top10gene.csv", 
          row.names = TRUE, quote = FALSE)

##########find markers##########
de.01_02.markers <- FindMarkers(seurat, ident.1 = 1, ident.2 = 2, only.pos = F, min.pct = 0.05)
de.01_02.markers %>%
  dplyr::filter(avg_log2FC > 2) %>%
  filter(p_val_adj < 0.01) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
write.csv(de.01_02.markers, "./results/competition/DEGs_res01_cluster_01.02.csv", 
          row.names = TRUE, quote = FALSE)
write.csv(top10, "./results/competition/DEGs_res01_top10gene_cluster_01.02.csv", 
          row.names = TRUE, quote = FALSE)

# setup threshold
log2FC = 2
padj = 0.01 

# Initial assignment
de.01_02.markers$threshold <- "ns"

# Define the conditions
up_condition <- which(de.01_02.markers$avg_log2FC > (log2FC) & de.01_02.markers$p_val_adj < padj)
down_condition <- which(de.01_02.markers$avg_log2FC < (-log2FC) & de.01_02.markers$p_val_adj < padj)

# Debugging: Print the number of rows meeting the conditions
cat("Number of rows for 'up' condition:", length(up_condition), "\n")
cat("Number of rows for 'down' condition:", length(down_condition), "\n")

# Assign 'up' and 'down' thresholds
de.01_02.markers[up_condition,]$threshold <- "up"
de.01_02.markers[down_condition,]$threshold <- "down"

# Convert threshold to factor
de.01_02.markers$threshold <- factor(de.01_02.markers$threshold, levels = c('down', 'ns', 'up'))

p1 <- ggplot(data=de.01_02.markers, aes(x=avg_log2FC, y=-log10(p_val_adj), color=threshold)) +
  geom_point(alpha=0.8, size=0.8) +
  geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color="grey")+
  geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+
  #labs(title= ifelse(""==title, "", paste("DEG:", title)))+
  xlab(bquote(Log[2]*FoldChange))+
  ylab(bquote(-Log[10]*italic(P.adj)) )+
  theme_classic(base_size = 14) +
  scale_color_manual('',labels=c(paste0("down(",table(de.01_02.markers$threshold)[[1]],')'),'ns',
                                 paste0("up(",table(de.01_02.markers$threshold)[[3]],')' )),
                     values=c("blue", "grey","red" ) )+
  guides(color=guide_legend(override.aes = list(size=3, alpha=1)))

# Subset data for labels
label_data <- subset(de.01_02.markers, de.01_02.markers$p_val_adj < padj & abs(de.01_02.markers$avg_log2FC) >= log2FC)
p2 <- p1 + geom_text_repel(
  data = label_data,
  aes(label = rownames(label_data)), 
  size = 3,
  box.padding = unit(0.5, "lines"),
  point.padding = unit(0.8, "lines"), 
  segment.color = "black", 
  show.legend = FALSE
)

cluster.01_02.de.genes <- rownames(de.01_02.markers)
write.csv(cluster.01_02.de.genes, file = "./results/competition/DEGs_res01_cluster_01.02.csv")
