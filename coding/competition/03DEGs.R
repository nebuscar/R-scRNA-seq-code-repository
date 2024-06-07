library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)

##########1.Load data##########
seurat <- readRDS("./tmp/competition/sc.combined.after_umap.cca~without QC.rds")
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
de.01_02.markers <- FindMarkers(seurat, ident.1 = 1, ident.2 = 2, only.pos = F, min.pct = 0.5)
de.01_02.markers %>%
  dplyr::filter(avg_log2FC > 2) %>%
  filter(p_val_adj < 0.01) %>%
  # slice_head(n = 10) %>%
  ungroup() -> top10
write.csv(de.01_02.markers, "./results/competition/DEGs_res0.1_cluster_01.02.csv", 
          row.names = TRUE, quote = FALSE)
write.csv(top10, "./results/competition/DEGs_res0.1_top10gene_cluster_01.02.csv", 
          row.names = TRUE, quote = FALSE)

# setup threshold
log2FC = 2
padj = 0.01 

# Initial assignment
de.01_02.markers$threshold <- "ns"

# Define the conditions
up_condition <- which(de.01_02.markers$avg_log2FC > (log2FC) & de.01_02.markers$p_val_adj < padj)
down_condition <- which(de.01_02.markers$avg_log2FC < (-log2FC) & de.01_02.markers$p_val_adj < padj)

# Assign 'up' and 'down' thresholds
de.01_02.markers[up_condition,]$threshold <- "up"
de.01_02.markers[down_condition,]$threshold <- "down"

# Convert threshold to factor
de.01_02.markers$threshold <- factor(de.01_02.markers$threshold, levels = c('down', 'ns', 'up'))

p1 <- ggplot(data = de.01_02.markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = threshold)) +
  geom_point(alpha = 0.8, size = 0.8) +
  geom_vline(xintercept = c(-log2FC, log2FC), linetype = 2, color = "grey") +
  geom_hline(yintercept = -log10(padj), linetype = 2, color = "grey") +
  xlab(bquote(Log[2] * FoldChange)) +
  ylab(bquote(-Log[10] * italic(P.adj))) +
  theme_classic(base_size = 14) +
  scale_color_manual(
    '',
    labels = c(paste0("down(", table(de.01_02.markers$threshold)[[1]], ')'), 'ns', paste0("up(", table(de.01_02.markers$threshold)[[3]], ')')),
    values = c("blue", "grey", "red")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(
    text = element_text(size = 20),
    axis.title.x = element_text(size = 24, face = "bold"),  # Increase axis title size
    axis.title.y = element_text(size = 24, face = "bold"),  # Increase axis title size
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.text = element_text(size = 15)
  )

# Subset data for labels
label_data <- subset(de.01_02.markers, de.01_02.markers$p_val_adj < padj & abs(de.01_02.markers$avg_log2FC) >= log2FC)

# Add labels to the plot
p2 <- p1 + ggrepel::geom_text_repel(
  data = label_data,
  aes(label = rownames(label_data)), 
  size = 5,
  box.padding = unit(0.5, "lines"),
  point.padding = unit(0.8, "lines"), 
  segment.color = "black", 
  show.legend = FALSE,
  max.overlaps = 10
) +
  theme(
    text = element_text(size = 20),
    axis.title.x = element_text(size = 24, face = "bold"),  # Increase axis title size
    axis.title.y = element_text(size = 24, face = "bold"),  # Increase axis title size
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.text = element_text(size = 15)
  )
pdf(file = "./results/competition/figure2/DEGs_plot_res0.1_IntegrateLayers_cca.pdf", width = 10, height = 10)
print(p2)
dev.off()

cluster.01_02.de.genes <- rownames(de.01_02.markers)
write.csv(cluster.01_02.de.genes, file = "./results/competition/DEGs_res01_cluster_01.02.csv")

