library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(Seurat)

sc.combined <- readRDS("./tmp/competition/sc.combined.after_umap.cca~without QC.rds")
# cell number in each cluster
recovery_status <- ifelse(sc.combined$orig.ident %in% c("Ctrl1", "Ctrl2"), "Ctrl", "STM")
sc.combined <- AddMetaData(sc.combined, metadata = recovery_status, col.name = "recovery_status")

cluster_counts <- table(Idents(sc.combined), sc.combined$recovery_status)
cluster_counts_df <- as.data.frame(cluster_counts)
names(cluster_counts_df) <- c("cluster", "recovery_status", "cell_count")

cluster_counts_wide <- reshape(cluster_counts_df, idvar = "cluster", timevar = "recovery_status", direction = "wide")
names(cluster_counts_wide) <- c("cluster", "Ctrl", "STM")
cluster_counts_wide[is.na(cluster_counts_wide)] <- 0

cluster_counts_long <- reshape::melt(cluster_counts_wide, id.vars = "cluster", variable.name = "recovery_status", value.name = "cell_count")
names(cluster_counts_long) <- c("cluster", "recovery_status", "cell_count")

bar_plot <- ggplot(cluster_counts_long, aes(x = factor(cluster), y = cell_count, fill = recovery_status)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Ctrl" = "#008c54", "STM" = "#ffa100")) +
  labs(title = "Cell Counts by Recovery Status Across Clusters", x = "Cluster", y = "Cell number") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 30), 
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    legend.position = "right",
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 30)
  ) +
  guides(fill = guide_legend(title = "Recovery Status", title.position = "top", title.theme = element_text(size = 30)))

pdf(file = "./results/competition/bar_plot_cluster_cell_counts_cca.pdf", width = 10, height = 6)
print(bar_plot)
dev.off()

# vln_plot
genes_of_interest <- c("CD2", "CD3G", "CD86")
expression_data <- FetchData(sc.combined, vars = c(genes_of_interest, "seurat_clusters"))
expression_data_long <- melt(expression_data, id.vars = "seurat_clusters", variable.name = "gene", value.name = "expression")
custom_colors <- c("#c9a1f2", "#1e77b3", "#6a3d9a", "#f47a79")

cluster_count <- length(unique(expression_data_long$seurat_clusters))
vln_plot <- ggplot(expression_data_long, aes(x = factor(seurat_clusters), y = expression, fill = factor(seurat_clusters))) +
  geom_violin(scale = "width", trim = TRUE) +
  facet_wrap(~ gene, scales = "free_y") +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Seurat Clusters", y = "Expression") +
  theme(
    text = element_text(size = 30),
    panel.background = element_rect(fill = "white", color = NA), 
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.title = element_blank()
  )
pdf(file = "./results/competition/violin_plot_gene_expression.pdf", width = 10, height = 6)
print(vln_plot)
dev.off()
