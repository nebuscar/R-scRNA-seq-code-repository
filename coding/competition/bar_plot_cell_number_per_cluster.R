library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)

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
  scale_fill_manual(values = c("Ctrl" = "red", "STM" = "green")) +
  labs(title = "Cell Counts by Recovery Status Across Clusters",
       x = "Cluster", y = "Cell number") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),  
    plot.background = element_rect(fill = "white"),    
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "right"
  )
ggsave(filename = "./results/competition/bar_plot_cluster_cell_counts_cca.png", plot = bar_plot, width = 10, height = 6)