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
  theme_minimal()+
  ggtitle("Gene Expression Violin Plot") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.title = element_blank()
  )
print(vln_plot)

ggsave(filename = "./results/competition/violin_plot_gene_expression.png", plot = vln_plot, width = 10, height = 6)