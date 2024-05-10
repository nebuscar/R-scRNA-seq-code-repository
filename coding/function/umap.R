library(dplyr)
library(monocle3)
library(Seurat)

setwd(".")
mtx <- readRDS("dataset4.rds")
mtx <- FindClusters(mtx)
mtx <- RunUMAP(mtx, dims = 1:30, reduction = "integrated.dr")

pdf("../results/dimplot.pdf", width = 20, height = 10, useDingbats = F)
pal <- as.data.frame(ggthemes::ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]$`Tableau 20`)[,2]
plot1 <- DimPlot(
  mtx, reduction = "umap", group.by = c("seurat_clusters"),
  pt.size = 1,
  combine = TRUE, label.size = 4,
  cols = pal
)
plot2 <- DimPlot(
  mtx, reduction = "umap", group.by = c("dataset"),
  pt.size = 1,
  combine = TRUE, label.size = 4,
  cols = c("#A4D38E", "#5D9BBE", "#E45A5F", "#F58135")
)
gridExtra::grid.arrange(plot1, plot2, ncol=2)
dev.off()

#umap <- cbind(mtx@reductions[["umap"]]@cell.embeddings, mtx@meta.data)
#ggplot(umap, 
#       aes(x = umap_1, y = umap_2, colour = dataset)) +
#  geom_point(size = 1, alpha = 0.5) +
#  ggthemes::scale_color_tableau() + theme_classic() +
#  xlab("PC1") + ylab("Timepoint") +
#  ggtitle("Cells ordered by first principal component")

# SCE
data <- GetAssayData(mtx, assay = "SCT", slot = "count")
cell_metadata <- mtx@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 10)

cds <- reduce_dimension(cds)
SingleCellExperiment::reducedDims(cds)[["UMAP"]] <- mtx@reductions[["umap"]]@cell.embeddings
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

get_earliest_principal_node <- function(cds, time_bin){
  cell_ids <- which(colData(cds)[, "dataset"] == time_bin)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] 
  root_pr_nodes
}

cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds, "293T"))

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

# DEG
track_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 6)
track <- track_genes[track_genes[, 6] < 1e-3, c(5,2,3,4,1,6)]
track_gene_top6 <- track %>% top_n(n = 6, morans_I) %>%
  pull(gene_short_name) %>% as.character()
plot_genes_in_pseudotime(cds[track_gene_top6,], color_cells_by = "dataset", 
                         min_expr=0.5, ncol = 3)

# subsample control cells
control.cells <- row.names(mtx@meta.data[mtx@meta.data[,9] == "293T" ,])
set.seed(seed = 14412)
sample.cells <- sample(x = 1:length(control.cells), size = 1000, replace = FALSE)
control.cells.remove <- control.cells[-sample.cells]
control.cells.keep <- row.names(mtx@meta.data)[-which(Cells(mtx) %in% control.cells.remove)]
mtx@meta.data[["cellID"]] <- Cells(mtx)
mtx.sub <- subset(mtx, subset = cellID %in% control.cells.keep)



subset1 <- Cells(subset(mtx.sub, subset = EGR1 > 1))
subset2 <- Cells(subset(mtx.sub, subset = IER3 > 1))
subset3 <- Cells(subset(mtx.sub, subset = JUND > 1))
subset4 <- Cells(subset(mtx.sub, subset = LIMD2 > 1))
subset5 <- Cells(subset(mtx.sub, subset = TNFAIP3 > 1))
subset6 <- Cells(subset(mtx.sub, subset = TRBC1 > 1))

subset <- c(subset1,subset2,subset3,subset4,subset5,subset6)
subset.keep <- table(subset)[table(subset)>1]
mtx.sub.sub <- subset(mtx.sub, subset = cellID %in% names(subset.keep))

pdf("../results/heatplot_identifocation.pdf", width = 10, height = 6, useDingbats = F)
DoHeatmap(mtx.sub.sub, features = c("EGR1", "IER3", "JUND", "LIMD2", "TNFAIP3", "TRBC1"),
          group.by = "dataset", 
          slot = "data",group.colors = c("#A4D38E", "#5D9BBE", "#E45A5F")) +
  ggplot2::scale_fill_gradient(low = "white", high = "#9E0142")
dev.off()

cell.num <- rbind(table(mtx@meta.data[["dataset"]]),
                  table(mtx.sub.sub@meta.data[["dataset"]]))
rownames(cell.num) <- c("PreFilter", "AfterFilter")
cell.num <- as.data.frame(t(cell.num))
cell.num$pct <- cell.num[,2]/cell.num[,1]*100

library(ggplot2)
pdf("../results/pct.pdf", width = 4, height = 6, useDingbats = F)
ggplot(cell.num, aes(x = rownames(cell.num), y = pct, fill = rownames(cell.num)))+
  geom_bar(stat = "identity", width = 0.6)+
  theme(  axis.text.x = element_text(angle = 0, size = 10, color="black"),
          axis.text.y = element_text(angle = 0, size = 10, color="black"),
          axis.title = element_text(size = 10, color = "black"),
          axis.ticks.length = unit(.2, "cm"),
          axis.ticks = element_line(colour = "black"),
          panel.grid.major.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position = "none",
          legend.text = element_text(size = 10),
          strip.background = element_blank(),
          
          #panel.border = element_rect(colour = "black", fill = NA, size = .5)
  )+
  scale_color_manual(values=c("#A4D38E", "#5D9BBE", "#E45A5F"))
dev.off()
