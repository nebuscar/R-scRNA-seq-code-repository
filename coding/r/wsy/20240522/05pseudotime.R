library(Seurat)
library(dplyr)
library(monocle3)
library(ggplot2)

##########
sc.combined <- readRDS("./data/sc.sub.rds")

##########
# SCE
counts <- GetAssayData(sc.combined, assay = "RNA", layer = "counts")
cell_metadata <- sc.combined@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(counts))
rownames(gene_annotation) <- rownames(counts)
cds <- new_cell_data_set(counts, cell_metadata = cell_metadata, gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 10)

cds <- reduce_dimension(cds)
SingleCellExperiment::reducedDims(cds)[["UMAP"]] <- sc.combined@reductions$umap@cell.embeddings
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

get_earliest_principal_node <- function(cds, time_bin){
  #  cell_ids <- which(colData(cds)[, "orig.ident"] == time_bin)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex))))] 
  root_pr_nodes
}

cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds))

pseudotime_plot <- plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=5,)
ggsave(filename = "./results/wsy/pseudotime_plot.png", plot = pseudotime_plot, width = 10, height = 6)
