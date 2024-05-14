
library(dplyr)
library(monocle3)
library(Seurat)

# SCE
tmp <- readRDS("data/downstream/markers/tcells_clutser_markers.rds")
DimPlot(tmp, split.by = "scissor")
data <- GetAssayData(tmp, assay = "RNA", layer = "counts")
cell_metadata <- tmp@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 10)

cds <- reduce_dimension(cds)
SingleCellExperiment::reducedDims(cds)[["UMAP"]] <- tmp@reductions[["umap"]]@cell.embeddings
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

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=5,)
