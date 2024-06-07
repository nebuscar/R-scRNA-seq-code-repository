library(Seurat)
library(dplyr)
devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
library(ggplot2)
library(showtext)

##########
sc.combined <- readRDS("./tmp/competition/sc.combined.after_umap.cca~without QC.rds")

##########
# SCE
counts <- GetAssayData(sc.combined, assay = "RNA", layer = "counts")
cell_metadata <- sc.combined@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(counts))
rownames(gene_annotation) <- rownames(counts)
cds <- new_cell_data_set(counts, cell_metadata = cell_metadata, gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 10)

cds <- reduce_dimension(cds)
SingleCellExperiment::reducedDims(cds)[["UMAP"]] <- sc.combined@reductions$umap.cca@cell.embeddings
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
           graph_label_size=8,)
showtext_auto()
pseudotime_plot <- pseudotime_plot &
  theme(
    axis.title.x = element_text(size = 20, family = "Arial"),
    axis.title.y = element_text(size = 20, family = "Arial"),
    axis.text.x = element_text(size = 20, family = "Arial"),
    axis.text.y = element_text(size = 20, family = "Arial")
  )

pseudotime_plot <- pseudotime_plot + theme(text = element_text(size = 20))
pdf(file = "./results/competition/figure2/pseudotime_plot_umap.cca.pdf", width = 10, height = 5)
print(pseudotime_plot)
dev.off()
showtext_auto(FALSE)
