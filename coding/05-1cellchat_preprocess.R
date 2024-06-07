library(Seurat)


# load raw data.tcells
sc.tcells <- readRDS("data/downstream/markers/tcells_clutser_markers.rds") #num_tcells:19434
sc.raw.meta <- (read.csv(gzfile("data/download/sc_dataset/Meta_GBM.txt.gz"),# metadata
                         sep = ",",
                         header = T))[-1,] #不可用read.csv
sc.raw.meta_sub <- sc.raw.meta[sc.raw.meta$Assignment %in% c("TCells", "Glioma", "Myeloid"),] # meta.interest.CellType 
sc.raw.meta.tcells <- sc.raw.meta[sc.raw.meta$Assignment %in% "TCells",]  # num_tcells:19570
sc.raw.meta.Glioma <- sc.raw.meta[sc.raw.meta$Assignment %in% "Glioma",]  # num_tcells:81764
sc.raw.meta.Myeloid <- sc.raw.meta[sc.raw.meta$Assignment %in% "Myeloid",]  # num_tcells:90954

# random subsut 20.pct of sc.raw.meta.Glioma$Myeloid
# num_Glioma_sub:16353 | num_Myeloid_sub:18191
set.seed(42)
datasets <- list(Glioma = sc.raw.meta.Glioma,
                 Myeloid = sc.raw.meta.Myeloid)
for (dataset_name in names(datasets)) {
  num_samples <- round(0.2 * nrow(datasets[[dataset_name]]))
  random_indices <- sample(1:nrow(datasets[[dataset_name]]), num_samples, replace = FALSE)
  assign(paste("sc.raw.meta.", dataset_name, "_subset", sep = ""), datasets[[dataset_name]][random_indices,])
}
sc.raw.meta_sub_sub <- rbind(sc.raw.meta.Glioma_subset, sc.raw.meta.Myeloid_subset, sc.raw.meta.tcells)

rm(sc.raw.meta, datasets, 
   sc.raw.meta_sub, 
   sc.raw.meta.Glioma, sc.raw.meta.Myeloid)
sc.raw.glioma <- readRDS("data/seurat/glioma_v5.rds")


# intersect tcells & meta | rename TCells Type |annotation_cluster
common_annotation.tcells.id <- intersect(sc.raw.meta_sub_sub$NAME, colnames(sc.tcells)) # num_tcells:19434
common_annotation.tcells <- sc.tcells[, colnames(sc.tcells) %in% common_annotation.tcells.id]
common_meta.tcells <- sc.raw.meta_sub_sub[sc.raw.meta_sub_sub$NAME %in% common_annotation.tcells.id, ]
common_annotation_cluster.tcells <- as.matrix(Idents(common_annotation.tcells))
sc.raw.meta_sub_sub$annotation_cluster <- sc.raw.meta_sub_sub$Assignment
sc.raw.meta_sub_sub$annotation_cluster[sc.raw.meta_sub_sub$annotation_cluster == "TCells"] <- common_annotation_cluster.tcells
saveRDS(sc.raw.meta_sub_sub, "tmp/sc.meta_sub_sub_annotation.rds")

# rename TCells Type |scissor_cells
sc.raw.meta_sub_sub <- readRDS("tmp/sc.meta_sub_sub_annotation.rds")
sc.tcells <- readRDS("data/downstream/markers/tcells_clutser_markers.rds") #num_tcells:19434

common_scissor.tcells.id <- intersect(sc.raw.meta_sub_sub$NAME, colnames(sc.tcells)) # num_scissor_cells:19434
common_scissor.tcells <- sc.tcells[ , colnames(sc.tcells) %in% common_scissor.tcells.id]
common_scissor_meta.tcells <- sc.raw.meta_sub_sub[sc.raw.meta_sub_sub$NAME %in% common_scissor.tcells.id, ]
common_scissor_tcells <- as.matrix(common_scissor.tcells$scissor)
sc.raw.meta_sub_sub$scissor_tcells <- sc.raw.meta_sub_sub$Assignment
sc.raw.meta_sub_sub$scissor_tcells[sc.raw.meta_sub_sub$scissor_tcells == "TCells"] <- common_scissor_tcells
saveRDS(sc.raw.meta_sub_sub, "tmp/sc.meta_sub_sub_annotation+scissor.rds")


# intersect glioma & meta
sc.raw.meta <- readRDS("tmp/sc.meta_sub_sub_annotation+scissor.rds")
sc.raw.glioma <- readRDS("data/seurat/glioma_v5.rds")
common_sc.glioma.id <- intersect(sc.raw.meta$NAME, colnames(sc.raw.glioma))
common_sc.glioma <- sc.raw.glioma[, common_sc.glioma.id]
common_sc.glioma <- GetAssayData(common_sc.glioma, assay = "RNA", layer = "counts")
common_sc.glioma <- normalizeData(common_sc.glioma)
saveRDS(common_sc.glioma, "tmp/sc.common_glioma_counts.matrix.rds")
