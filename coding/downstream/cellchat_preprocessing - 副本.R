devtools::install_github("jinworks/CellChat")

library(dplyr)
library(tidyr)
library(CellChat)
library(Seurat)


# load raw data.tcells
sc.tcells <- readRDS("data/downstream/markers/tcells_clutser_markers.rds") #num_tcells:19434
sc.raw.meta <- (read.csv(gzfile("data/download/sc_dataset/Meta_GBM.txt.gz"),# metadata
                           sep = ",",
                           header = T))[-1,] #不可用read.csv
sc.raw.meta_sub <- sc.raw.meta[sc.raw.meta$Assignment %in% c("TCells", "Glioma", "Myeloid"),] # meda.interest.CellType 
sc.raw.meta.tcells <- sc.raw.meta[sc.raw.meta$Assignment %in% "TCells",]  # num_tcells:19570
sc.raw.meta.Glioma <- sc.raw.meta[sc.raw.meta$Assignment %in% "Glioma",]  # num_tcells:81764
sc.raw.meta.Myeloid <- sc.raw.meta[sc.raw.meta$Assignment %in% "Myeloid",]  # num_tcells:90954
sc.raw.glioma <- readRDS("data/seurat/glioma_v5.rds")


# intersect tcells & meta | rename TCells Type |annotation_cluster
common_annotation.tcells.id <- intersect(sc.raw.meta_sub$NAME, colnames(sc.tcells)) # num_tcells:19434
common_annotation.tcells <- sc.tcells[, colnames(sc.tcells) %in% common_annotation.tcells.id]
common_meta.tcells <- sc.raw.meta_sub[sc.raw.meta_sub$NAME %in% common_annotation.tcells.id, ]
common_annotation_cluster.tcells <- as.matrix(Idents(common_annotation.tcells))
sc.raw.meta$annotation_cluster <- sc.raw.meta$Assignment
sc.raw.meta$annotation_cluster[sc.raw.meta$annotation_cluster == "TCells"] <- common_annotation_cluster.tcells
saveRDS(sc.raw.meta, "tmp/sc.meta_sub_annotation.rds")

# rename TCells Type |scissor_cells
sc.raw.meta_sub <- readRDS("tmp/sc.meta_sub_annotation.rds")
sc.tcells <- readRDS("data/downstream/markers/tcells_clutser_markers.rds") #num_tcells:19434

common_scissor.tcells.id <- intersect(sc.raw.meta_sub$NAME, colnames(sc.tcells)) # num_scissor_cells:19434
common_scissor.tcells <- sc.tcells[ , colnames(sc.tcells) %in% common_scissor.tcells.id]
common_scissor_meta.tcells <- sc.raw.meta_sub[sc.raw.meta_sub$NAME %in% common_scissor.tcells.id, ]
common_scissor_tcells <- as.matrix(common_scissor.tcells$scissor)
sc.raw.meta_sub$scissor_tcells <- sc.raw.meta_sub$Assignment
sc.raw.meta_sub$scissor_tcells[sc.raw.meta_sub$scissor_tcells == "TCells"] <- common_scissor_tcells
saveRDS(sc.raw.meta_sub, "tmp/sc.meta_sub_annotation+scissor.rds")


# intersect glioma & meta
sc.raw.meta <- readRDS("tmp/sc.meta_sub_annotation+scissor.rds")
sc.raw.glioma <- readRDS("data/seurat/glioma_v5.rds")
common_sc.glioma.id <- intersect(sc.raw.meta$NAME, colnames(sc.raw.glioma))
common_sc.glioma <- sc.raw.glioma[, common_sc.glioma.id]
common_sc.glioma <- GetAssayData(common_sc.glioma, assay = "RNA", layer = "counts")
common_sc.glioma <- normalizeData(common_sc.glioma)
saveRDS(common_sc.glioma, "tmp/sc.common_glioma_counts.matrix.rds")


################################################################################
# Load dataset
################################################################################
sc.raw <- readRDS("tmp/sc.common_glioma_counts.matrix.rds")
sc.raw.meta <- readRDS("tmp/sc.meta_sub_annotation+scissor.rds")
###################################### optical: randomly sample
# random sub-sample
set.seed(100)
random <- sample(1:ncol(sc.raw), 60000, replace = FALSE)
sc.raw.subset <- sc.raw[, random]
common_sc.raw.id <- intersect(colnames(sc.raw.subset), sc.raw.meta$NAME)
sc.meta.subset <- sc.raw.meta[sc.raw.meta$NAME %in% common_sc.raw.id, ]

# CellChat
cellchat <- createCellChat(sc.raw.subset, 
                           sc.meta.subset, 
                           # group.by = "annotation_cluster",
                           group.by = "scissor_tcells",
                           )
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#future::plan("multisession", workers = 10)                                    # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) # num_highly variable lr pars used for signaling inference:946
#options(future.globals.maxSize = 400000000000000000000000000000000)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, "tmp/cellchat_scissor_tcells.rds")
# saveRDS(cellchat, "tmp/cellchat_annotation.rds")
#rm(data.input, data.raw, meta, sc.raw.meta, sc.raw)