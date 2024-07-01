devtools::install_github("jinworks/CellChat")

library(dplyr)
library(tidyr)
library(CellChat)
library(Seurat)

################################################################################
# Load dataset
################################################################################
sc.raw <- readRDS("./tmp/sc.common_glioma_counts.matrix.rds")
sc.raw.meta <- readRDS("tmp/sc.meta_sub_sub_annotation+scissor.rds")
###################################### optical: randomly sample
# # random sub-sample
# set.seed(42)
# random <- sample(1:ncol(sc.raw), 60000, replace = FALSE)
# sc.raw.subset <- sc.raw[, random]
# common_sc.raw.id <- intersect(colnames(sc.raw.subset), sc.raw.meta$NAME)
# sc.meta.subset <- sc.raw.meta[sc.raw.meta$NAME %in% common_sc.raw.id, ]

# CellChat
cellchat <- createCellChat(sc.raw, 
                           sc.raw.meta, 
                           group.by = "annotation_cluster",
                           # group.by = "scissor_tcells",
                           )
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#future::plan("multisession", workers = 10)                                    # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) # num_highly variable lr pars used for signaling inference:946
#options(future.globals.maxSize = 400000000000000000000000000000000)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# saveRDS(cellchat, "tmp/cellchat_scissor_tcells.rds")
saveRDS(cellchat, "tmp/cellchat_annotation.rds")
#rm(data.input, data.raw, meta, sc.raw.meta, sc.raw)