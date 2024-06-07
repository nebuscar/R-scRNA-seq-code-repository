devtools::install_github("jinworks/CellChat")

library(dplyr)
library(tidyr)
library(CellChat)
setwd("F:/Projects/Forstudents/NASH")


sc.raw <- read.table(gzfile("raw/GSE166504_cell_raw_counts.20220204.txt.gz"),   # count matrix 
                    header = TRUE)
sc.raw <- as.matrix(sc.raw)
sc.raw <- as(sc.raw, "dgCMatrix")
saveRDS(sc.raw, "tmp/sc_raw.RDS")

################################################################################
# Load dataset
################################################################################
sc.raw <- readRDS("tmp/sc_raw.RDS")
sc.raw.meta <- read.table(gzfile("raw/GSE166504_cell_metadata.20220204.tsv.gz"),# metadata 
                          sep = "\t", header = T)
sc.raw.meta$CellID <- paste(sc.raw.meta$FileName, sc.raw.meta$CellID, sep = "_")
sc.raw.meta <- sc.raw.meta %>%
  separate(FileName, into = c("Group", "Time", "Donor", "Capture"), sep = "_")

# CellChat
CellChatDB.use <- CellChatDB.mouse

cellchat.list <- list()
for (dataset in unique(sort(sc.raw.meta$Time))){
  cell.use <- sc.raw.meta$CellID[sc.raw.meta$Time == dataset]
  data.input <- sc.raw[, cell.use]
  meta <- sc.raw.meta[sc.raw.meta$CellID %in% cell.use,]
  cellchat.list[[dataset]] <- createCellChat(object = data.input, 
                                             meta = meta, 
                                             group.by = "CellType")
  cellchat.list[[dataset]]@DB <- CellChatDB.use                                                     
  cellchat.list[[dataset]] <- subsetData(cellchat.list[[dataset]])              # This step is necessary even if using the whole database
  future::plan("multisession", workers = 10)                                    # do parallel
  cellchat.list[[dataset]] <- identifyOverExpressedGenes(cellchat.list[[dataset]])
  cellchat.list[[dataset]] <- identifyOverExpressedInteractions(cellchat.list[[dataset]])
  cellchat.list[[dataset]] <- computeCommunProb(cellchat.list[[dataset]], type = "triMean")
  cellchat.list[[dataset]] <- filterCommunication(cellchat.list[[dataset]], min.cells = 10)
  cellchat.list[[dataset]] <- computeCommunProbPathway(cellchat.list[[dataset]])
  cellchat.list[[dataset]] <- aggregateNet(cellchat.list[[dataset]])
}
saveRDS(cellchat.list, "tmp/cellchat.list.rds")
#rm(data.input, data.raw, meta, sc.raw.meta, sc.raw)