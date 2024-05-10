##################CGGA-693##################
library(Seurat)
library(Matrix)
data.bulk <- read.table("data/download/sc.data.CGGA/mRNAseq_693/bulkdata/CGGA.mRNAseq_693.RSEM-genes.20200506.txt", 
                        header = T, 
                        sep = "",
)
data.bulk <- `rownames<-`(data.bulk[,-1], data.bulk$Gene_Name)
data.bulk <- as.matrix(data.bulk)

data.clinical <- read.table("data/download/sc.data.CGGA/mRNAseq_693/clinical/CGGA.mRNAseq_693_clinical.20200506.txt",
                            header = T,
                            sep = "\t",
)
data.clinical.gbm <- data.clinical[grepl("GBM$", data.clinical$Histology), ]
data.clinical.gbm <- `rownames<-`(data.clinical.gbm[,-1], data.clinical.gbm$CGGA_ID)
data.clinical.gbm <- data.clinical.gbm[,c("OS","Censor..alive.0..dead.1.")]
colnames(data.clinical.gbm) <- c("time","status")
data.clinical.gbm <- na.omit(data.clinical.gbm)

common <- intersect(colnames(data.bulk), rownames(data.clinical.gbm))
data.bulk <- data.bulk[,common]
data.clinical.gbm <- data.clinical.gbm[common,]

# data.sc <- readRDS("data/seurat/tcells_v5.rds")
# data.sc <- GetAssayData(object = data.sc, assay = "RNA", slot = "counts") # get count matrix data
# data.sc <- Scissor::Seurat_preprocessing(data.sc, verbose = F)
# # 重命名Assay，切换v3版本
# data.sc <- RenameAssays(object = data.sc, assay.name = "RNA", new.assay.name = "RNA5")
# data.sc[["RNA"]] <- as(object = data.sc[["RNA5"]], Class = "Assay")
# DefaultAssay(data.sc) <- "RNA"
# ##sc_dataset@assays$RNA <- NULL

data.sc <- readRDS("data/downstream/scissor_data/CGGA.693/data.sc.rds")
# Scissor pct.scissor:2.213% | neg cells:268 | pos cells:165
info.693 <- Scissor::Scissor(data.bulk, data.sc, data.clinical.gbm, alpha = 0.05,
                             family = "cox")
saveRDS(data.bulk, file = "data/downstream/scissor_data/CGGA.693/data.bulk.rds")
saveRDS(data.clinical.gbm, file = "data/downstream/scissor_data/CGGA.693/data.clinical.rds")
saveRDS(data.sc, file = "data/downstream/scissor_data/CGGA.693/data.sc.rds")
saveRDS(info.693, file = "data/downstream/scissor_data/CGGA.693/scissor.info.693.rds")

data.bulk <- readRDS("data/downstream/scissor_data/CGGA.693/data.bulk.rds")
data.clinical.gbm <- readRDS("data/downstream/scissor_data/CGGA.693/data.clinical.rds")
data.sc <- readRDS("data/downstream/scissor_data/CGGA.693/data.sc.rds")
info.693 <- readRDS("data/downstream/scissor_data/CGGA.693/scissor.info.693.rds")

# add metadata.scissor
Scissor_select <- rep("Background cells",ncol(data.sc))
names(Scissor_select) <- colnames(data.sc)
Scissor_select[info.693$Scissor_pos] <- "Scissor+ cell"
Scissor_select[info.693$Scissor_neg] <- "Scissor- cell"
data.sc <- AddMetaData(data.sc,metadata = Scissor_select,col.name = "scissor")
DimPlot(data.sc,reduction = 'umap',group.by = 'scissor',cols = c('grey','royalblue','indianred1'),pt.size = 2, order =T)
saveRDS(data.sc,file = "data/downstream/scissor_data/CGGA.693/data.sc.scissor.rds")

# match data.sc.tcells和data.sc.scissor
data.sc.tcells <- readRDS("data/downstream/markers/tcells_clutser_markers.rds")
data.sc.scissor <- readRDS("data/downstream/scissor_data/CGGA.693/data.sc.scissor.rds")
scissor_sub <- data.sc.scissor@meta.data
scissor_sub <- scissor_sub[match(colnames(data.sc.tcells), rownames(scissor_sub)),]
data.sc.tcells@meta.data$scissor <- scissor_sub$scissor

saveRDS(data.sc.tcells, file = "data/downstream/scissor_data/CGGA.693/data.sc.scissor.tcells")

##################pct.693.scissor.neg##################
library(Seurat)
data.sc <- readRDS("data/downstream/scissor_data/CGGA.693/data.sc.scissor.tcells")
sc.scissor.neg <- subset(data.sc, subset = scissor == "Scissor- cell")
# 计算所有细胞中scissor- cell占比：4.81 | 5.52 |1.37
print(paste("所有细胞中Scissor- cell 占比：", 
            total.scissor.neg.percent <- ncol(sc.scissor.neg) / ncol(data.sc) * 100, "%"))

# 计算naive中scissor- cell占比:18.89 |15.10
sc.naive <- subset(data.sc, subset = seurat_clusters == "3")
sc.naive.scissor.neg <- subset(sc.naive, subset = scissor == "Scissor- cell")
print(paste("naive中Scissor- cell 占比：", 
            naive.scissor.neg.percent <- ncol(sc.naive.scissor.neg) / ncol(sc.naive) * 100, "%"))

# 统计不同样本中scissor- cell占比
data.sc <- readRDS("data/downstream/scissor_data/CGGA.693/data.sc.scissor.tcells")
metadata <- read.csv("data/download/sc_dataset/Meta_GBM.txt.gz")
metadata <- metadata[metadata$Assignment %in% "TCells",] #meta.tcells
metadata <- metadata[metadata$NAME %in% rownames(data.sc@meta.data), ] # meta.matched
#table(metadata$donor_id)
# sample.id <- metadata[metadata$donor_id %in% "human_ndGBM-11",][, 1]
# sample <- subset(data.sc, barcode %in% sample.id)
# sample.scissor.neg <- subset(sample, subset = scissor == "Scissor- cell")
# print(paste("样本中Scissor- cell 占比：", 
#             total.scissor.neg.percent <- ncol(sample.scissor.neg) / ncol(sample) * 100, "%"))

data.sc@meta.data$patient <- metadata$Patient
cells.693 <- data.sc@meta.data[data.sc@meta.data$seurat_clusters == 1, c( "scissor", "patient")]
pct.scissor <- table(data.sc@meta.data[data.sc@meta.data$seurat_clusters == 1, c( "scissor", "patient")])
pct.scissor <- rbind(pct.scissor, colSums(pct.scissor))
pct.scissor <- t(pct.scissor)
pct.scissor <- as.data.frame(pct.scissor)
pct.scissor$pct <- pct.scissor$`Scissor- cell`/pct.scissor$V4
write.csv(pct.scissor, "results/pct.693.tcells.scissor.neg.csv", row.names = F)
saveRDS(cells.693, file = "data/downstream/scissor_data/CGGA.693/scissor.cells.693.rds")
