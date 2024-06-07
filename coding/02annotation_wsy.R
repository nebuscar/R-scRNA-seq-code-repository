workdir = "F:/Projects/Skull_scRNAseq"  # where you place the data and results
setwd(workdir)
library(dplyr)
library(Seurat)

sk.scData <- readRDS("./tmp/sk.scData.rds")

# find markers, report only the positive
sk.scData.markers <- FindAllMarkers(sk.scData, only.pos = TRUE)
sk.scData.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)
sk.scData.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
write.csv(sk.scData.markers, "./result/res01_findallmarkers.csv", 
          row.names = TRUE, quote = FALSE)
write.csv(top10, "./result/res01_top10gene.csv", 
          row.names = TRUE, quote = FALSE)

pdf("./result/heatmap_top10.pdf", width = 20, height = 10, useDingbats = F)
DoHeatmap(sk.scData, features = top10$gene) + NoLegend()
dev.off()

# Annotated by ChatGPT3.5
new.cluster.ids <- c("T lymphocytes", "Eosinophils_1", "Neutrophils_1", 
                     "Neutrophils_2", "Macrophages", "Erythrocytes_1",
                     "B lymphocytes", "Erythrocytes_2", "Neutrophils_3",
                     "Plasma B cells", "pDC")
names(new.cluster.ids) <- levels(sk.scData)
sk.scData <- RenameIdents(sk.scData, new.cluster.ids)
pdf("./result/dimplot_integrated_gpt.pdf", width = 6, height = 6, useDingbats = F)
DimPlot(sk.scData, reduction = "umap", raster = FALSE, label = TRUE, label.size = 4)+ NoLegend()
dev.off()

pdf("./result/dimplot_integrated_donor.pdf", width = 18, height = 6, useDingbats = F)
DimPlot(sk.scData, reduction = "umap", split.by = "orig.ident", 
        label = F, label.size = 4)+ NoLegend()
dev.off()
################################################################################
# TCR
#devtools::install_github("ncborcherding/scRepertoire")

file.list <- list.files("rawData/", pattern = "matrix$")

sk.tcr <- list()
sk.sample <- c()
for (file in file.list){
  file.t <- list.files(paste("rawData/", file, sep = ""), pattern = "csv")
  sk.tcr[[file]] <- read.csv(paste("rawData", file, file.t, sep = "/"))
  sk.sample <- c(sk.sample, strsplit(file, "-")[[1]][1])
  rm(file.t)
}
combined.TCR <- scRepertoire::combineTCR(sk.tcr, 
                           samples = sk.sample,
                           removeNA = FALSE, 
                           removeMulti = FALSE, 
                           filterMulti = FALSE)

sce <- Seurat::as.SingleCellExperiment(scRep_example)

# QC
sk.tcr <- do.call(rbind.data.frame, combined.TCR)
#sk.tcr <- na.omit(sk.tcr)
sk.scData@meta.data$isTcell <- 0
sk.scData.tmpMeta <- sk.scData@meta.data %>%
  mutate(isTcell = case_when(colnames(sk.scData) %in% sk.tcr$barcode ~ "1",
                             .default = as.character(0)))
sk.scData@meta.data$isTcell <- sk.scData.tmpMeta$isTcell
rm(sk.scData.tmpMeta)

pdf("./result/dimplot_integrated_tcr.pdf", width = 24, height = 6, useDingbats = F)
DimPlot(sk.scData, reduction = "umap", group.by = c("isTcell"), split.by = "orig.ident",
        cols = c("grey", "red"), pt.size = .01, order = TRUE)+ NoLegend()
dev.off()


sk.08 <- subset(sk.scData, subset = RNA_snn_res.0.1 == "8")
sk.tcr.08 <- sk.tcr[sk.tcr$barcode %in% colnames(sk.08), ]

