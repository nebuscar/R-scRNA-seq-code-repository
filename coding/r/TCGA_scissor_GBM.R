BiocManager::install("Scissor")
BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")

sc_dataset <- readRDS("data/seurat/glioma.rm.batch.rds")
library(Matrix)
sc_dataset <- as(tmp, "sparseMatrix") 
sc_dataset <- Scissor::Seurat_preprocessing(sc_dataset, verbose = F)


################################################################################
##################data.bulk##################
query <- TCGAbiolinks:: GDCquery(project = "TCGA-GBM",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")
tmp <- query[[1]][[1]]
tmp2 <- tmp[substr(tmp$sample.submitter_id, 14, 15) < 10, ]
tmp2 <- tmp2[order(tmp2$cases.submitter_id, tmp2$sample.submitter_id, decreasing = F),]
tmp2 <- tmp2[!duplicated(tmp2$cases.submitter_id),]
query[[1]][[1]] <- tmp2
TCGAbiolinks::GDCdownload(query)
exp <- TCGAbiolinks::GDCprepare(query = query)
exp_tpm <- SummarizedExperiment::assay(exp, "tpm_unstrand")
write.table(exp_tpm, "data/raw/TCGA_GBM.tsv", sep = "\t", row.names = T, col.names = T, quote = F)


# covert gene name
gtf <- readRDS("data/gtf.rds")
exp_tpm <- read.table("data/raw/TCGA_GBM.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
exp_median <- data.frame(median = apply(exp_tpm, 1, median))
exp_median$gene_id <- substr(rownames(exp_median), 1, 15)
gtf_match <- gtf[match(exp_median$gene_id, gtf$gene_id), ]
exp_median <- cbind(exp_median, gtf_match$gene_name)
colnames(exp_median)[3] <- "gene_name"

exp_median <- exp_median[!is.na(exp_median$gene_name),]
#exp_median <- exp_median[!is.na(exp_median$median),]


table(duplicated(exp_median$gene_name))
exp_median <- exp_median[order(exp_median$median, decreasing = TRUE), ]
exp_median <- exp_median[!duplicated(exp_median$gene_name), ]
exp_tpm <- exp_tpm[match(rownames(exp_median), rownames(exp_tpm)), ]
rownames(exp_tpm) <- exp_median$gene_name
colnames(exp_tpm) <- substr(colnames(exp_tpm), 1, 12)

# colnames convert "." to "-"
colnames(exp_tpm) <- gsub("\\.", "-", colnames(exp_tpm))
write.table(exp_tpm, "data/raw/TCGA_GBM_TPM_geneName.tsv", sep = "\t", row.names = T, col.names = T, quote = F)

##################data.clinical##################
clinic <- TCGAbiolinks::GDCquery_clinic(project = "TCGA-GBM", type = "clinical")
clinical <- clinic[, c("bcr_patient_barcode", 
                       "vital_status", 
                       "days_to_death", 
                       "days_to_last_follow_up",
                       "race",
                       "age_at_diagnosis",
                       "gender",
                       "ajcc_pathologic_stage")]
# rownames(clinical) <- NULL
rownames(clinical) <- clinical$bcr_patient_barcode
clinical[,3][is.na(clinical[,3])] = 0
clinical[,4][is.na(clinical[,4])] = 0
clinical$time <- (as.numeric(clinical[,3]) + as.numeric(clinical[,4]))
clinical$event <- ifelse(clinical$vital_status == "Alive", 0, 1)

saveRDS(clinical,file = "data/raw_TCGA/clinical_GBM.rds")

##################合并match bulk与clinicala##################
bulk_survival <- readRDS("data/raw_TCGA/clinical.rds") # n=617
bulk_dataset <- read.table("data/raw_TCGA/TCGA_GBM_TPM_geneName.tsv") # n=163
# colnames convert "." to "-"
colnames(bulk_dataset) <- gsub("\\.", "-", colnames(bulk_dataset))

match_order <- match(colnames(bulk_dataset),rownames(bulk_survival))
##tmp1 <- bulk_survival[match_order,]
bulk_survival <- bulk_survival[match_order,]
##tmp2 <- bulk_dataset[,match_order]
common <- intersect(rownames(bulk_survival),colnames(bulk_dataset))
bulk_dataset <- bulk_dataset[,common]
bulk_survival <- bulk_survival[common,]

all(colnames(bulk_dataset) == bulk_survival$bcr_patient_barcode)
phenotype <- bulk_survival[,c("time","event")]
colnames(phenotype) <- c("time","status")
saveRDS(phenotype,file = "data/phenotype.rds")
tmp <- as.matrix(bulk_dataset)
saveRDS(tmp,file = "data/bulk_dataset.rds")

##################data.sc##################
load("data/tcells_sct.rdata")
tcells_sub <- tcells
GBM_sc_dataset <- GetAssayData(object = tcells_sub, assay = "RNA", slot = "counts") # get count matrix data
sc_dataset <- Scissor::Seurat_preprocessing(GBM_sc_dataset, verbose = F)
saveRDS(sc_dataset,file = "data/sc_dataset.rds")

##################Scissor##################
bulk_dataset <- readRDS("data/downstream/scissor_data/bulk_dataset.rds")
sc_dataset <- readRDS("data/downstream/scissor_data/sc_dataset.rds")
phenotype <- readRDS("data/downstream/scissor_data/phenotype.rds")

info <- Scissor::Scissor(bulk_dataset, 
                         sc_dataset, 
                         phenotype, 
                         alpha = 0.05,
                         family = "cox")
Scissor_select <- rep("Background cells",ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos1$Scissor_pos] <- "Scissor+ cell"
Scissor_select[infos1$Scissor_neg] <- "Scissor- cell"

# add meta.scissor
sc_dataset <- AddMetaData(sc_dataset,metadata = Scissor_select,col.name = "scissor")
Dim_plot1 <- DimPlot(sc_dataset,reduction = 'umap',group.by = 'scissor',cols = c('grey','royalblue','indianred1'),pt.size = 2)
Dim_plot1
saveRDS(sc_dataset,file = "data/sc_scissor.rds")

