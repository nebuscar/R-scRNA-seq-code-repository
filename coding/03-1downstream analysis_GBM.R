setwd("E:/GitHub/scRNA-seq_GBM")

library(Seurat)
library(dplyr)

# load data
sc_dataset <- readRDS("data/seurat/tcells_rm_batch-1.rds")
tmp <- sc_dataset

# re-join layers after integration
tmp[["RNA"]] <- JoinLayers(tmp[["RNA"]])
tmp <- FindNeighbors(tmp, reduction = "integrated.rpca", dims = 1:30)
tmp <- FindClusters(tmp, resolution = 0.8)
tmp <- RunUMAP(tmp, dims = 1:30, reduction = "integrated.rpca")
DimPlot(tmp, reduction = "umap", split.by = "scissor", label.size = 6, label = T, pt.size = 1)
saveRDS(tmp, file = "data/downstream/tcells_resolution=0.8.rds")

# find markers
tmp <- readRDS("data/downstream/tcells_resolution=0.8.rds")
tcells.markers <- FindAllMarkers(tmp, only.pos = T)
tcells.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

#DoHeatmap(tmp, features = top10$gene) + NoLegend()
saveRDS(tcells.markers, file = "data/downstream/markers/tcells.markers.rds")
saveRDS(top10, file = "data/downstream/markers/top10.rds")

##################细胞注释##################
tmp <- readRDS("data/downstream/markers/tcells_resolution=0.8.rds")
tcells.markers <- readRDS("data/downstream/markers/tcells.markers.rds")
top10 <- readRDS("data/downstream/markers/top10.rds")
#genes_cluster_0 <- top10$gene[top10$cluster == "0"]
DoHeatmap(tmp, features = top10$gene)
data <- readxl::read_excel("data/downstream/markers/golden_key_markers_total.csv")


DefaultAssay(tmp)
data.sc <- data$`Feature set`
data.sc <- data.sc[data.sc != "ND"]
data.sc <- unique(sort(data.sc))
DoHeatmap(tmp, 
          #data.sc,
          #slot = "counts",
          features = c("NKG7", "CD8A", "CD4", "GNLY", "KLRB1", "IL17A"))

DotPlot(tmp, features = c("CXCR3", "CCR6", "CD4", "RORC"))
FeaturePlot(tmp, features = c("CXCR6", "RORA"))

new.cluster.ids <- c("0.CD8", "1.CD4", 
                     
                     "2.CD8", "3.Naive", 
                     
                     "4.Tregs", "5.CD8", 
                     
                     "6.Innate", "7.CD4 CTL", 
                     
                     "8.CD8", "9.γδ T", 
                     
                     "10.Innate", "11.NKT",
                     
                     "12.unknown", "13.CD4", 
                     
                     "14.CD8", "15.CD4")
names(new.cluster.ids) <- levels(tmp)
tmp <- RenameIdents(tmp, new.cluster.ids)
saveRDS(tmp, file = "data/downstream/markers/tcells_clutser_markers.rds")
DimPlot(tmp, reduction = "umap", split.by = "scissor", label = T, pt.size = 1, label.size = 6)

##################DEseq##################
tmp <- readRDS("data/downstream/markers/tcells_clutser_markers.rds")

# Find DE features
cluster1.de.markers <- FindMarkers(tmp, ident.1 = "1.Naive", ident.2 = NULL, only.pos = T, min.pct = 0.5)
cluster10.de.markers <- FindMarkers(tmp, ident.1 = "10.Innate", ident.2 = NULL, only.pos = T, min.pct = 0.25)
Cluster3_10.de.markers <- FindMarkers(tmp, ident.1 = "3.CD4", ident.2 = "1.Innate", only.pos = T, min.pct = 0.25)

library(ggplot2)
library(ggrepel)
# 读取聚类标记数据
cluster_markers <- readRDS("data/downstream/markers/tcells_clutser_markers.rds")
cluster1.de.markers <- FindMarkers(cluster_markers, ident.1 = "3.Naive", ident.2 = NULL, only.pos = F, min.pct = 0.25)

#自定义阈值
log2FC = 1
padj = 0.01 

cluster1.de.markers$threshold="ns";
cluster1.de.markers[which(cluster1.de.markers$avg_log2FC  > log2FC & cluster1.de.markers$p_val_adj <padj),]$threshold="up";
cluster1.de.markers[which(cluster1.de.markers$avg_log2FC  < (-log2FC) & cluster1.de.markers$p_val_adj < padj),]$threshold="down";
cluster1.de.markers$threshold=factor(cluster1.de.markers$threshold, levels=c('down','ns','up'))

p1 <- ggplot(data=cluster1.de.markers, aes(x=avg_log2FC, y=-log10(p_val_adj), color=threshold)) +
  geom_point(alpha=0.8, size=0.8) +
  geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color="grey")+
  geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+
  #labs(title= ifelse(""==title, "", paste("DEG:", title)))+
  xlab(bquote(Log[2]*FoldChange))+
  ylab(bquote(-Log[10]*italic(P.adj)) )+
  theme_classic(base_size = 14) +
  scale_color_manual('',labels=c(paste0("down(",table(cluster1.de.markers$threshold)[[1]],')'),'ns',
                                 paste0("up(",table(cluster1.de.markers$threshold)[[3]],')' )),
                     values=c("blue", "grey","red" ) )+
  guides(color=guide_legend(override.aes = list(size=3, alpha=1)))

##添加基因名标记###
data = subset(cluster1.de.markers, cluster1.de.markers$p_val_adj < padj & abs(cluster1.de.markers$avg_log2FC) >= log2FC)
p2 <- p1 + geom_text_repel(
  data = subset(cluster1.de.markers, cluster1.de.markers$p_val_adj < padj & abs(cluster1.de.markers$avg_log2FC) >= log2FC),
  aes(label = rownames(data)), size = 3,
  box.padding = unit(0.5, "lines"),
  point.padding = unit(0.8,"lines"), segment.color = "black", show.legend = FALSE )
p2

cluster1.de.genes <- rownames(cluster1.de.markers)
write.csv(cluster1.de.genes, file = "data/downstream/enrichment/cluster1.de.genes.csv")

##################Enrichment##################
library(stringi)
library(ggplot2)
library(dplyr)

# KEGG
downgokegg<-read.delim("data/downstream/enrichment/KEGG.txt")
enrich<-downgokegg
enrich_signif=enrich[which(enrich$PValue<0.05),]
enrich_signif=enrich_signif[,c(1:3,5)]
head(enrich_signif)
enrich_signif=data.frame(enrich_signif)
KEGG=enrich_signif
KEGG$Term<-stri_sub(KEGG$Term,10,100)
ggplot(KEGG,aes(x=Count,y=Term))+geom_point(aes(color=PValue,size=Count))+scale_color_gradient(low='slateblue4',high='firebrick3')+theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())

# GO_CC
GO_CC<-read.delim("data/downstream/enrichment/GO_CC.txt")
GO_CC_signif=GO_CC[which(GO_CC$PValue<0.05),]
GO_CC_signif=GO_CC[,c(1:3,5)]
head(GO_CC_signif)
GO_CC_signif=data.frame(GO_CC_signif)
GO_CC_signif$Term<-stri_sub(GO_CC_signif$Term,12,100)

# GO_BP
GO_BP<-read.delim("data/downstream/enrichment/GO_BP.txt")
GO_BP_signif=GO_BP[which(GO_BP$PValue<0.05),]
GO_BP_signif=GO_BP_signif[,c(1:3,5)]
head(GO_BP_signif)
GO_BP_signif=data.frame(GO_BP_signif)
GO_BP_signif$Term<-stri_sub(GO_BP_signif$Term,12,100)

# GO_MF
GO_MF<-read.delim("data/downstream/enrichment/GO_MF.txt")
GO_MF_signif=GO_MF[which(GO_MF$PValue<0.05),]
GO_MF_signif=GO_MF_signif[,c(1:3,5)]
head(GO_MF_signif)
GO_MF_signif=data.frame(GO_MF_signif)
GO_MF_signif$Term<-stri_sub(GO_MF_signif$Term,12,100)
enrich_signif=rbind(GO_BP_signif,rbind(GO_CC_signif,GO_MF_signif))
go=enrich_signif
go=arrange(go,go$Category,go$PValue)

##图例名称设置
m=go$Category
m=gsub("TERM","",m)
m=gsub("_DIRECT","",m)
go$Category=m
GO_term_order=factor(as.integer(rownames(go)),labels = go$Term)
COLS<-c("#66C3A5","#8DA1CB","#FD8D62")

###开始画图
ggplot(data=go,aes(x=GO_term_order,y=Count,fill=Category))+
  geom_bar(stat = "identity",width = 0.8)+
  scale_fill_manual(values = COLS)+
  theme_bw()+
  xlab("Terms")+
  ylab("Gene_counts")+
  labs()+
  theme(axis.text.x = element_text(face = "bold",color = "black",angle = 90,vjust = 1,hjust = 1)) 


##################统计比例##################
library(Seurat)
library(stringi)
library(rtracklayer)
source("coding/Deconvolution_functions.R")
dataSC <- readRDS("data/downstream/markers/tcells_clutser_markers.rds")
dataSC.count <- GetAssayData(dataSC, assay = "RNA", layer = "counts") # get count matrix data
dataSC.label <- dataSC@meta.data$seurat_clusters
dataBulk <- readRDS("data/downstream/scissor_data/bulk_dataset.rds")

Signature<-buildSignatureMatrixMAST (scdata=dataSC.count, id=dataSC.label, 
                                     path="data/downstream/results", diff.cutoff=0.5, pval.cutoff=0.01)
saveRDS(Signature, file = "data/downstream/results/signature.rds")

#trim signature and bulk data to contain the same differentially expressed genes
tr<-trimData(Sig,dataBulk)
#estimate using dampened weighted least squares
solDWLS <- data.frame()
for (i in 1:ncol(tr$bulk)){
  tmp<-solveDampenedWLS(tr$sig,tr$bulk[,i])
  tmp <- as.data.frame(tmp)
  colnames(tmp) <- colnames(dataBulk)[i]
  if(i == 1){
    solDWLS <- tmp
  }else{
    solDWLS <- cbind(solDWLS, tmp)
  }
}
saveRDS(solDWLS, file = "data/downstream/results/solDWLS.rds")


##################COX生存分析-DWLS##################
library(survival)
library(dplyr)
library(survminer)
solDWLS <- readRDS("data/downstream/results/solDWLS.rds")
#naive_percent <- solDWLS["1", ]
clinical <- readRDS("data/raw_TCGA/clinical.rds")
clinical <- clinical[, c("age_at_diagnosis", "gender", "time", "event")]
#clinical <- clinical[complete.cases(clinical), ]
names(clinical) <- c("age", "sex", "time", "status")
#clinical <- rename(clinical, age = age_at_diagnosis, sex = gender, status = event, time = time)

# 1. 转置 solDWLS 数据框
solDWLS_transposed <- t(solDWLS)
colnames(solDWLS_transposed) <- rownames(solDWLS)

# 2. 匹配操作，将两个数据框中具有相同ID的样本匹配起来
matched_samples <- intersect(rownames(clinical), rownames(solDWLS_transposed))
clinical_matched <- clinical[matched_samples, ]
solDWLS_transposed_matched <- solDWLS_transposed[matched_samples, ]

# 3. 对匹配后的数据进行排序，确保两个数据框中的ID顺序一致
clinical_matched_sorted <- clinical_matched[order(rownames(clinical_matched)), ]
solDWLS_transposed_matched_sorted <- solDWLS_transposed_matched[order(rownames(solDWLS_transposed_matched)), ]

# 4. 合并临床数据和Naive细胞类型比例数据
merged_data <- data.frame(clinical_matched_sorted, solDWLS_transposed_matched_sorted)
rownames(merged_data) <- rownames(clinical_matched_sorted)

# 选择所需的列进行分析，分析Naive细胞类型对预后的影响
naive_survival_data <- merged_data[, c("age", "time", "status", "X3")]
#naive_survival_data <- merged_data[, c("age", "time", "status", "X10")]
colnames(naive_survival_data)[colnames(naive_survival_data) == "X3"] <- "naive"
# 对naive_survival_data中的naive列进行从高到低的排序
naive_survival_data <- naive_survival_data[order(naive_survival_data$naive, decreasing = TRUE), ]
# 创建高低两组标签
naive_survival_data$naive_group <- ifelse(rank(naive_survival_data$naive) <= nrow(naive_survival_data) / 2, "low", "high")

# 创建 Cox 比例风险模型
cox_model_naive <- coxph(Surv(time, status) ~ naive_group+age, data = naive_survival_data)
summary(cox_model_naive)

# 绘制Kaplan-Meier曲线
fit <- survfit(Surv(time, status) ~ naive_group, data = naive_survival_data)
ggsurvplot(fit, data = naive_survival_data, risk.table = TRUE, palette = "jco",
           title = "Kaplan-Meier Survival Curve by Naive Group",
           xlab = "Time (days)", ylab = "Survival Probability")

# 创建线性回归模型
lm_age_model <- lm(naive ~ age, data = naive_survival_data)

# 查看模型摘要
summary(lm_age_model)

# 创建散点图和拟合线性回归线
plot(naive_survival_data$age, naive_survival_data$naive, xlab = "Age", ylab = "Naive proportion", main = "Relationship between Age and Naive proportion")
abline(lm_age_model, col = "red")

##################lasso-cox##################
library(Seurat)
library(glmnet)
library(survival)
library(DESeq2)
# 读取聚类标记数据
cluster_markers <- readRDS("data/downstream/markers/tcells_clutser_markers.rds")
cluster1.de.markers <- FindMarkers(cluster_markers, ident.1 = "3.Naive", ident.2 = NULL, only.pos = F, min.pct = 0.25)

#自定义阈值
log2FC <-  0.5
padj <-  0.01 

cluster1.de.markers$threshold="ns";
cluster1.de.markers[which(cluster1.de.markers$avg_log2FC  > log2FC & cluster1.de.markers$p_val_adj <padj),]$threshold="up";
cluster1.de.markers[which(cluster1.de.markers$avg_log2FC  < (-log2FC) & cluster1.de.markers$p_val_adj < padj),]$threshold="down";
cluster1.de.markers$threshold=factor(cluster1.de.markers$threshold, levels=c('down','ns','up'))
de_genes <-  subset(cluster1.de.markers, cluster1.de.markers$p_val_adj < padj & abs(cluster1.de.markers$avg_log2FC) >= log2FC)

# 匹配差异基因与bulk，得到差异基因表达谱
bulkdata <- readRDS("data/downstream/scissor_data/bulk_dataset.rds")
#match_order <- match(rownames(bulkdata),rownames(de_genes))
#bulkdata_matched <- bulkdata[match_order, , drop = FALSE]  # 使用 drop = FALSE 以保持数据框的结构
#bulkdata_matched <- na.omit(bulkdata_matched)  # 删除包含缺失值的行
common_genes <- intersect(rownames(de_genes),rownames(bulkdata_matched))
bulkdata.de.genes <- bulkdata_matched[common_genes,]

# 加载TCGA临床数据
clinical <- readRDS("data/raw_TCGA/clinical.rds")
clinical <- clinical[, c("age_at_diagnosis", "gender", "time", "event")]
names(clinical) <- c("age", "sex", "time", "status")
#clinical <- clinical[complete.cases(clinical), ]
clinical$status <- as.double(clinical$status)
clinical$time <- as.double(clinical$time)
clinical <- na.omit(clinical)

# 匹配差异基因表达谱与clinical
#match_order_names <- match(colnames(bulkdata.de.genes),rownames(clinical))
#bulkdata.de.genes <- bulkdata.de.genes[match_order_names]
#clinical.de.genes <- clinical[match_order_names]
common_names <- intersect(rownames(clinical), colnames(bulkdata.de.genes))
clinical.de.genes <- clinical[common_names, ]
bulkdata.de.genes <- bulkdata.de.genes[ , common_names]
bulkdata.de.genes <- t(bulkdata.de.genes)
all(rownames(bulkdata.de.genes) == rownames(clinical.de.genes))
x <- bulkdata.de.genes
y <- data.matrix(Surv(time = clinical.de.genes$time, 
                      event = clinical.de.genes$status))
rownames(y) <- rownames(clinical.de.genes)

fit = glmnet(x,y,
             alpha = 1,family = "cox") 
plot(fit,xvar= 'lambda',label = TRUE)

fit <- glmnet(x, y, family = 'cox', type.measure = "deviance", nfolds = 10)
#十折交叉检验筛选最佳lambda：
set.seed(007)
lasso_fit <- cv.glmnet(x, y, family = 'cox', type.measure = 'deviance', nfolds = 10)
plot(lasso_fit)

#提取最佳λ值(这里选择1se对应lambda)：
lambda.1se <- lasso_fit$lambda.1se
lambda.1se

#使用1se的lambda重新建模：
model_lasso_1se <- glmnet(x, y, family = 'cox', type.measure = 'deviance', nfolds = 10,lambda = lambda.1se)
#拎出建模使用基因：
gene_1se <- rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]#as.numeric后"."会转化为0
gene_1se #筛选出7个

##################Lasso + Cox 生存分析模式##################

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("glmnet")
library(glmnet)

### 1. 载入数据
data("CoxExample")
class(CoxExample) # [1] "list"

dim(CoxExample$x)
dim(CoxExample$y)

## 加上假基因名，样本名
rownames(CoxExample$x) <- paste0('sample',1:dim(CoxExample$x)[1])
colnames(CoxExample$x) <- paste0('gene',1:dim(CoxExample$x)[2])

rownames(CoxExample$y) <- paste0('sample',1:dim(CoxExample$x)[1])


head(CoxExample$x) # 行：患者/样本，列：特征，可以是特征基因的表达谱
head(CoxExample$y) # 生存数据，行：患者/样本，列：time, status

### 2. Lasso Cox模型
fit = glmnet(CoxExample$x,CoxExample$y,
             alpha = 1,family = "cox") 
## alpha = 1，用于变量选择，去除系数为零的变量
## lambda =0, 没有正则项，等于coxph_fit
plot(fit,xvar= 'lambda',label = TRUE)

## 交叉验证
# 默认 nfolds = 10
cvfit <- cv.glmnet(CoxExample$x,CoxExample$y,
                   alpha = 1,family = "cox")
# type.measure :deviance

# cvfit <- cv.glmnet(CoxExample$x,CoxExample$y,
#                    alpha = 1,family = "cox",
#                    type.measure= "C")

plot(cvfit)

### 3. 筛选特征
coefficient <- coef(cvfit, s= cvfit$lambda.min)
selected_index <- which(as.numeric(coefficient) != 0)
selected_features <- names(coefficient[selected_index,])

### 4. 预测生存时间
predict(cvfit, newx = CoxExample$x[1:10, ], 
        type = "response")
# type:“link”, “response”, “coefficients”, “nonzero”其中的一个

### 5. 计算risk score 
# 一般选择risk score计算公式

risk_score <- apply(CoxExample$x[1:10,selected_features],1,
                    function(x) sum(x*coefficient[selected_index,]))

# sum(CoxExample$x[10,selected_features]*coefficient[selected_index,])

##################COX生存分析-BayesPrism##################
library(survival)
library(dplyr)
library(survminer)
cluster.percent <- read.csv("data/downstream/thera.first.state.csv")
clinical <- readRDS("data/raw_TCGA/clinical.rds")
clinical <- clinical[, c("age_at_diagnosis", "gender", "time", "event")]
names(clinical) <- c("age", "sex", "time", "status")
clinical <- na.omit(clinical)

# 取临床数据与细胞比例数据的交集
matched_samples <- intersect(rownames(clinical), cluster.percent$X)
clinical <- clinical[matched_samples, ]
rownames(cluster.percent) <- cluster.percent$X
cluster.percent <- cluster.percent[matched_samples, ]

# 对匹配的数据进行排序
clinical <- clinical[order(rownames(clinical)), ]
cluster.percent <- cluster.percent[order(rownames(cluster.percent)), ]

# 合并数据
merged_data <- data.frame(clinical, cluster.percent)
rownames(merged_data) <- rownames(clinical)

# 选择所需的列进行分析，分析Naive细胞类型对预后的影响
naive_survival_data <- merged_data[, c("age", "time", "status", "sex", "X3")]
#naive_survival_data <- merged_data[, c("age", "time", "status", "sex", "X10")]
colnames(naive_survival_data)[colnames(naive_survival_data) == "X3"] <- "naive"
# 对naive_survival_data中的naive列进行从高到低的排序
naive_survival_data <- naive_survival_data[order(naive_survival_data$naive, decreasing = TRUE), ]
# 创建高低两组标签
naive_survival_data$naive_group <- ifelse(rank(naive_survival_data$naive) <= nrow(naive_survival_data) / 2, "low", "high")

# 创建 Cox 比例风险模型，并绘制森林图
cox_model_naive <- coxph(Surv(time, status) ~ naive_group + age + sex, data = naive_survival_data)
summary(cox_model_naive)
ggforest(cox_model_naive, data = naive_survival_data)

# 绘制Kaplan-Meier曲线
fit <- survfit(Surv(time, status) ~ naive_group, data = naive_survival_data)
surv_pvalue(fit)$pval.txt

ggsurvplot(fit, data = naive_survival_data, risk.table = TRUE, palette = "jco",
                   title = "Kaplan-Meier Survival Curve by Naive Group",
                   xlab = "Time (days)", ylab = "Survival Probability",)
#          pval = surv_pvalue(fit)$pval.txt)
# 保存图片
ggsave("results/image/survival_curve.png")

fit <- coxph(Surv(time, status)) ~ sex + age

coxph(formula = Surv(time, status) ~ sex + age +naive_group, data = naive_survival_data)

# 创建线性回归模型
lm_age_model <- lm(naive ~ age, data = naive_survival_data)
summary(lm_age_model)

# 创建散点图和拟合线性回归线
plot(naive_survival_data$age, naive_survival_data$naive, xlab = "Age", ylab = "Naive proportion", main = "Relationship between Age and Naive proportion")
abline(lm_age_model, col = "red")

##################统计scissor.neg占比##################
library(Seurat)
data.sc <- readRDS("data/downstream/markers/tcells_clutser_markers.rds")
sc.scissor.neg <- subset(data.sc, subset = scissor == "Scissor- cell")
# 计算所有细胞中scissor- cell占比：4.81
print(paste("所有细胞中Scissor- cell 占比：", 
            total.scissor.neg.percent <- ncol(sc.scissor.neg) / ncol(data.sc) * 100, "%"))

# 计算naive中scissor- cell占比:18.89
sc.naive <- subset(data.sc, subset = seurat_clusters == "3")
sc.naive.scissor.neg <- subset(sc.naive, subset = scissor == "Scissor- cell")
print(paste("naive中Scissor- cell 占比：", 
            naive.scissor.neg.percent <- ncol(sc.naive.scissor.neg) / ncol(sc.naive) * 100, "%"))

# 统计不同样本中scissor- cell占比
data.sc <- readRDS("data/downstream/markers/tcells_clutser_markers.rds")
metadata <- read.csv("data/download/sc_dataset/Meta_GBM.txt.gz")
metadata <- metadata[metadata$Assignment %in% "TCells",] #meta.tcells
metadata <- metadata[metadata$NAME %in% rownames(data.sc@meta.data), ] # meta.matched
# table(metadata$donor_id)
# sample.id <- metadata[metadata$donor_id %in% "human_ndGBM-11",][, 1]
# sample <- subset(data.sc, barcode %in% sample.id)
# sample.scissor.neg <- subset(sample, subset = scissor == "Scissor- cell")
# print(paste("样本中Scissor- cell 占比：", 
#             total.scissor.neg.percent <- ncol(sample.scissor.neg) / ncol(sample) * 100, "%"))

data.sc@meta.data$patient <- metadata$Patient
cells.tcga <- data.sc@meta.data[data.sc@meta.data$seurat_clusters == 1, c( "scissor", "patient")]
pct.scossor <- table(data.sc@meta.data[data.sc@meta.data$seurat_clusters == 1, c( "scissor", "patient")])
pct.scossor <- rbind(pct.scossor, colSums(pct.scossor))
pct.scossor <- t(pct.scossor)
pct.scossor <- as.data.frame(pct.scossor)
pct.scossor$pct <- pct.scossor$`Scissor- cell`/pct.scossor$V4

write.csv(pct.scossor, "results/pct.tcells.scissor.neg.csv", row.names = F)
saveRDS(cells.orig, file = "data/downstream/scissor_data/scissor.cells.orig.rds")

##################for循环##################
data.sc <- readRDS("data/downstream/markers/tcells_clutser_markers.rds") # sc.tcells
#data.sc <- readRDS("data/seurat/glioma_v5.rds") #sc.total.glioma
metadata <- read.csv("data/download/sc_dataset/Meta_GBM.txt.gz") #meta.total
#metadata <- metadata[metadata$Assignment %in% "TCells",] #meta.tcells

# 获取所有样本的 donor_id
sample_ids <- unique(metadata$donor_id)
sample_ids <- sample_ids[sample_ids != "group"]
numeric_suffix <- as.numeric(gsub("^.*-(\\d+)$", "\\1", sample_ids))
sample_ids <- sample_ids[order(numeric_suffix)]
# 初始化一个空数据框用于保存结果
results <- data.frame(sample_id = character(),
                      scissor_neg_percent = numeric(),
                      stringsAsFactors = FALSE)
# 循环遍历每个样本
for (sample_id in sample_ids) {
  # 从元数据中选择当前样本的 donor_id
  # sample_id = "human_LGG-03"
  # sample_id = "human_LGG-04"
  
  # sample_id = " human_rGBM-01"
  # sample_id = " human_rGBM-02"
  # sample_id = " human_rGBM-03"
  # sample_id = " human_rGBM-04"
  # sample_id = " human_rGBM-05"
  
  # sample_id = "human_ndGBM-01"
  # sample_id = "human_ndGBM-02"
  # sample_id = "human_ndGBM-03"
  # sample_id = "human_ndGBM-04"
  # sample_id = "human_ndGBM-05"
  # sample_id = "human_ndGBM-06"
  # sample_id = "human_ndGBM-07"
  # sample_id = "human_ndGBM-08"
  # sample_id = "human_ndGBM-09"
  # sample_id = "human_ndGBM-10"
  # sample_id = "human_ndGBM-11"

  current_sample.id <- metadata[metadata$donor_id == sample_id, ][, 1]
  
  # 提取当前样本数据
  current_sample <- subset(data.sc, barcode %in% current_sample.id)
  
  # 从当前样本中提取scissor- cell
  current_sample.scissor.neg <- subset(current_sample, subset = scissor == "Scissor- cell")
  
  # 计算当前样本中 scissor- cell 的占比
  scissor_neg_percent <- ncol(current_sample.scissor.neg) / ncol(current_sample) * 100
  
  # 将结果添加到结果数据框中
  results <- rbind(results, data.frame(sample_id = sample_id,
                                       scissor_neg_percent = scissor_neg_percent))
  results <- unique(results)
}
write.csv(results, file = "results/scissor.neg.percent_total.csv", row.names = FALSE)
##################循环遍历三组样本##################
data.sc <- readRDS("data/downstream/markers/tcells_clutser_markers.rds") # sc.tcells
metadata <- read.csv("data/download/sc_dataset/Meta_GBM.txt.gz") #meta.total
#metadata <- metadata[metadata$Assignment %in% "TCells",] #meta.tcells

sample_ids <- unique(metadata$donor_id)
sample_ids <- sample_ids[sample_ids != "group"]

# 提取末尾数字,按照末尾数字排序
numeric_suffix <- as.numeric(gsub("^.*-(\\d+)$", "\\1", sample_ids))
sample_ids <- sample_ids[order(numeric_suffix)]
# 分成三组
LGG_samples <- sample_ids[grep("LGG", sample_ids)]
rGBM_samples <- sample_ids[grep("rGBM", sample_ids)]
ndGBM_samples <- sample_ids[grep("ndGBM", sample_ids)]

results <- data.frame(sample_id = character(),
                      scissor_neg_percent = numeric(),
                      stringsAsFactors = FALSE)
for (group in list(LGG_samples, rGBM_samples, ndGBM_samples)) {
  for (sample_id in group) {
    
    # 从元数据中选择当前样本的 donor_id
    current_sample.id <- metadata[metadata$donor_id == sample_id, ][, 1]
    
    # 提取当前样本数据
    current_sample <- subset(data.sc, barcode %in% current_sample.id)
    
    # 从当前样本中提取 scissor- cell
    current_sample.scissor.neg <- subset(current_sample, subset = scissor == "Scissor- cell")
    
    # 计算当前样本中 scissor- cell 的占比
    scissor_neg_percent <- ncol(current_sample.scissor.neg) / ncol(current_sample) * 100
    
    # 将结果添加到结果数据框中
    results <- rbind(results, data.frame(sample_id = sample_id,
                                         scissor_neg_percent = scissor_neg_percent))
    results <- unique(results)
  }
}
write.csv(results, file = "results/scissor.neg.percent_total.csv", row.names = FALSE)
##################glioma_统计total.scissor.neg占比##################
library(Seurat)
data.sc.glioma <- readRDS("data/seurat/glioma_v5.rds")

# count_dataset <- table(data.sc.glioma@meta.data$orig.ident)
# count_few <- count_dataset[count_dataset < 30]
# tmp <- subset(tmp, subset = orig.ident %in% names(count_few), invert = T)
data.sc.glioma[["RNA"]] <- split(data.sc.glioma[["RNA"]], f = data.sc.glioma$orig.ident)
tmp <- data.sc.glioma
tmp <- NormalizeData(tmp)
tmp <- FindVariableFeatures(tmp)
tmp <- ScaleData(tmp)
tmp <- RunPCA(tmp)

options(future.globals.maxSize = 5 * 1024^3)  # 设置future全局变量大小阈值为5 GiB
tmp <- IntegrateLayers(
  object = tmp, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = F, k.weight = 32)
saveRDS(tmp, file = "data/seurat/glioma.rm.batch.rds")

tmp <- readRDS("data/seurat/glioma.rm.batch.rds")
# re-join layers after integration
tmp[["RNA"]] <- JoinLayers(tmp[["RNA"]])
tmp <- FindNeighbors(tmp, reduction = "integrated.rpca", dims = 1:30)
tmp <- FindClusters(tmp, resolution = 0.8)
tmp <- RunUMAP(tmp, dims = 1:30, reduction = "integrated.rpca")
saveRDS(tmp, file = "data/seurat/glioma/glioma.integrated.rpca.rds")

# Scissor
tmp <- readRDS("data/seurat/glioma/glioma.integrated.rpca.rds")
#tmp <- readRDS("data/seurat/glioma_v5.rds")
options(Seurat.object.assay.version = 'v4')
GBM_sc_dataset <- GetAssayData(object = tmp, assay = "RNA", slot = "counts") # get count matrix data
# sc_dataset <- as(GBM_sc_dataset, "sparseMatrix")
sc_dataset <- Scissor::Seurat_preprocessing(GBM_sc_dataset, verbose = T)
saveRDS(sc_dataset,file = "data/downstream/scissor_data/glioma.sc_dataset.rds")


bulk_dataset <- readRDS("data/downstream/scissor_data/bulk_dataset.rds")
phenotype <- readRDS("data/downstream/scissor_data/phenotype.rds")
sc_dataset <- readRDS("data/downstream/scissor_data/glioma.sc_dataset.rds")

# sc_dataset[["RNA3"]] <- as(object = sc_dataset[["RNA"]], Class = "Assay")
# DefaultAssay(sc_dataset) <- 'RNA3'
# sc_dataset[['RNA']]<-NULL
# sc_dataset <- RenameAssays(sc_dataset,assay.name = 'RNA3',new.assay.name = 'RNA')

infos <- Scissor::Scissor(bulk_dataset, sc_dataset, phenotype, alpha = 0.05,
                           family = "cox")

DimPlot(tmp, reduction = "umap", split.by = "scissor", label.size = 6, label = T, pt.size = 1)
saveRDS(tmp, file = "data/downstream/tcells_resolution=0.8.rds")

# find markers
tmp <- readRDS("data/downstream/tcells_resolution=0.8.rds")
tcells.markers <- FindAllMarkers(tmp, only.pos = T)
tcells.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

#DoHeatmap(tmp, features = top10$gene) + NoLegend()
saveRDS(tcells.markers, file = "data/downstream/markers/tcells.markers.rds")
saveRDS(top10, file = "data/downstream/markers/top10.rds")
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

data.sc <- readRDS("data/seurat/tcells_v5.rds")
data.sc <- GetAssayData(object = data.sc, assay = "RNA", slot = "counts") # get count matrix data
data.sc <- Scissor::Seurat_preprocessing(data.sc, verbose = F)
# 重命名Assay，切换v3版本
data.sc <- RenameAssays(object = data.sc, assay.name = "RNA", new.assay.name = "RNA5")
data.sc[["RNA"]] <- as(object = data.sc[["RNA5"]], Class = "Assay")
DefaultAssay(data.sc) <- "RNA"
##sc_dataset@assays$RNA <- NULL

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
DimPlot(data.sc,reduction = 'umap',group.by = 'scissor',cols = c('grey','royalblue','indianred1'),pt.size = 2)
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
##################CGGA-325##################
library(Seurat)
library(Matrix)
data.bulk <- read.table("data/download/sc.data.CGGA/mRNAseq_325/bulkdata/CGGA.mRNAseq_325.RSEM-genes.20200506.txt", 
                        header = T, 
                        sep = "",
)
data.bulk <- `rownames<-`(data.bulk[,-1], data.bulk$Gene_Name)
data.bulk <- as.matrix(data.bulk)

data.clinical <- read.table("data/download/sc.data.CGGA/mRNAseq_325/clinical/CGGA.mRNAseq_325_clinical.20200506.txt",
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
# data.sc <- GetAssayData(data.sc, assay = "RNA", slot = "counts") # get count matrix data
# data.sc <- Scissor::Seurat_preprocessing(data.sc, verbose = F)
# # 重命名Assay
# data.sc <- RenameAssays(data.sc, assay.name = "RNA", new.assay.name = "RNA5")
# data.sc[["RNA"]] <- as(data.sc[["RNA5"]], Class = "Assay")
# DefaultAssay(data.sc) <- "RNA"
# ##sc_dataset@assays$RNA <- NULL

data.sc <- readRDS("data/downstream/scissor_data/CGGA.325/data.sc.rds")
# Scissor  pct.scissor:15.856% | neg cells:1295 | pos cells:1808
info.325 <- Scissor::Scissor(data.bulk, data.sc, data.clinical.gbm, alpha = 0.05,
                         family = "cox")

saveRDS(data.bulk, file = "data/downstream/scissor_data/CGGA.325/data.bulk.rds")
saveRDS(data.clinical.gbm, file = "data/downstream/scissor_data/CGGA.325/data.clinical.rds")
saveRDS(data.sc, file = "data/downstream/scissor_data/CGGA.325/data.sc.rds")
saveRDS(info.325, file = "data/downstream/scissor_data/CGGA.325/scissor.info.325.rds")

data.bulk <- readRDS("data/downstream/scissor_data/CGGA.325/data.bulk.rds")
data.clinical <- readRDS("data/downstream/scissor_data/CGGA.325/data.clinical.rds")
data.sc <- readRDS("data/downstream/scissor_data/CGGA.325/data.sc.rds")
info.325 <- readRDS("data/downstream/scissor_data/CGGA.325/scissor.info.325.rds")

# add metadata.scissor
Scissor_select <- rep("Background cells",ncol(data.sc))
names(Scissor_select) <- colnames(data.sc)
Scissor_select[info.325$Scissor_pos] <- "Scissor+ cell"
Scissor_select[info.325$Scissor_neg] <- "Scissor- cell"
data.sc <- AddMetaData(data.sc,metadata = Scissor_select,col.name = "scissor")
DimPlot(data.sc,reduction = 'umap',group.by = 'scissor',cols = c('grey','royalblue','indianred1'),pt.size = 2)

saveRDS(data.sc,file = "data/downstream/scissor_data/CGGA.325/data.sc.scissor.rds")

# match data.sc.tcells和data.sc.scissor
data.sc.tcells <- readRDS("data/downstream/markers/tcells_clutser_markers.rds")
data.sc.scissor <- readRDS("data/downstream/scissor_data/CGGA.325/data.sc.scissor.rds")
scissor_sub <- data.sc.scissor@meta.data
scissor_sub <- scissor_sub[match(colnames(data.sc.tcells), rownames(scissor_sub)),]
data.sc.tcells@meta.data$scissor <- scissor_sub$scissor

saveRDS(data.sc.tcells, file = "data/downstream/scissor_data/CGGA.325/data.sc.scissor.tcells")
##################pct.325.scissor.neg##################
library(Seurat)
data.sc <- readRDS("data/downstream/scissor_data/CGGA.325/data.sc.scissor.tcells")
sc.scissor.neg <- subset(data.sc, subset = scissor == "Scissor- cell")
# 计算所有细胞中scissor- cell占比：4.81 | 6.59
print(paste("所有细胞中Scissor- cell 占比：", 
            total.scissor.neg.percent <- ncol(sc.scissor.neg) / ncol(data.sc) * 100, "%"))

# 计算naive中scissor- cell占比:18.89 |10.07
sc.naive <- subset(data.sc, subset = seurat_clusters == "3")
sc.naive.scissor.neg <- subset(sc.naive, subset = scissor == "Scissor- cell")
print(paste("naive中Scissor- cell 占比：", 
            naive.scissor.neg.percent <- ncol(sc.naive.scissor.neg) / ncol(sc.naive) * 100, "%"))

# 统计不同样本中scissor- cell占比
data.sc <- readRDS("data/downstream/scissor_data/CGGA.325/data.sc.scissor.tcells")
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
cells.325 <- data.sc@meta.data[data.sc@meta.data$seurat_clusters == 1, c( "scissor", "patient")]
pct.scissor <- table(data.sc@meta.data[data.sc@meta.data$seurat_clusters == 1, c( "scissor", "patient")])
pct.scissor <- rbind(pct.scissor, colSums(pct.scissor))
pct.scissor <- t(pct.scissor)
pct.scissor <- as.data.frame(pct.scissor)
pct.scissor$pct <- pct.scissor$`Scissor- cell`/pct.scissor$V4
write.csv(pct.scissor, "results/pct.325.tcells.scissor.neg.csv", row.names = F)
saveRDS(cells.325, file = "data/downstream/scissor_data/CGGA.325/scissor.cells.325.rds")

##################intersect##################
#cells <- data.sc@meta.data[data.sc@meta.data$seurat_clusters == 1, c( "scissor", "patient")]
cells.tcga <- readRDS("data/downstream/scissor_data/scissor.cells.tcga.rds")
cells.693 <- readRDS("data/downstream/scissor_data/CGGA.693/scissor.cells.693.rds")
#cells.325 <- readRDS("data/downstream/scissor_data/CGGA.325/scissor.cells.325.rds")

#intersect(rownames(cells.tcga[cells.tcga$scissor == "Scissor- cell",]), 
#          rownames(cells.693[cells.693$scissor == "Scissor- cell",]))
intersect(rownames(cells.tcga[cells.tcga$scissor == "Scissor- cell",]), 
          rownames(cells.693[cells.693$scissor == "Scissor- cell",]))


# 按样本ID
scissor_cells_orig <- cells.tcga$patient[cells.tcga$scissor == "Scissor- cell"]
scissor_cells_693 <- cells.693$patient[cells.693$scissor == "Scissor- cell"]
scissor_cells_325 <- cells.325$patient[cells.325$scissor == "Scissor- cell"]
intersection.by.patient.id <- intersect(intersect(scissor_cells_orig, scissor_cells_693), scissor_cells_325)

# 按细胞ID（行名）
scissor_cells_orig <- rownames(cells.tcga)[cells.tcga$scissor == "Scissor- cell"]
scissor_cells_693 <- rownames(cells.693)[cells.693$scissor == "Scissor- cell"]
scissor_cells_325 <- rownames(cells.325)[cells.325$scissor == "Scissor- cell"]
intersection.by.cell.id <- intersect(intersect(scissor_cells_orig, scissor_cells_693), scissor_cells_325)


