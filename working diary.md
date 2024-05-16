# 2024/05/13-2024/0519

## 2024/05/16

### 整理基于R语言的单细胞数据分析--preprocessing部分

#### 1.实验目的

2.1熟悉基于R的seurat标准预处理分析流程

2.2理解seurat obeject数据结构

2.3为后续细胞注释分析工作作准备

#### 3.实验材料：

##### 3.1数据来源：

Single-cell analysis of human glioma and immune cells identifies S100A4 as an immunotherapy target(GSE182109)

Single Cell Portal数据，人，脑胶质瘤GBM，T cells（subset）

链接：https://singlecell.broadinstitute.org/single_cell/study/SCP1985/single-cell-analysis-of-human-glioma-and-immune-cells-identifies-s100a4-as-an-immunotherapy-target-gse182109

##### 3.2软件与平台：

R（v4.3.3）；RStudio；Seurat（v5.0.3，https://github.com/satijalab/seurat）；Seurat tutorial（https://satijalab.org/seurat/）

#### 4.实验步骤

简单介绍预处理流程

##### 4.1Load raw data

1.read10X / readhd5 / readRDS

2.subset my interested cell type：T cells

##### 4.2Quality control

1.add percent.mt

2.subset according to percent.mt , nFeature_RNA and nCount_RNA

##### 4.3Standard work flow

1.Normalization

2.FindVariableFeatures

3.Scale

4.PCA

5.FindNeighbours

6.FindClusters (we can setup resolution e.g. res=0.1, 0.8)

7.UMAP

7.Finally, we will get a clustering map

Remember: if the data comes from muti datasets, do integration after PCA(e.g. CCA-Integration, RPCA-Integration)
