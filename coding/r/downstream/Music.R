devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")

workdir = "."  # where you place the data and results
setwd(workdir)
library(dplyr)

####################################### Meta data
# load sc.data for T cells, and drop non-selected T cells from the primary dataset
sc.tcells <- readRDS("tmp/dataSC.rds")
sc.tcell.cellId <- colnames(sc.tcells)
#rm(sc.tcells)

# Load meta data recording cluster of per cell
sc.meta <- read.csv(gzfile("raw/download/Meta_GBM.txt.gz"), sep = ",", header = TRUE)
sc.meta <- sc.meta[-1, ]
sc.meta.tcell <- sc.meta[sc.meta$Assignment == "TCells",]
sc.meta.tcell.drop <- sc.meta.tcell[-match(sc.tcell.cellId, sc.meta.tcell$NAME), 1]
rm(sc.meta.tcell)
# Filter dropped T cells
sc.meta <- sc.meta[!sc.meta$NAME %in% sc.meta.tcell.drop, ]
# get subclusters for T cells
sc.tcell.meta <- sc.tcells@meta.data
sc.tcell.meta$NAME <- row.names(sc.tcell.meta)

###################################### scRNA data
sc.data <- Seurat::Read10X("./raw/download/sc_dataset/")
sc.meta.cluster <- merge(sc.meta, sc.tcell.meta, by = "NAME", all = TRUE)
sc.meta.cluster <- sc.meta.cluster %>%
  mutate(Assignment = case_when(Assignment == "TCells" ~ paste("X", seurat_clusters, sep = ""),
                                .default = as.character(Assignment)))
# subset sc.data
sc.data <- sc.data[, colnames(sc.data) %in% sc.meta.cluster$NAME]
# re-order sc.meta.cluster
sc.meta.cluster <- sc.meta.cluster[match(sc.meta.cluster$NAME, colnames(sc.data)), ]

rm(sc.tcell.meta, sc.tcells, random, sc.tcell.cellId, sc.meta.tcell.drop, sc.meta)
saveRDS(sc.data, "./tmp/sc.data.filter.RDS")
write.csv(sc.meta.cluster, "./tmp/sc.meta.cluster.csv", row.names = FALSE, quote = FALSE)

###################################### optical: randomly sample
# random sub-sample
sc.data <- readRDS("./tmp/sc.data.filter.RDS")
set.seed(100)
random <- sample(1:ncol(sc.data), 80000, replace = FALSE)
sc.data.subset <- sc.data[, random]

################################################################################
# bulk RNAseq data
bulk.data <- readRDS("tmp/dataBulk.rds")
bulk.data <- t(bulk.data)

sc.data.subset <- MatrixExtra::t_shallow(sc.data.subset)
sc.data.subset <- MatrixExtra::as.csc.matrix(sc.data.subset, binary = FALSE, logical = FALSE, sort = FALSE)

sc.meta.cluster <- read.csv("./tmp/sc.meta.cluster.csv", header = TRUE)
rownames(sc.meta.cluster) <- sc.meta.cluster$NAME
sc.meta.subset <- sc.meta.cluster[match(rownames(sc.data.subset), sc.meta.cluster$NAME), "Assignment"]
sc.meta.state <- sc.meta.subset
sc.meta.subset[grepl("^X", sc.meta.subset)] <- "TCell" 


sc.dat.filtered <- BayesPrism::cleanup.genes(input=sc.data.subset,
                                  input.type="count.matrix",
                                  species="hs", 
                                  gene.group=c("Rb","Mrp","other_Rb","chrM",
                                               "MALAT1","chrX","chrY") ,
                                  exp.cells=5)

sc.dat.filtered.pc <- BayesPrism::select.gene.type(sc.dat.filtered, 
                                                   gene.type = "protein_coding")# Subset protein coding genes

diff.exp.stat <- BayesPrism::get.exp.stat(
  sc.dat=sc.dat.filtered.pc[,colSums(as.matrix(sc.dat.filtered.pc)>0)>3],       # filter genes to reduce memory use
  cell.type.labels=sc.meta.subset,
  cell.state.labels=sc.meta.state,
  pseudo.count=0.1,                                                             # a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
  cell.count.cutoff=50,                                                         # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
  n.cores=5                                                                     # number of threads
)

sc.dat.filtered.pc.sig <- BayesPrism::select.marker(sc.dat=sc.dat.filtered.pc,
                                         stat=diff.exp.stat,
                                         pval.max=0.01,
                                         lfc.min=0.1)

sc.dat.filtered.pc.sig <- MatrixExtra::as.csc.matrix(sc.dat.filtered.pc.sig, 
                                                 binary = FALSE, logical = FALSE, sort = FALSE)
myPrism <- BayesPrism::new.prism(
  reference=sc.dat.filtered.pc.sig, 
  mixture=bulk.data,
  input.type="count.matrix", 
  cell.type.labels = sc.meta.subset, 
  cell.state.labels = sc.meta.state,
  key="Glioma",
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

bp.res <- BayesPrism::run.prism(prism = myPrism, n.cores=10)
saveRDS(bp.res, "./tmp/bp.res.rds")

theta <- BayesPrism::get.fraction (bp=bp.res,
                       which.theta="first",
                       state.or.type="state")
write.csv(theta, "./result/thera.first.state.csv", row.names = TRUE, quote = TRUE)
