library(CellChat)

# Load dataset
# cellchat <- readRDS("tmp/cellchat_scissor_tcells.rds")
cellchat <- readRDS("tmp/cellchat_annotation.rds")

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

weight.max <- getMaxWeight(cc.list, attribute = c("idents","count"))
par(mfrow = c(2,2), xpd=TRUE)


# Signaling: receive
mat <- cellchat@net$count
groupSize <- as.numeric(table(cellchat@idents))
pdf(paste("tmp/results/cell", "receive.pdf", sep = "_"),
    width = 20, height = 12, useDingbats = F)
par(mfrow = c(3,5), xpd = TRUE)
for(i in 1:ncol(mat)){
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i,] <- mat[,i]
  netVisual_circle(mat2, vertex.weight = groupSize,
                   weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
}
dev.off()
  
# Signaling: sent
pdf(paste("./result/", "sent.pdf", sep = "_"), 
    width = 20, height = 12, useDingbats = F)
par(mfrow = c(3,5), xpd = TRUE)
for(i in 1:nrow(mat)){
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2, vertex.weight = groupSize,
                   weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
}
dev.off()


