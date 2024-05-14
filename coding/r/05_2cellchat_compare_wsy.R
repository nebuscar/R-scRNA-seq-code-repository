library(CellChat)
setwd("F:/Projects/Forstudents/NASH")

# Load dataset
cc.list <- readRDS("tmp/cellchat.list.rds")

object.list <- list(group1 = cc.list[["Chow"]], group2 = cc.list[["15weeks"]])
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

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
for (i in 1:length(cc.list)) {
  groupSize <- as.numeric(table(cc.list[[i]]@idents))
  netVisual_circle(cc.list[[i]]@net$count, vertex.weight = groupSize,
                   weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(cc.list)[i]))
}

par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:length(cc.list)) {
  groupSize <- as.numeric(table(cc.list[[i]]@idents))
  netVisual_circle(cc.list[[i]]@net$weight, vertex.weight = groupSize,
                   weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(cc.list)[i]))
}

# Signaling: receive
for (time in names(cc.list)){
  mat <- cc.list[[time]]@net$count
  groupSize <- as.numeric(table(cc.list[[time]]@idents))
  pdf(paste("./result/cellchat/cellchat", time, "receive.pdf", sep = "_"), 
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
  pdf(paste("./result/cellchat/cellchat", time, "sent.pdf", sep = "_"), 
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
}


