library(stringi)
library(ggplot2)
library(dplyr)

##########KEGG##########
downgokegg<-read.delim("./results/competition/KEGG.txt")
enrich<-downgokegg
enrich_signif=enrich[which(enrich$PValue<0.05),]
enrich_signif=enrich_signif[,c(1:3,5)]
head(enrich_signif)
enrich_signif=data.frame(enrich_signif)
KEGG=enrich_signif
KEGG$Term<-stri_sub(KEGG$Term,10,100)
kegg_plot <- ggplot(KEGG, aes(x = Count, y = Term)) +
  geom_point(aes(color = PValue, size = Count)) +
  scale_color_gradient(low = 'slateblue4', high = 'firebrick3') +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
ggsave(filename = "./results/competition/KEGG_enrichment_plot_res1_IntegrateData_cca.png", plot = kegg_plot, width = 10, height = 6)

##########GO##########
# GO_CC
GO_CC<-read.delim("./results/wsy/GO_CC_cluster06.txt")
GO_CC_signif=GO_CC[which(GO_CC$PValue<0.05),]
GO_CC_signif=GO_CC[,c(1:3,5)]
head(GO_CC_signif)
GO_CC_signif=data.frame(GO_CC_signif)
GO_CC_signif$Term<-stri_sub(GO_CC_signif$Term,12,100)

# GO_BP
GO_BP<-read.delim("./results/wsy/GO_BP_cluster06.txt")
GO_BP_signif=GO_BP[which(GO_BP$PValue<0.05),]
GO_BP_signif=GO_BP_signif[,c(1:3,5)]
head(GO_BP_signif)
GO_BP_signif=data.frame(GO_BP_signif)
GO_BP_signif$Term<-stri_sub(GO_BP_signif$Term,12,100)

# GO_MF
GO_MF<-read.delim("./results/wsy/GO_MF_cluster06.txt")
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

# 生成条形图
go_plot <- ggplot(data = go, aes(x = GO_term_order, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = COLS) +
  theme_bw() +
  xlab("Terms") +
  ylab("Gene counts") +
  labs() +
  coord_flip() +
  theme(
    axis.text.y = element_text(face = "bold", color = "black", angle = 0, vjust = 0.5, hjust = 1)
  )
ggsave(filename = "./results/wsy/go_enrichment_plot_cluster06.png", plot = go_plot, width = 15, height = 15)

