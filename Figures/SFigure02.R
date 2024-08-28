#####R version: 3.2.1#####
library(scales)
library(ggplot2)

#############################################
## S. Fig 2B
#############################################
pdf("FigureS2B.pdf", width=6, height=4, useDingbats = FALSE)
  pF_4vs4 <- read.csv("../Data/DEG/DEG_pF_4vs4.csv", header=T, stringsAsFactors=FALSE,row.names=1)
  pF_4vs4 <- pF_4vs4[!is.na(pF_4vs4$pvalue),]
  colorVec <- rep("Black", nrow(pF_4vs4))
  colorVec[pF_4vs4$log2FoldChange>0&pF_4vs4$padj<=0.1] <- "#1E90FF"
  colorVec[pF_4vs4$log2FoldChange<0&pF_4vs4$padj<=0.1] <- "#DC143C"
  
  ggplot(data=pF_4vs4, aes(x=-log10(padj), y=log2FoldChange, color=colorVec)) + 
    geom_point(alpha = 0.8) +
    xlim(0,3) +
    ylim(-5,5) +
    ggtitle("Young vs Old pFIB") + 
    theme_bw() + 
    theme(plot.title = element_text(size=25, face="bold", hjust = 0.5),
          axis.text=element_text(size=16),
          axis.title.x=element_text(size=20,face="bold"),
          axis.title.y=element_text(size=20,face="bold")) +
    xlab(expression(paste(-Log[10] ,'-Log10(FDR)',sep=""))) +
    ylab(expression(paste(Log[2] ,'Log2FC(Young/Old)',sep=""))) +
    guides(alpha = "none", color = "none") +
    scale_color_manual(values=c("#1E90FF", "#DC143C", "black"))+
    annotate(geom="text", x=2.5, y=-4.8, label="466 genes", color="#DC143C", size=8) +
    annotate(geom="text", x=2.5, y=4.8, label="402 genes", color="#1E90FF", size=8)
  
  iPS_3vs3 <- read.csv("../Data/DEG/DEG_iPS_3vs3.csv", header=T, stringsAsFactors=FALSE,row.names=1)
  colorVec <- rep("Black", nrow(iPS_3vs3))
  colorVec[iPS_3vs3$log2FoldChange>0&iPS_3vs3$padj<=0.1] <- "#1E90FF"
  colorVec[iPS_3vs3$log2FoldChange<0&iPS_3vs3$padj<=0.1] <- "#DC143C"
  
  ggplot(data=iPS_3vs3, aes(x=-log10(padj), y=log2FoldChange, color=colorVec)) + 
    geom_point(alpha = 0.8) +
    xlim(0,3) +
    ylim(-5,5) +
    ggtitle("Young vs Old iPSC") + 
    theme_bw() + 
    theme(plot.title = element_text(size=25, face="bold", hjust = 0.5),
          axis.text=element_text(size=16),
          axis.title.x=element_text(size=20,face="bold"),
          axis.title.y=element_text(size=20,face="bold")) +
    xlab(expression(paste(-Log[10] ,'-Log10(FDR)',sep=""))) +
    ylab(expression(paste(Log[2] ,'Log2FC(Young/Old)',sep=""))) +
    guides(alpha = "none", color = "none") +
    scale_color_manual(values=c("#1E90FF", "#DC143C", "black"))+
    annotate(geom="text", x=2.5, y=-4.8, label="1 genes", color="#DC143C", size=8) +
    annotate(geom="text", x=2.5, y=4.8, label="11 genes", color="#1E90FF", size=8)
  
  iF_3vs3 <- read.csv("../Data/DEG/DEG_iF_3vs3.csv", header=T, stringsAsFactors=FALSE,row.names=1)
  colorVec <- rep("Black", nrow(iF_3vs3))
  colorVec[iF_3vs3$log2FoldChange>0&iF_3vs3$padj<=0.1] <- "#1E90FF"
  colorVec[iF_3vs3$log2FoldChange<0&iF_3vs3$padj<=0.1] <- "#DC143C"
  
  ggplot(data=iF_3vs3, aes(x=-log10(padj), y=log2FoldChange, color=colorVec)) + 
    geom_point(alpha = 0.8) +
    xlim(0,3) +
    ylim(-5,5) +
    ggtitle("Young vs Old iPSC-FIB") + 
    theme_bw() + 
    theme(plot.title = element_text(size=25, face="bold", hjust = 0.5),
          axis.text=element_text(size=16),
          axis.title.x=element_text(size=20,face="bold"),
          axis.title.y=element_text(size=20,face="bold")) +
    xlab(expression(paste(-Log[10] ,'-Log10(FDR)',sep=""))) +
    ylab(expression(paste(Log[2] ,'Log2FC(Young/Old)',sep=""))) +
    guides(alpha = "none", color = "none") +
    scale_color_manual(values=c("#1E90FF", "#DC143C", "black"))+
    annotate(geom="text", x=2.5, y=-4.8, label="23 genes", color="#DC143C", size=8) +
    annotate(geom="text", x=2.5, y=4.8, label="53 genes", color="#1E90FF", size=8)
dev.off()
