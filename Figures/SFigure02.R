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

#############################################
## Figure S2D
#############################################
scoreCard <- read.csv("../Data/RNAgeScore/FiugreS2_3_Score.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)
perCard <- read.csv("../Data/RNAgeScore/FiugreS2_3_Percentage.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)
scoreCard <- scoreCard[,-3]
perCard <- perCard[,-3]
dist<-as.data.frame(scoreCard)
dist <- cbind(`id`=rownames(dist), dist)
df1 = melt(dist)
dist<-as.data.frame(perCard)
dist <- cbind(`id`=rownames(dist), dist)
df2 = melt(dist)
df = cbind(df1, df2$value)
colnames(df) <- c('Markers', 'Group', 'AgingScore', 'Percentage')

pdf("FigureS02C.pdf", width = 6, height = 2, useDingbats = FALSE)
  ggplot(df, aes(Markers, Group)) +
    geom_point(aes(size=Percentage, fill=AgingScore), colour="Grey", stroke = 1, shape=21) +
    scale_size(range = c(0, 10)) +
    theme_bw() + scale_y_discrete(limits=rev(colnames(scoreCard))) + scale_x_discrete(limits=rownames(scoreCard)) +
    theme(panel.grid.major = element_blank(), panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_gradient2(low="#00B0F0", high="#FF0000", mid="white")
dev.off()
