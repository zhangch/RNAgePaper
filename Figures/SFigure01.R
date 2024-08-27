#####R version: 3.6.2#####
library(ggplot2)
library(vegan)

#############################################
## Fibroblast - Figure S1 B.1 & C.1
#############################################
count.matrix <- read.csv("../Data/RawCounts/pF_Gencode_count.csv", header=T, stringsAsFactors=FALSE, row.names=1, check.names=FALSE)
sample.table1 <- read.csv(".../Data//RawCounts/pF_metadata.csv", header=T, stringsAsFactors=FALSE,row.names=1)
count.matrix <- count.matrix[rowSums(count.matrix)>=5,]
res.FC <- read.csv("../Data/DEG/DEG_Total_pF_YvsO.csv", header=T, stringsAsFactors=FALSE,row.names=1)

pdf("FigureS1B_PF.pdf", width=4, height=4, useDingbats = FALSE)
  pc <- prcomp(t(count.matrix))
  percentVar <- pc$sdev^2 / sum( pc$sdev^2)
  dataset <- cbind(Sample=rownames(sample.table), t(count.matrix), sample.table)
  
  distance.matrix <- vegdist(t(count.matrix), method = "bray")
  mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE)
  mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)
  mds.values <- mds.stuff$points
  mds.data <- data.frame(Sample=rownames(mds.values),
                         X=mds.values[,1],
                         Y=mds.values[,2])
  dataset <- cbind(mds.data, sample.table)
  ggplot(data=dataset, aes(x=X, y=Y, color=group)) +
    geom_point(shape=19,size=6,alpha = 0.75) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text=element_text(size=16),
          axis.title=element_text(size=20,face="bold")) +
    xlab(paste("PC1 - ", mds.var.per[1], "%", sep="")) +
    ylab(paste("PC2 - ", mds.var.per[2], "%", sep="")) +
    guides(alpha = "none", color = "none") +
    scale_color_manual(values=c("#DC143C", "#1E90FF"))
dev.off()

colorVec <- rep("Black", nrow(res.FC))
colorVec[res.FC$log2FoldChange>0&res.FC$padj<=0.1] <- "#1E90FF"
colorVec[res.FC$log2FoldChange<0&res.FC$padj<=0.1] <- "#DC143C"
pdf("Figure1B_PF_Volcano_v3.pdf", width=4, height=4, useDingbats = FALSE)
  ggplot(data=res.FC, aes(x=-log10(padj), y=log2FoldChange, color=colorVec)) + 
    geom_point(alpha = 0.8) +
    xlim(0,12) +
    ylim(-6,6) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text=element_text(size=16),
          axis.title=element_text(size=20,face="bold")) +
    xlab(expression(paste(-Log[10] ,'-log10(FDR)',sep=""))) +
    ylab(expression(paste(Log[2] ,'FC(Young/Old)',sep=""))) +
    guides(alpha = "none", color = "none") +
    scale_color_manual(values=c("#1E90FF", "#DC143C", "black"))+
    annotate(geom="text", x=3.4, y=5.8, label="748 genes", color="#1E90FF", size=8) +
    annotate(geom="text", x=3.4, y=-5.8, label="554 genes", color="#DC143C", size=8)
dev.off()
