#####R version: 3.2.1#####
require(DESeq2)
library(ggplot2)

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

par(cex.lab=1.5, cex.axis=1.2, cex=1.5)
#############################################
## Fibroblast - Figure3B
#############################################
count.matrix <- read.csv("../../1_Paper_Data/RawCount/13814_B_Gencode_count.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)
sample.table <- read.csv("../../1_Paper_Data//RawCount/13814_B_Metadata.csv", header=T, stringsAsFactors=FALSE,row.names=1)
sample.table <- sample.table[sample.table$Cell_line=="Neonatal_Fib",]
count.matrix <- count.matrix[,match(rownames(sample.table), colnames(count.matrix))]
count.matrix <- count.matrix[rowSums(count.matrix)>10*ncol(count.matrix),]

dse <- DESeqDataSetFromMatrix(count.matrix, sample.table, ~ Condition)
design(dse) <- formula(~ Condition)
dse <- DESeq(dse)
dds <- estimateSizeFactors(dse)
normalized_counts <- counts(dds, normalized=TRUE)
cdsFullBlind <- estimateDispersions(dds)
vsdFull <- varianceStabilizingTransformation(cdsFullBlind)
pcaDatayoungpf <- plotPCA(vsdFull, intgroup=c("Condition"), ntop = 5000)
mds.data <- data.frame(Sample=pcaDatayoungpf$data$name,
                       X=pcaDatayoungpf$data$PC1,
                       Y=pcaDatayoungpf$data$PC2)
dataset <- cbind(mds.data, sample.table)
mds.var.per <- pcaDatayoungpf$data$PC1

pdf("Figure3B.pdf", width=3, height=3, useDingbats = FALSE)
ggplot(data=dataset, aes(x=X, y=Y, color=Condition)) +
  geom_point(shape=19, size=4, show.legend = TRUE) + 
  scale_color_manual(values = c("#ffff33", "#f781bf", "Black", "#a65628")) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     legend.title = element_blank(),
                     legend.spacing.y = unit(0, "mm"), 
                     panel.border = element_rect(colour = "black", fill=NA),
                     aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = "black")) +
  xlab(pcaDatayoungpf$labels$x) +
  ylab(pcaDatayoungpf$labels$y) +
  ggtitle("Neonatal Fibroblast")
dev.off()


#############################################
## neuron
#############################################
count.matrix <- read.csv("../../1_Paper_Data/RawCount/13963_B_Gencode_count.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)
sample.table <- read.csv("../../1_Paper_Data//RawCount/13963_B_Metadata.csv", header=T, stringsAsFactors=FALSE,row.names=1)
count.matrix <- count.matrix[,match(rownames(sample.table), colnames(count.matrix))]
count.matrix <- count.matrix[rowSums(count.matrix)>10*ncol(count.matrix),]

dse <- DESeqDataSetFromMatrix(count.matrix, sample.table, ~ Condition)
design(dse) <- formula(~ Condition)
dse <- DESeq(dse)
dds <- estimateSizeFactors(dse)
normalized_counts <- counts(dds, normalized=TRUE)
cdsFullBlind <- estimateDispersions(dds)
vsdFull <- varianceStabilizingTransformation(cdsFullBlind)
pcaDatayoungpf <- plotPCA(vsdFull, intgroup=c("Condition"), ntop = 5000)
mds.data <- data.frame(Sample=pcaDatayoungpf$data$name,
                       X=pcaDatayoungpf$data$PC1,
                       Y=pcaDatayoungpf$data$PC2)
dataset <- cbind(mds.data, sample.table)
mds.var.per <- pcaDatayoungpf$data$PC1

mds.data <- data.frame(Sample=pcaDataAllneuron$data$name,
                       X=pcaDataAllneuron$data$PC1,
                       Y=pcaDataAllneuron$data$PC2)
dataset <- cbind(mds.data, sample.table)

pdf("Figure3D.pdf", width=3, height=3, useDingbats = FALSE)
  ggplot(data=dataset, aes(x=X, y=Y, color=Condition)) +
    geom_point(shape=19, size=4, show.legend = TRUE) + #geom_text(aes(label=Sample), hjust=0.5, vjust=1.7) +
    scale_color_manual(values = c("#e41a1c", "Black", "#377eb8", "#4daf4a")) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                       legend.title = element_blank(),
                       legend.spacing.y = unit(0, "mm"), 
                       panel.border = element_rect(colour = "black", fill=NA),
                       aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
                       legend.background = element_blank(),
                       legend.box.background = element_rect(colour = "black")) +
    xlab(pcaDataAllneuron$labels$x) +
    ylab(pcaDataAllneuron$labels$y) +
    ggtitle("Neuron")
dev.off()



#############################################
## Fibroblast - Figure3C & 3E
#############################################
library(RColorBrewer)
library(tidyverse)

make_gradient <- function(leftPanel = 500, rightPanel = 500) {
  cols <- c(rev(colorRampPalette(brewer.pal(9, "Blues"))(300)[1:leftPanel]),
            "#FFFFFF",
            colorRampPalette(brewer.pal(9, "Reds"))(300)[1:rightPanel]
  )
  mat <- t(rep(cols, each=500))
  grid::rasterGrob(
    image = mat,
    width = unit(1, "npc"),
    height = unit(1, "npc"),
    interpolate = TRUE
  )
}

par(cex.lab=1.5, cex.axis=1.2, cex=1.5)
#################### Figure 3C
KO.pf.as <- t(read.csv("../Data/RNAgeScore/13814_B_AgingScore_5.csv", header=F, stringsAsFactors=FALSE, row.names = 1))
KO.pf.per <- t(read.csv("../Data/RNAgeScore/13814_B_Percentage_5.csv", header=F, stringsAsFactors=FALSE, row.names = 1))
colnames(KO.pf.as) <- c("Group", "Pan markers", "pF markers", "Neuronal markers", "FC markers", "SN markers")
colnames(KO.pf.per) <- c("Group", "Pan markers", "pF markers", "Neuronal markers", "FC markers", "SN markers")
KO.pf.as <- KO.pf.as[c(4:6), ]
KO.pf.per <- KO.pf.per[c(4:6), ]
KO.pf.as[,1] <- c("DMSO vs 5-iodotubercidin", "DMSO vs Alvocidib", "DMSO vs Mitoxantrone")
KO.pf.per[,1] <- c("DMSO vs 5-iodotubercidin", "DMSO vs Alvocidib", "DMSO vs Mitoxantrone")

df <- full_join(
  as.data.frame(KO.pf.as) %>% 
    pivot_longer(cols = 2:6, names_to ="Cell type", values_to = "AgingScore"),
  as.data.frame(KO.pf.per) %>% 
    pivot_longer(cols = 2:6, names_to ="Cell type", values_to = "Percentage"),
  by = c('Group', 'Cell type'))

df$AgingScore <- as.numeric(df$AgingScore)
df$Percentage <- as.numeric(df$Percentage)

g <- make_gradient(42, 134)

Figure3C <-   
  ggplot(df, aes(x = `Cell type`, y = AgingScore))+
  annotation_custom(
    grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) + 
  geom_point(aes(color = Group, size=Percentage), alpha = 0.7) + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        legend.text=element_text(size=15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15))+
  scale_size(range=c(0,10)) +
  scale_x_discrete(limits=rev(colnames(KO.pf.per[,-1]))) +
  scale_colour_manual(values = c("#ffff33", "#f781bf", "#a65628")) +
  guides(color = guide_legend(override.aes = list(size=4))) +
  scale_size(range=c(0,12),breaks=c(20,40,60,80), labels=c("20","40","60","80"), guide="legend", limits = c(0, 100)) +
  labs(
    title = "",
    x    = "",
    y    = "RNAge Score",
    fill = ""
  ) +
  coord_flip()

#################### Figure 3E
KO.pf.as <- t(read.csv("../Data/RNAgeScore/13963_B_AgingScore_5.csv", header=F, stringsAsFactors=FALSE, row.names = 1))
KO.pf.per <- t(read.csv("../Data/RNAgeScore/13963_B_Percentage_5.csv", header=F, stringsAsFactors=FALSE, row.names = 1))
colnames(KO.pf.as) <- c("Group", "Pan markers", "pF markers", "Neuronal markers", "FC markers", "SN markers")
colnames(KO.pf.per) <- c("Group", "Pan markers", "pF markers", "Neuronal markers", "FC markers", "SN markers")
KO.pf.as[,1] <- c("DMSO vs AGK-2", "DMSO vs Fludarabine", "DMSO vs Vinorelbine")
KO.pf.per[,1] <- c("DMSO vs AGK-2", "DMSO vs Fludarabine", "DMSO vs Vinorelbine")

df <- full_join(
  as.data.frame(KO.pf.as) %>% 
    pivot_longer(cols = 2:6, names_to ="Cell type", values_to = "AgingScore"),
  as.data.frame(KO.pf.per) %>% 
    pivot_longer(cols = 2:6, names_to ="Cell type", values_to = "Percentage"),
  by = c('Group', 'Cell type'))

df$AgingScore <- as.numeric(df$AgingScore)
df$Percentage <- as.numeric(df$Percentage)

g <- make_gradient(64, 129)

Figure3E <-  
  ggplot(df, aes(x = `Cell type`, y = AgingScore))+
  annotation_custom(
    grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) + 
  geom_point(aes(color = Group, size=Percentage), alpha = 0.7) + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        legend.text=element_text(size=15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15))+
  scale_x_discrete(limits=rev(colnames(KO.pf.per[,-1]))) +
  scale_colour_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) +
  guides(color = guide_legend(override.aes = list(size=4))) +
  scale_size(range=c(0,12),breaks=c(20,40,60,80), labels=c("20","40","60","80"), guide="legend", limits = c(0, 100)) +
  labs(
    title = "",
    x    = "",
    y    = "RNAge Score",
    fill = ""
  ) +
  coord_flip()

library(grid)
library(gridExtra)
library(gtable) 
# Get the gtables
gA <- ggplotGrob(Figure3C)
gB <- ggplotGrob(Figure3E)

# Set the widths
gB$widths <- gA$widths

pdf("Figure3CE.pdf", width = 11, height = 9, useDingbats = FALSE)
  grid.newpage()
  grid.arrange(gA, gB, nrow = 2)
dev.off()

