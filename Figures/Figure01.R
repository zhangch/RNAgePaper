library(ggplot2)
library(reshape2)

orignial.avg <- read.csv("Training_AgingScore_13.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)
orignial.per <- read.csv("Training_Percentage_13.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)

Fleischer.avg <- read.csv("Fleischer_GSE113957_AgingScore_5.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)
Fleischer.per <- read.csv("Fleischer_GSE113957_Percentage_5.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)

GTEx.avg <- read.csv("GTEx_AgingScore_5.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)
GTEx.per <- read.csv("GTEx_Percentage_5.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)

Holland.avg <- read.csv("Holland_GSE36192_AgingScore_5.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)
Holland.per <- read.csv("Holland_GSE36192_Percentage_5.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)

TabulaMuris.avg <- read.csv("TabulaMuris_Brain_AgingScore_5.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)
TabulaMuris.per <- read.csv("TabulaMuris_Brain_Percentage_5.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)

Mertens.avg <- read.csv("Mertens_E-MTAB-10352_AgingScore_5.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)
Mertens.per <- read.csv("Mertens_E-MTAB-10352_Percentage_5.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)

OurHGPS.avg <- read.csv("Our_HGPS_AgingScore_5.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)
OurHGPS.per <- read.csv("Our_HGPS_Percentage_5.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)

Fathi.avg <- read.csv("Fathi_SLO_AgingScore_5.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)
Fathi.per <- read.csv("Fathi_SLO_Percentage_5.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)

Miller.avg <- read.csv("Miller_GSE52431_AgingScore_5.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)
Miller.per <- read.csv("Miller_GSE52431_Percentage_5.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)

Riessland.avg <- read.csv("Riessland_SATB1_AgingScore_5.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)
Riessland.per <- read.csv("Riessland_SATB1_Percentage_5.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)

scoreCard <- cbind(orignial.avg[1:5,c(1:3,6:8)], 
                   Fleischer.avg[,1], Holland.avg[,1], GTEx.avg[,c(4,2,5)], TabulaMuris.avg[,1], 
                   GTEx.avg[,6], Holland.avg[,4],
                   orignial.avg[1:5,c(9,10,5)],
                   Fleischer.avg[,4], OurHGPS.avg[,2], 
                   Mertens.avg[,1],
                   Miller.avg[,1], Riessland.avg[,1], Fathi.avg)

perCard <- cbind(orignial.per[1:5,c(1:3,6:8)], 
                 Fleischer.per[,1], Holland.per[,1], GTEx.per[,c(4,2,5)], TabulaMuris.per[,1], 
                 GTEx.per[,6], Holland.per[,4],
                 orignial.per[1:5,c(9,10,5)],
                 Fleischer.per[,4], OurHGPS.per[,2], 
                 Mertens.per[,1],
                 Miller.per[,1], Riessland.per[,1], Fathi.per)


colnames(scoreCard) <- c("Young vs Old (Total pF)", "Young vs Old (PolyA pF)", "Young vs Old (Gage pF)", "Young vs Old (Total FC)", "Young vs Old (Gage FC)", "Young vs Old (Total SN)",
                         "Young vs Old (Fleischer pF)", "Young vs Old (Holland FC)", "Young vs Old (GTEx FC)", "Young vs Old (GTEx Cortex)", "Young vs Old (GTEx SN)", "Young vs Old (Mouse Brain)",
                         "Young vs Old (GTEx Stomach)", "Young vs Old (Holland Cerebellum)",
                         "iF vs pF (Old Donors)", "iF vs pF (Young Donors)", "Young vs Old (iF)",
                         "Young vs HGPS (Fleischer pF)", "Young vs HGPS (PolyA pF)",
                         "iPSCiN vs iN (Mertens)",
                         "WT vs Progerin (Miller)", "WT vs SATB1_KO (Riessland)", "WT vs SLO (Fathi)"
)
colnames(perCard) <- colnames(scoreCard)

rownames(scoreCard) <- c("Pan markers", "pF markers", "Neuronal markers", "FC markers", "SN markers")
rownames(perCard) <- rownames(scoreCard)

dist<-as.data.frame(scoreCard)
dist <- cbind(`id`=rownames(dist), dist)
df1 = melt(dist)
dist<-as.data.frame(perCard)
dist <- cbind(`id`=rownames(dist), dist)
df2 = melt(dist)
df = cbind(df1, df2$value)
colnames(df) <- c('Markers', 'Group', 'AgingScore', 'Percentage')

pdf("Bubble_Plot_All_by_figure.pdf", width = 6, height = 10, useDingbats = FALSE)
ggplot(df, aes(Markers, Group)) +
  geom_point(aes(size=Percentage, fill=AgingScore), colour="Grey", stroke = 1, shape=21) +
  #geom_text(aes(label=Percentage),hjust=0.5,vjust=0.5) +
  scale_size(range = c(0, 10)) +
  theme_bw() + scale_y_discrete(limits=rev(colnames(scoreCard))) + scale_x_discrete(limits=rownames(scoreCard)) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))+
  #theme(panel.grid.major = element_blank(), panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), 
  #      axis.text.y = element_text(colour = rev(c(rep("blue",3), rep("orange",2), rep("purple",1))))) + 
  #scale_color_gradientn(colours = rev(brewer.pal(11, "RdYlBu")))
  scale_fill_gradient2(low="#00B0F0", high="#FF0000", mid="white")
#scale_fill_gradient2(low="#00B0F0", high="#FF0000", mid="white",#colors in the scale
#midpoint=mean(rng),    #same midpoint for plots (mean of the range)
#breaks=seq(-100,100,4), #breaks in the scale bar
#limits=c(floor(min(df$AgingScore)), ceiling(max(df$AgingScore)))) #same limits for plots
dev.off()

