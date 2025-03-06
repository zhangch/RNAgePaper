library(ggplot2)
library(reshape2)

scoreCard <- read.csv("../Data/RNAgeScore/Fiugre1_Score.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)
perCard <- read.csv("../Data/RNAgeScore/Fiugre1_Percentage.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)

dist<-as.data.frame(scoreCard)
dist <- cbind(`id`=rownames(dist), dist)
df1 = melt(dist)
dist<-as.data.frame(perCard)
dist <- cbind(`id`=rownames(dist), dist)
df2 = melt(dist)
df = cbind(df1, df2$value)
colnames(df) <- c('Markers', 'Group', 'AgingScore', 'Percentage')

pdf("Figure01.pdf", width = 6, height = 10, useDingbats = FALSE)
  ggplot(df, aes(Markers, Group)) +
    geom_point(aes(size=Percentage, fill=AgingScore), colour="Grey", stroke = 1, shape=21) +
    scale_size(range = c(0, 10)) +
    theme_bw() + scale_y_discrete(limits=rev(colnames(scoreCard))) + scale_x_discrete(limits=rownames(scoreCard)) +
    theme(panel.grid.major = element_blank(), panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_gradient2(low="#00B0F0", high="#FF0000", mid="white")
dev.off()
