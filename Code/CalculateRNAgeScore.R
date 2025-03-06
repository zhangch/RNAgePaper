number = 100

source("CommonFunc.R")
###############
## markers
###############
loadingMarkers()

###############
## plugin dataset
###############
scoreCard <- NULL
percentageCard <-NULL

Total_pF.9vs9 <- read.csv("../Data/TPM/pF_TPM.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)
Total_pF.9vs9 <- cbind(`Description`=rownames(Total_pF.9vs9), Total_pF.9vs9)
result <- getAgingScores(Total_pF.9vs9,9)
scoreCard <- cbind(scoreCard, `Young vs Old\n(Total pF)`=result[,1])
percentageCard <- cbind(percentageCard, `Young vs Old\n(Total pF)`=result[,2])





rownames(scoreCard) <- c("Pan markers","pF markers","Neuronal markers", "FC markers", "SN markers")
rownames(percentageCard) <- rownames(scoreCard)

write.csv(scoreCard.edgington, file="FiugreS1_Score.csv")
write.csv(percentageCard.edgington, file="FiugreS1_Percentage.csv")







