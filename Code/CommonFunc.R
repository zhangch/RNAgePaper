loadingMarkers <- function(){
  L1.markers <<- read.csv("../Data/AgingMarkers/L1_Pan_Markers.csv", header=T, stringsAsFactors=FALSE, row.names=1)
  L2.markers <<- read.csv("../Data/AgingMarkers/L2_pF_Markers.csv", header=T, stringsAsFactors=FALSE, row.names=1)
  L6.markers <<- read.csv("../Data/AgingMarkers/L6_Neuronal_Markers.csv", header=T, stringsAsFactors=FALSE, row.names=1)
  L7.markers <<- read.csv("../Data/AgingMarkers/L7_FC_Markers.csv", header=T, stringsAsFactors=FALSE, row.names=1)
  L8.markers <<- read.csv("../Data/AgingMarkers/L8_SN_Markers.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)
  colnames(L8.markers)[2] <- "AvgLogFC"
  
  L1.markers <<- L1.markers[order(L1.markers$CombinedP), ]
  L2.markers <<- L2.markers[order(L2.markers$CombinedP), ]
  L6.markers <<- L6.markers[order(L6.markers$CombinedP), ]
  L7.markers <<- L7.markers[order(L7.markers$CombinedP), ]
  L8.markers <<- L8.markers[order(L8.markers$pvalue), ]
}

getAgingScores <- function(data, num1){
  return(rbind(tstatistic(L1.markers, data, num1),
               tstatistic(L2.markers, data, num1),
               tstatistic(L6.markers, data, num1),
               tstatistic(L7.markers, data, num1),
               tstatistic(L8.markers, data, num1)))
}

tstatistic <- function(markers, data, num1){
  if(nrow(markers)==0) {
    return(0)
  } else {
    markers <- markers[rownames(markers) %in% data$Description,]
    data.markers <- data[match(rownames(markers), data$Description),]
    markers <- markers[match(data.markers$Description, rownames(markers)),]
    if(nrow(markers)>=100) {
      n=100
    } else {
      n=nrow(markers)
    }
    data.markers <- data.markers[1:n,-1]
    markers <- markers[1:n,]
    
    dist <- apply(data.markers, 1, function(x) t.test(x[1:num1], x[(1+num1):ncol(data.markers)])$statistic)
    dist <- dist*sign(markers$AvgLogFC)
    dist[is.na(dist)] <- 0
    
    results <- sign(dist)
    results[results<0]<-0
    sum(results)
    return(c(sum(dist)/n, sum(results)))
  }
}
