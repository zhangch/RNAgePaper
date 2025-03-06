loadingMarkers <- function(){
  L1.markers <<- read.csv("../02_Block/L1_common_all_pvals.csv", header=T, stringsAsFactors=FALSE, row.names=1)
  L2.markers <<- read.csv("../02_Block/L2_common_all_pF_pvals.csv", header=T, stringsAsFactors=FALSE, row.names=1)
  L6.markers <<- read.csv("../02_Block/L6_pan_Brain_pvals.csv", header=T, stringsAsFactors=FALSE, row.names=1)
  L7.markers <<- read.csv("../02_Block/L7_FC_pvals.csv", header=T, stringsAsFactors=FALSE, row.names=1)
  L8.markers <<- read.csv("../02_Block/L8_SN_pvals.csv", header=T, stringsAsFactors=FALSE,row.names=1, check.names=FALSE)
  colnames(L8.markers)[2] <- "AvgLogFC"
  
  L1.markers <<- L1.markers[order(L1.markers$edgington), ]
  L2.markers <<- L2.markers[order(L2.markers$edgington), ]
  L6.markers <<- L6.markers[order(L6.markers$edgington), ]
  L7.markers <<- L7.markers[order(L7.markers$edgington), ]
  L8.markers <<- L8.markers[order(L8.markers$padj), ]
}

getAgingScores <- function(data, num1){
  #data = GTEx_Stomach.5vs5
  #num1 = 9
  #tstatistic(L1.markers, data, num1)
  return(rbind(tstatistic(L1.markers, data, num1),
               tstatistic(L2.markers, data, num1),
               tstatistic(L6.markers, data, num1),
               tstatistic(L7.markers, data, num1),
               tstatistic(L8.markers, data, num1)))
}

tstatistic <- function(markers, data, num1){
  #markers = L8.markers
  #data = GTEx_FC.9vs9
  #num1 = 9
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
    write.csv(dist, "L1_markers.csv")
    dist <- dist*sign(markers$AvgLogFC)
    dist[is.na(dist)] <- 0
    
    results <- sign(dist)
    results[results<0]<-0
    sum(results)
    return(c(sum(dist)/n, sum(results)))
  }
}
