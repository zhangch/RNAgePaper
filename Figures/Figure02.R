#####R version: 3.2.1#####
library(tidyverse)
library(ggdist)
library(tidyquant)

tag = "mean"

markerNames = c("Pan markers", "pF markers", "Neuronal markers", "FC markers")

pdf(paste0("Figure02.pdf"), width = 12, height = 4, useDingbats = FALSE)

detail.1HAE <- read.csv("../Data/L1000/1HAE_allData.csv", header=T, stringsAsFactors=FALSE)
colnames(detail.1HAE)[4:7] <- markerNames
if(tag=="mean") {
  result<- detail.1HAE %>% 
    dplyr::select(2, 4:7) %>% 
    group_by(pert_id) %>% summarise(across(everything(), mean))
} else {
  result<- detail.1HAE %>% 
    dplyr::select(2, 4:7) %>% 
    group_by(pert_id) %>% summarise(across(everything(), median))
}

selectedMarkers1 <- result[order(result[,-1]$`Pan markers`),]$`Pan markers`
selectedMarkers1 <- selectedMarkers1[c(1:4, 345:348)]
selectedMarkers2 <- result[order(result[,-1]$`pF markers`),]$`pF markers`
selectedMarkers2 <- selectedMarkers2[c(1:4, 345:348)]

result[,-1] %>% 
  pivot_longer(cols = 1:4, names_to ="Type", values_to = "Score")  %>% 
  ggplot(aes(x = factor(Type), y = Score, fill = factor(Type))) +
  # add half-violin from (ggdist} package
  ggdist::stat_halfeye(
    width = 0.6,
    ## custom bandwidth
    adjust = 0.5,
    ## move geom to the right
    justification = -.3,
    ## remove slab interval
    .width = 0,
    point_colour = NA
  ) +
  scale_x_discrete(limits=rev(markerNames)) +
  ylim(c(-0.8,1.75)) +
  geom_boxplot(
    width = .2,
    ## remove outliers
    outlier.color = NA,
    alpha = 0.5
  ) +
  # Adjust theme
  scale_fill_tq() +
  theme_tq(base_size = 20) +
  theme(legend.position="right",legend.text=element_text(size=15)) +
  scale_fill_discrete(limits = markerNames)+
  labs(
    title = "1HAE (pF)",
    x    = "",
    y    = "RNAge",
    fill = ""
  ) +
  geom_point(aes(x = 4, y = selectedMarkers1[1])) +
  geom_point(aes(x = 4, y = selectedMarkers1[2])) +
  geom_point(aes(x = 4, y = selectedMarkers1[3])) +
  geom_point(aes(x = 4, y = selectedMarkers1[4])) +
  geom_point(aes(x = 4, y = selectedMarkers1[5])) +
  geom_point(aes(x = 4, y = selectedMarkers1[6])) +
  geom_point(aes(x = 4, y = selectedMarkers1[7])) +
  geom_point(aes(x = 4, y = selectedMarkers1[8])) +
  geom_point(aes(x = 3, y = selectedMarkers2[1])) +
  geom_point(aes(x = 3, y = selectedMarkers2[2])) +
  geom_point(aes(x = 3, y = selectedMarkers2[3])) +
  geom_point(aes(x = 3, y = selectedMarkers2[4])) +
  geom_point(aes(x = 3, y = selectedMarkers2[5])) +
  geom_point(aes(x = 3, y = selectedMarkers2[6])) +
  geom_point(aes(x = 3, y = selectedMarkers2[7])) +
  geom_point(aes(x = 3, y = selectedMarkers2[8])) +
  coord_flip()



detail.NEU <- read.csv("../Data/L1000/NEU_allData.csv", header=T, stringsAsFactors=FALSE)
colnames(detail.NEU)[4:7] <- markerNames
if(tag=="mean") {
  result<- detail.NEU %>% 
    dplyr::select(2, 4:7) %>% 
    group_by(pert_id) %>% summarise(across(everything(), mean))
} else {
  result<- detail.NEU %>% 
    dplyr::select(2, 4:7) %>% 
    group_by(pert_id) %>% summarise(across(everything(), median))
}

selectedMarkers1 <- result[order(result[,-1]$`Pan markers`),]$`Pan markers`
selectedMarkers1 <- selectedMarkers1[c(1:4, 3965:3968)]
selectedMarkers2 <- result[order(result[,-1]$`Neuronal markers`),]$`Neuronal markers`
selectedMarkers2 <- selectedMarkers2[c(1:4, 3965:3968)]

result[,-1] %>% 
  pivot_longer(cols = 1:4, names_to ="Type", values_to = "Score")  %>% 
  ggplot(aes(x = factor(Type), y = Score, fill = factor(Type))) +
  # add half-violin from (ggdist} package
  ggdist::stat_halfeye(
    width = 0.6,
    ## custom bandwidth
    adjust = 0.5,
    ## move geom to the right
    justification = -.3,
    ## remove slab interval
    .width = 0,
    point_colour = NA
  ) +
  scale_x_discrete(limits=rev(markerNames)) +
  ylim(c(-0.9,3.6)) +
  geom_boxplot(
    width = .2,
    ## remove outliers
    outlier.color = NA,
    alpha = 0.5
  ) +
  # Adjust theme
  scale_fill_tq() +
  theme_tq(base_size = 20) +
  theme(legend.position="right",legend.text=element_text(size=15)) +
  scale_fill_discrete(limits = markerNames)+
  labs(
    title = "NEU",
    x    = "",
    y    = "RNAge",
    fill = ""
  ) +
  geom_point(aes(x = 4, y = selectedMarkers1[1])) +
  geom_point(aes(x = 4, y = selectedMarkers1[2])) +
  geom_point(aes(x = 4, y = selectedMarkers1[3])) +
  geom_point(aes(x = 4, y = selectedMarkers1[4])) +
  geom_point(aes(x = 4, y = selectedMarkers1[5])) +
  geom_point(aes(x = 4, y = selectedMarkers1[6])) +
  geom_point(aes(x = 4, y = selectedMarkers1[7])) +
  geom_point(aes(x = 4, y = selectedMarkers1[8])) +
  geom_point(aes(x = 2, y = selectedMarkers2[1])) +
  geom_point(aes(x = 2, y = selectedMarkers2[2])) +
  geom_point(aes(x = 2, y = selectedMarkers2[3])) +
  geom_point(aes(x = 2, y = selectedMarkers2[4])) +
  geom_point(aes(x = 2, y = selectedMarkers2[5])) +
  geom_point(aes(x = 2, y = selectedMarkers2[6])) +
  geom_point(aes(x = 2, y = selectedMarkers2[7])) +
  geom_point(aes(x = 2, y = selectedMarkers2[8])) +
  coord_flip()

dev.off()
