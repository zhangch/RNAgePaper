#####R version: 3.2.1#####
library(tidyverse)
library(ggdist)
library(tidyquant)

#############################################
## S. Fig 3A
#############################################
tag = "mean"

markerNames = c("Pan markers", "pF markers", "Neuronal markers", "FC markers")

pdf(paste0("Figure6B_V4_", tag, ".pdf"), width = 12, height = 4, useDingbats = FALSE)
detail.1HAE <- read.csv("../Data/L1000/NPC_allData.csv", header=T, stringsAsFactors=FALSE)

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
  ylim(c(-0.8,1.6)) +
  geom_boxplot(
    width = .2,
    ## remove outliers
    outlier.color = NA,
    alpha = 0.5
  ) +
  # Adjust theme
  scale_fill_tq() +
  theme(legend.position="right",legend.text=element_text(size=15)) +
  theme_tq(base_size = 20) +
  scale_fill_discrete(limits = markerNames)+
  labs(
    title = "NPC",
    x    = "",
    y    = "RNAge",
    fill = ""
  ) +
  coord_flip()

dev.off()
