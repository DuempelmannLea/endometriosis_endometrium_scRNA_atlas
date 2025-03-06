## ---------------------------------- ##
##  Setup
## ---------------------------------- ##

# load librairies
library(ggplot2)
library(rlang)
library(stringr)
library(dplyr)
library(tidyr)
library(cowplot)
library(gridExtra)
library(patchwork)
library(tibble)
library(tinter)
library(UpSetR)
library(forcats)

## ---------------------------------- ##
##  Supplemental Figure 5.a
## ---------------------------------- ##

#Supplementary Data Figure 5a @Lea to do

## ---------------------------------- ##
##  Supplemental Figure 5.b-d
## ---------------------------------- ##

#Figure 3c and Supplementary Figure 5b
gk1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gk1 
ggsave(plot = gk1,
       filename = paste0(dir_out, 'output_comparison_Refined/rankNet.pdf'), width =10,height=12)


#Load CellChat data
cellchat <- readRDS(paste0(dir_out,"cellchat_CntrlvsEndo_Symphony_Refined_Final.rds"))

#Supplementary Data Figure 5d TNF
netVisual_bubble(cellchat, comparison = c(1, 2),signaling = c("TNF"), remove.isolate = FALSE) 
ggsave(plot = last_plot(),
       filename = paste0(dir_out, 'output_comparison_Refined/netVisual_bubble_TNF.pdf'), width =25,height=5)

