# Aim : Code used for the generation of figure 3
# Author: Lea Duempelmann

## ---------------------------------- ##
## Setup
## ---------------------------------- ##

# load librairies
library(Seurat)
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggbreak) 
library(tidyr)
library(patchwork)
library(ComplexUpset)
library(ggrepel)
library(Polychrome)
library(ComplexHeatmap)
library(gplots)
library(RColorBrewer)
library(randomcoloR)

## ---------------------------------- ##
## Figure 3.a
## ---------------------------------- ##

##Figure 3a
#ENDO/E044/A033/PaperFigures.Rmd
#Clean up!! use DEG and GO table as basis for plots
#Barplot unique genes per main cell type
RNAlistRDS005_df %>% 
  select(Symphony_Global_main5, gene, stage) %>%
  distinct() %>%
  ggplot(aes(x = stage, fill = factor(Symphony_Global_main5, levels = c("immune", "endothelial", "epithelial", "stromal")))) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.3))
ggsave(paste0(dir_out,"BarPlot_DEGs_main4.pdf"), 
       plot = last_plot(), 
       height = 4, width = 10)

## ---------------------------------- ##
## Figure 3.b
## ---------------------------------- ##

##Figure 3b
#Clean up!! use DEG and GO table as basis for plots
#ENDO/E044/A031/ExtractInterestingGenesFromGProfilerA031_RNA_AveExpr18_53202.Rmd
gProfiler_df %>% 
  #dplyr::select(p_value, term_size, query_size, intersection_size, term_id, term_name) %>% 
  dplyr::filter(term_size > 25 & term_size < 750) %>% #larger terms are likely to be excessively broad and prone to false positives
  distinct(cluster_id,analysis,assay,grade,term_name, .keep_all = TRUE) %>%
  #filter(grade == "all")  %>%
  filter(!grepl(pattern = "HPA|HP|CORUM|MIRNA", x = term_id)) %>%
  dplyr::arrange(p_value) %>% 
  #dplyr::slice_head(n = 20) %>% 
  ggplot2::ggplot() + 
  ggplot2::aes(y = reorder(term_name, -log(p_value)), x = -log(p_value), color = grade, size = intersection_size/term_size * 100) + #, shape = ordered_query
  ggplot2::geom_point() + 
  #ggplot2::scale_color_continuous(name = "-Log(adj. p-value)") + 
  ggplot2::scale_size_continuous("Intersection size / Term size * 100 (%)") + 
  ggplot2::theme_classic() + 
  ggplot2::xlab("-Log(adj. p-value)") + 
  ggplot2::ylab("enriched terms") + 
  ggplot2::ggtitle("all cells")
ggsave(paste0(dir_out, "GOterms_all_severe.pdf"),
       plot = last_plot(),
       width = 14, height = 14)


## ---------------------------------- ##
## Figure 3.c
## ---------------------------------- ##

CellChat
E044A076
#cellchat <-readRDS(paste0(dir_out,"cellchat_CntrlvsEndo_Symphony_Refined_Final.rds"))
dir_out <- "/home/common/data/output/projects/ENDO/E044/A076/"

#Figure 3c and Supplementary Figure 5b
gk1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gk1 
ggsave(plot = gk1,
       filename = paste0(dir_out, 'output_comparison_Refined/rankNet.pdf'), width =10,height=12)

## ---------------------------------- ##
## Figure 3.d
## ---------------------------------- ##

#Figure 3d (or was it done separately? xxx)
setwd(paste0(dir_out,"output_comparison_Refined"))
###circle plot
pathways.show <- c("TNF") 
object.list <- list(Endo = cellchat.Endo, Cntrl = cellchat.Cntrl)
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

# Save plot
pdf(file = "EndoCtl_TNF_circle_plot_edge10.pdf")
par(mfrow = c(1, length(object.list)), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

# Save plot
pdf(file = "EndoCtl_TNF_circle_plot_edge15.pdf")
par(mfrow = c(1, length(object.list)), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

## ---------------------------------- ##
## Figure 3.e
## ---------------------------------- ##
#Load CellChat data
cellchat <- readRDS(paste0(dir_out,"cellchat_CntrlvsEndo_Symphony_Refined_Final.rds"))

#Figure 3e ICAM
netVisual_bubble(cellchat, comparison = c(1, 2),signaling = c("ICAM"), remove.isolate = FALSE) 
ggsave(plot = last_plot(),
       filename = paste0(dir_out, 'output_comparison_Refined/netVisual_bubble_ICAM.pdf'), width =25,height=5)


## ---------------------------------- ##
## Figure 3.f
## ---------------------------------- ##

#Figure 3f CXCL
#Load data
netVisual_bubble(cellchat, comparison = c(1, 2),signaling = c("CXCL"), remove.isolate = FALSE) 
ggsave(plot = last_plot(),
       filename = paste0(dir_out, 'output_comparison_Refined/netVisual_bubble_CXCL.pdf'), width =25,height=7)

## ---------------------------------- ##
## Figure 3.g
## ---------------------------------- ##

#Figure 3g VEGF
netVisual_bubble(cellchat, comparison = c(1, 2),signaling = c("VEGF"), remove.isolate = FALSE) # sources.use = 4, targets.use = c(5:11),
ggsave(plot = last_plot(),
       filename = paste0(dir_out, 'output_comparison_Refined/netVisual_bubble_VEGF.pdf'), width =15,height=5)

## ---------------------------------- ##
## Figure 3.h
## ---------------------------------- ##

#Figure 3h IGFBP
netVisual_bubble(cellchat, comparison = c(1, 2),signaling = c("IGFBP"), remove.isolate = FALSE,targets.use = c(5:11)) 
ggsave(plot = last_plot(),
       filename = paste0(dir_out, 'output_comparison_Refined/netVisual_bubble_IGFBP.pdf'), width =15,height=5)