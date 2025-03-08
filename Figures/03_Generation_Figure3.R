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
library(CellChat)

dir_out <- "../_Data/Figures/"

## ---------------------------------- ##
## Figure 3.a
## ---------------------------------- ##

#Barplot unique genes per main cell type
SupTable4 <- read.csv("./Tables/SupplementaryTable4.csv") 
SupTable4$AnnotationMain <- gsub("myeloid|lymphoid", "immune", SupTable4$AnnotationMain)
SupTable4 %>% 
  select(AnnotationMain, gene, stage) %>%
  distinct() %>%
  ggplot(aes(x = stage, fill = factor(AnnotationMain, levels = c("immune", "endothelial", "epithelial", "mesenchymal")))) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.3))
ggsave(paste0(dir_out,"BarPlot_DEGs_main4.pdf"), 
       plot = last_plot(), 
       height = 4, width = 10)

## ---------------------------------- ##
## Figure 3.b
## ---------------------------------- ##

#Read in SupplementaryTable5
SupTable5 <- read.table("./Tables/SupplementaryTable5.txt", sep = "\t", header = TRUE)

#Filter and plot GO terms
SupTable5 %>% 
  #dplyr::select(p_value, term_size, query_size, intersection_size, term_id, term_name) %>% 
  dplyr::filter(term_size > 25 & term_size < 750) %>% #larger terms are likely to be excessively broad and prone to false positives
  distinct(cluster_id,assay,grade,term_name, .keep_all = TRUE) %>%
  #filter(grade == "all")  %>%
  filter(!grepl(pattern = "HPA|HP|CORUM|MIRNA", x = term_id)) %>%
  dplyr::arrange(p_value) %>% 
  ggplot2::ggplot() + 
  ggplot2::aes(y = reorder(term_name, -log(p_value)), x = -log(p_value), color = grade, size = intersection_size/term_size * 100) + #, shape = ordered_query
  ggplot2::geom_point() + 
  ggplot2::scale_size_continuous("Intersection size / Term size * 100 (%)") + 
  ggplot2::theme_classic() + 
  ggplot2::xlab("-Log(adj. p-value)") + 
  ggplot2::ylab("enriched terms") + 
  ggplot2::ggtitle("all cells")
ggsave(paste0(dir_out, "FigGOterms.pdf"),
       plot = last_plot(),
       width = 14, height = 14)


## ---------------------------------- ##
## Figure 3.c
## ---------------------------------- ##

cellchat <- readRDS("../_Data/07_Ligand_receptor_analysis/cellchat_CntrlvsEndo.rds")
#without previous cell filtering

Fig3c <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
Fig3c
ggsave(plot = Fig3c,
       filename = paste0(dir_out, 'Fig3c_rankNet.pdf'), width =10,height=12) 

## ---------------------------------- ##
## Figure 3.d
## ---------------------------------- ##

#Load data
cellchat.Endo <- readRDS("../_Data/07_Ligand_receptor_analysis/Endo_CellChat_AnnotationRefined.rds")
cellchat.Cntrl <- readRDS("../_Data/07_Ligand_receptor_analysis/Control_CellChat_AnnotationRefined.rds")

#circle plot
pathways.show <- c("TNF") 
object.list <- list(Endo = cellchat.Endo, Cntrl = cellchat.Cntrl)
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 15, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#Save plot
pdf(file = paste0(dir_out,"Fig3d_EndoCtl_TNF_circle_plot_edge15.pdf"))
par(mfrow = c(1, length(object.list)), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 15, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


## ---------------------------------- ##
## Figure 3.e
## ---------------------------------- ##

#Load CellChat data
cellchat <- readRDS("../_Data/07_Ligand_receptor_analysis/cellchat_CntrlvsEndo.rds")

#Figure 3e ICAM
netVisual_bubble(cellchat, comparison = c(1, 2),signaling = c("ICAM"), remove.isolate = FALSE) 
ggsave(plot = last_plot(),
       filename = paste0(dir_out, 'Fig3e_netVisual_bubble_ICAM.pdf'), width =25,height=5)


## ---------------------------------- ##
## Figure 3.f
## ---------------------------------- ##

#Load CellChat data
cellchat <- readRDS("../_Data/07_Ligand_receptor_analysis/cellchat_CntrlvsEndo.rds")

#Figure 3f CXCL
#Load data
netVisual_bubble(cellchat, comparison = c(1, 2),signaling = c("CXCL"), remove.isolate = FALSE) 
ggsave(plot = last_plot(),
       filename = paste0(dir_out, 'Fig3f_netVisual_bubble_CXCL.pdf'), width =25,height=7)

## ---------------------------------- ##
## Figure 3.g
## ---------------------------------- ##

#Load CellChat data
cellchat <- readRDS("../_Data/07_Ligand_receptor_analysis/cellchat_CntrlvsEndo.rds")

#Figure 3g VEGF
netVisual_bubble(cellchat, comparison = c(1, 2),signaling = c("VEGF"), remove.isolate = FALSE) # sources.use = 4, targets.use = c(5:11),
ggsave(plot = last_plot(),
       filename = paste0(dir_out, 'Fig3g_netVisual_bubble_VEGF.pdf'), width =15,height=5)

## ---------------------------------- ##
## Figure 3.h
## ---------------------------------- ##

#Load CellChat data
cellchat <- readRDS("../_Data/07_Ligand_receptor_analysis/cellchat_CntrlvsEndo.rds")

#Figure 3h IGFBP
netVisual_bubble(cellchat, comparison = c(1, 2),signaling = c("IGFBP"), remove.isolate = FALSE,targets.use = c(5:11)) 
ggsave(plot = last_plot(),
       filename = paste0(dir_out, 'Fig3h_netVisual_bubble_IGFBP.pdf'), width =15,height=5)
