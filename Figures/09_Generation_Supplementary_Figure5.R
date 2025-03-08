## ---------------------------------- ##
##  Setup
## ---------------------------------- ##

# load librairies
library(ggplot2)
library(CellChat)
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

dir_out <- "../_Data/Figures/"

## ---------------------------------- ##
##  Supplemental Figure 5.a
## ---------------------------------- ##

SupTable4 <- read.csv("../Tables/SupplementaryTable4.csv") 
SupTable4$stage <- gsub("-", "_", SupTable4$stage)
SupTable5 <- read.table("../Tables/SupplementaryTable5.txt", sep = "\t", header = TRUE)

#Mutate so each row contains only a single gene in the "gene" column
SupTable5 <- SupTable5 %>%
  mutate(gene = strsplit(gene.intersection.between.DEGs.and.terms, ",")) %>%
  unnest(gene) %>%
  dplyr::rename(stage = grade) 

#join tables
dfFig5a <- left_join(SupTable5, SupTable4[,c("gene","logFC","AveExpr","p_adj.loc","cluster_id","stage")], by = c("gene","cluster_id","stage"))

#Save vector with term names of interest
term_name_oi <- c("Immune System","VEGFA-VEGFR2 signaling","focal adhesion","Focal adhesion","extracellular space","Chemokine signaling pathway","NF-kappa B signaling pathway","TGF-beta receptor signaling")

#plot and save
dfFig5a  %>%
  dplyr::filter(term_name %in% term_name_oi) %>%
  ggplot2::ggplot() + 
  ggplot2::aes(y = reorder(gene, -log(logFC)), x = -log(p_adj.loc), color = stage, size = logFC) + 
  ggplot2::geom_point() + 
  ggplot2::scale_size_continuous("logFC") + 
  ggplot2::scale_x_continuous(limits = c(1, 6)) +
  ggplot2::theme_classic() +
  ggplot2::facet_wrap(~term_name, scales = "free_y", ncol = 13) +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
ggsave(filename = paste0(dir_out, 'DotPlot_DEGcategories_DEGs_logFC_pval.pdf'), width =25,height=12)


## ---------------------------------- ##
##  Supplemental Figure 5.b
## ---------------------------------- ##

# Supplementary Figure 5b
cellchat <- readRDS("../_Data/07_Ligand_receptor_analysis/cellchat_CntrlvsEndo.rds")
#without previous cell filtering

SFig5b <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
SFig5b
ggsave(plot = SFig5b,
       filename = paste0(dir_out, 'SFig5b_rankNet.pdf'), width =10,height=12) 


## ---------------------------------- ##
##  Supplemental Figure 5.c
## ---------------------------------- ##

CellChat_visualization <- function(CONDITION){
  #Load cellchat object
  cellchat <- readRDS(paste0(dir_out,CONDITION,"_CellChat_AnnotationRefined.rds"))
  
  #Define label order in plots
  labels.levels <- c("VSMC","Prv-STEAP4","Prv-CCL19","Prv-MYH11","fib C7","fib C7-SFRP2","eF4-CXCL14", 
                     "dS1-myofibroblast", "dS2","eF1","eF3","eF2", 
                     "mesothelial","glandular", "TP63+/KRT5+", "lumenal", "lumenal 1",  "lumenal 2","MUC5B+", "ciliated", 
                     "EC-aPCV", "EC-PCV","EC-capillary","EC-tip", "EC-HEV", "EC-artery", "LEC", 
                     "monocytes-CD16+", "monocytes-CD16-", "mast cells","pDC", "mDC","cDC1","pre-cDC2","cDC2","DC3", 
                     "M$\\Phi$1-LYVE1", "M$\\Phi$2-peritoneal", "M$\\Phi$3-APOE", "M$\\Phi$4-infiltrated", "M$\\Phi$5-activated",
                     "B cell", "plasma",  "pNK","NK1","NK2", "NK3", "T$_{Reg}$","T$_N$/T$_{CM}$", "CD4 T$_{RM}$","CD8 T$_{RM}$",  "T$_{EM}$", "CTL",  "ILC")
  cellchat<- updateClusterLabels(cellchat, new.order = labels.levels)
  table(cellchat@meta[["labels"]])
  
  ## Hierarchy plot of all pathways significant in Endo
  pathways_list <- "TGFb"
  for (PATHWAY in pathways_list) {
    pathways.show <- PATHWAY
    netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
    pdf(file = paste0(dir_out,"SFig5c_",CONDITION,"_",PATHWAY,"_circle_plot.pdf"))
    netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
    dev.off() 
  }
}


## ---------------------------------- ##
##  Supplemental Figure 5.d
## ---------------------------------- ##

#Load CellChat data
cellchat <- readRDS("../_Data/07_Ligand_receptor_analysis/cellchat_CntrlvsEndo.rds")

#Figure 3e ICAM
netVisual_bubble(cellchat, comparison = c(1, 2),signaling = c("TNF"), remove.isolate = FALSE) 
ggsave(plot = last_plot(),
       filename = paste0(dir_out, 'SFig5d_netVisual_bubble_TNF.pdf'), width =25,height=5)
