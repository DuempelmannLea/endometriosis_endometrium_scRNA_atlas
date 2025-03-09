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

# input path 
input_path<- '../Tables' 
  
# source common functions 
source('./common_functions.R')

# set dir_out
dir_out <- "../_Data/Figures/"

## ---------------------------------- ##
##  Supplementary Data Figure 2
## ---------------------------------- ##

#Load data
EndoAtlas <- readRDS("../_Data/EndoAtlas.rds")

# Set SCT as default and remove RNA to save space
DefaultAssay(EndoAtlas) <- "SCT"
EndoAtlas@assays$RNA <- NULL
gc()

# Remove data without annotation (no NAs!)
Idents(EndoAtlas) <- EndoAtlas$AnnotationMain
EndoAtlas <- subset(EndoAtlas, subset = AnnotationMain %in% levels(as.factor(EndoAtlas$AnnotationMain)))
gc()

### AnnotationMain markers
genes_to_plot <- c("EPCAM","CDH1","PDGFRA","COL1A1", "VWF","PECAM1","PTPRC","CD68","CD14","CD2")
pops_to_plot <- c("stromal","epithelial","myeloid","lymphocyte","endothelial")
# DotPlot
p_main <- Seurat::DotPlot(
  object = EndoAtlas,
  scale = T,
  features = genes_to_plot,
  idents = pops_to_plot,
) + theme(axis.text.x = element_text(angle = 90))
# Save DotPlot
ggsave(plot = p_main, "../_Data/Figures/SFig2_AnnotationMain.pdf", width = 10, height = 5)


###Set identity to AnnotationRefined
Idents(EndoAtlas) <- EndoAtlas$AnnotationRefined

### DC markers (from Tan2022 Fig. 5a)
genes_to_plot <- c("CD1C", "CLEC9A", "LAMP3")
pops_to_plot <- c("mDC", "cDC1", "DC3", "cDC2", "pre-cDC2")
# reorder clusters
pops_seurat <- subset(EndoAtlas, idents = pops_to_plot)
pops_seurat@active.ident <- factor(pops_seurat@active.ident, 
                                   levels=pops_to_plot)
# DotPlot
p_DC <- Seurat::DotPlot(
  object = pops_seurat,
  scale = T,
  features = genes_to_plot,
  idents = pops_to_plot,
) + theme(axis.text.x = element_text(angle = 90))
# Save DotPlot
ggsave(plot = p_DC, "../_Data/Figures/SFig2_DCs_SCT.pdf", width = 5, height = 7)


### Lymphocytes markers (from Tan2022 Fig. 4b)
genes_to_plot <- c("FCGR3A","KLRD1","ITGA1","LILRB1","ITGAX","KLRC1","KLRC2","EPAS1","KIR2DL4","IGFBP2","CD160","GZMH","KLRF1","FGFBP2","TBX21","CD8A","CD8B","CRTAM","CCR5","NCR3","RORC","IL7R","TCF7","CD4","CCR7","SELL","LEF1","FOXP3","CTLA4","IL2RA","CXCR6","ZNF683","ITGAE","XCL1","SLC4A10","DPP4","GZMK","IGKC", "MZB1", "DERL3", "CD79A", "CD79B", "MS4A1", "BANK1")
pops_to_plot <- c(EndoAtlas@meta.data %>% filter(AnnotationMain == "lymphocyte") %>% pull(AnnotationRefined) %>% unique())

# reorder clusters
pops_seurat <- subset(EndoAtlas, idents = pops_to_plot)
pops_seurat@active.ident <- factor(pops_seurat@active.ident, 
                                   levels=c("B cell", "plasma", "CD8 T$_{RM}$","CD4 T$_{RM}$","T$_{Reg}$", "T$_N$/T$_{CM}$" , "ILC" , "T$_{EM}$","CTL", "pNK", "NK3","NK2","NK1"))
# DotPlot
p_lymph <- Seurat::DotPlot(
  object = pops_seurat,
  scale = T,
  features = genes_to_plot,
  idents = pops_to_plot,
) + theme(axis.text.x = element_text(angle = 90))
# Save DotPlot
ggsave(plot = p_lymph, "../_Data/Figures/SFig2_Lymphocytes_SCT.pdf", width = 17, height = 7)

### epithelial markers (from Tan2022 Sup Fig. 7c)
genes_to_plot <- c("TPPP3", "PIFO", "FOXJ1", "TP63", "KRT5", "SOX9", "LGR5", "MUC5B", "TFF3", "RUNX3", "SAA1", "SIX1", "PROM1", "PROM2", "ANPEP", "KIAA1324", "SPDEF", "SMAD9", "CD36", "HSD17B2", "SCGB2A2", "GDA", "FGFR2", "WNT7A", "FGF9", "PTGS1", "MSLN", "UPK3B", "CALB2", "PGR", "ESR1", "MMP7", "TIMP1", "IDO1", "PSAT1", "ENPP3", "GNG11", "CREB3L1", "IHH", "TOP2A", "MT1F", "MT1E", "MT1X", "SERPINA5", "SPP1", "PAEP", "CXCL14", "DPP4")
pops_to_plot <- c("mesothelial", "ciliated", "mid-secretory", "lumenal 2", "lumenal 1", "lumenal", "glandular early-secretory", "glandular", "TP63+/KRT5+", "MUC5B+")
# reorder clusters
pops_seurat <- subset(EndoAtlas, idents = pops_to_plot)
pops_seurat@active.ident <- factor(pops_seurat@active.ident, 
                                   levels=pops_to_plot)
#DotPlot
p_epi <- Seurat::DotPlot(
  object = pops_seurat,
  scale = T,
  features = genes_to_plot,
  idents = pops_to_plot,
) + theme(axis.text.x = element_text(angle = 90))
# Save DotPlot
ggsave(plot = p_epi, "../_Data/Figures/SFig2_epithelial_SCT.pdf", width = 17, height = 7)


#myeloid markers ( from Tan2022 Fig. 4b)
genes_to_plot <- c("FLOR2", "LYVE1", "MRC1", "CX3CR1", "ICAM2", "CLEC5A", "CCR2", "VEGFA", "FN1", "SPARC", "APOE", "SPP1")
pops_to_plot <- c(EndoAtlas@meta.data %>% filter(AnnotationMain == "myeloid") %>% pull(AnnotationRefined) %>% unique())
pops_to_plot <- c("M$\\Phi$5-activated","M$\\Phi$4-infiltrated","M$\\Phi$3-APOE","M$\\Phi$2-peritoneal","M$\\Phi$1-LYVE1")
# reorder clusters
pops_seurat <- subset(EndoAtlas, idents = pops_to_plot)
pops_seurat@active.ident <- factor(pops_seurat@active.ident, 
                                   levels=pops_to_plot)
# DotPlot
p_my <- Seurat::DotPlot(
  object = pops_seurat,
  scale = T,
  features = genes_to_plot,
  idents = pops_to_plot,
) + theme(axis.text.x = element_text(angle = 90))
# Save DotPlot
ggsave(plot = p_my, "../_Data/Figures/SFig2_macrophages_SCT.pdf", width = 10, height = 7)


#Prv and VSMC markers (from Tan2022 Fig. 3b)
genes_to_plot <- c("RERGL", "MCAM", "MYH11", "STEAP4","GGT5", "RGS5")
pops_to_plot <- c("Prv-MYH11", "Prv-CCL19", "Prv-STEAP4", "VSMC")
# reorder clusters
pops_seurat <- subset(EndoAtlas, idents = pops_to_plot)
pops_seurat@active.ident <- factor(pops_seurat@active.ident, 
                                   levels=c("VSMC", "Prv-STEAP4", "Prv-CCL19","Prv-MYH11"))
# DotPlot
p_Prv <- Seurat::DotPlot(
  object = pops_seurat,
  scale = T,
  features = genes_to_plot,
  idents = pops_to_plot,
) + theme(axis.text.x = element_text(angle = 90))
# Save DotPlot
ggsave(plot = p_Prv, "../_Data/Figures/SFig2_Prv_VSMC_SCT.pdf", width = 7, height = 7)
