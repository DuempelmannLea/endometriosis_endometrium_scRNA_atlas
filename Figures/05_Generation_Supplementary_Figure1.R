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

# set dir_out
dir_out <- "../_Data/Figures/"

## ---------------------------------- ##
##  Figure 1.a & b 
## ---------------------------------- ##

#Load data
EndoAtlas <- readRDS("../_Data/EndoAtlas.rds")
EndoAtlas_meta <- EndoAtlas@meta.data

#Plot and save
VlnPlot(EndoAtlas,
        features = c("nFeature_RNA", "nCount_RNA","nCount_SCT","nFeature_SCT"),
        group.by='EndometriosisStatus', ncol = 1, pt.size=0)
ggsave(paste0(dir_out, "VlnPlot_GroupByEndometriosisStatus_nFeature_nCounts_RNA_SCT_percMT.pdf"),
       plot = last_plot(),
       width = 4, height = 13)

#Percentage of main 5
table <- as.data.frame(table(EndoAtlas_meta$AnnotationMain))
table$Percentage <- round((table$Freq / sum(table$Freq)) * 100, 1)

#Number of cells for CTL and Endo
table(EndoAtlas_meta$EndometriosisStatus)

#Number of patients for CTL and Endo
EndoAtlas_meta %>% select(sample,EndometriosisStatus) %>%
  distinct() %>%
  group_by(EndometriosisStatus) %>%
  summarize(number_of_rows = n())

#mean number of cells per sample
mean(table(EndoAtlas_meta$sample))

#median count/feature
median_summary <- EndoAtlas_meta %>%
  select(nCount_RNA, nCount_SCT, nFeature_RNA, nFeature_SCT) %>%
  summarize(
    median_nCount_RNA = median(nCount_RNA),
    median_nCount_SCT = median(nCount_SCT),
    median_nFeature_RNA = median(nFeature_RNA),
    median_nFeature_SCT = median(nFeature_SCT)
  )

## ---------------------------------- ##
##  Supplementary Data Figure 1.c
## ---------------------------------- ##

#Load data
EndoAtlas <- readRDS("../_Data/EndoAtlas.rds")

#Plot and save
SFig1c <- FeaturePlot(EndoAtlas, features = "CycleDay", raster = FALSE)
ggsave(filename = paste0(dir_out, 'SFig1c.pdf'),
       plot = SFig1c,
       width =20, height=9)

## ---------------------------------- ##
##  Supplementary Data Figure 1.d
## ---------------------------------- ##

#Load data
EndoAtlas <- readRDS("../_Data/EndoAtlas.rds")

#Plot and save
SFig1d <- FeaturePlot(EndoAtlas, features = "ProgesteroneSerum", raster = FALSE)
ggsave(filename = paste0(dir_out, 'SFig1d.pdf'),
       plot = SFig1d,
       width =20, height=9)

## ---------------------------------- ##
##  Supplementary Data Figure 1.e
## ---------------------------------- ##

#Load data
p <- readRDS(paste0(dir_ENDO,"_Data/03_PCA_analysis/p_EndoAtlas_DEGsamples.rds"))

#Plot and save
pdf(paste0(dir_out,"SFig1e_eigencorplot_final_BH.pdf"), width = 10, height = 5)
eigencorplot(p,
             metavars = c("sample","EndometriosisStatus","EndometriosisGrade","CycleDay","ProgesteroneSerum","MenstrualCyclePhase","MenstrualCyclePhase_main","sequencing.batch", "X10x.capturing.batch"),
             col = c('darkorange3','darkorange', 'white','cadetblue2','cadetblue'),
             corMultipleTestCorrection = "BH"
)
dev.off()

## ---------------------------------- ##
##  Supplementary Data Figure 1.f
## ---------------------------------- ##

#Define subsets of main annotation
subsets <- c("myeloid", "endothelial", "epithelial", "lymphocyte","stromal")

#Plot eigencorplots for the main annotations
plotting_PCA_main <- function(SUBSET) {
  p <- readRDS(paste0(dir_ENDO,"_Data/03_PCA_analysis/p_EndoAtlas_",SUBSET,".rds"))
  
  pdf(paste0(dir_out,"SFig1f_eigencorplot_final_BH_",SUBSET,".pdf"), width = 10, height = 5)
  eigencorplot(p,
               metavars = c("sample","EndometriosisStatus","EndometriosisGrade","CycleDay","ProgesteroneSerum","MenstrualCyclePhase","MenstrualCyclePhase_main","sequencing.batch", "X10x.capturing.batch"),
               col = c('darkorange3','darkorange', 'white','cadetblue2','cadetblue'),
               corMultipleTestCorrection = "BH"
  )
  dev.off()
}

lapply(subsets, plotting_PCA_main)
