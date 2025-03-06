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
##  Figure 1.a & b 
## ---------------------------------- ##

#Supplementary Data Figure 1a and 1b
#CODE OK

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
##  Figure 1.c
## ---------------------------------- ##

#Supplementary Data Figure 1c
#CODE OK

SFig1c <- FeaturePlot(EndoAtlas, features = "CycleDay", raster = FALSE)

## ---------------------------------- ##
##  Figure 1.d
## ---------------------------------- ##

#Supplementary Data Figure 1d
#CODE OK

SFig1d <- FeaturePlot(EndoAtlas, features = "ProgesteroneSerum", raster = FALSE)


## ---------------------------------- ##
##  Figure 1.e
## ---------------------------------- ##

dir_out <- paste0(dir_ENDO,"_Data/03_PCA_analysis/")

p <- readRDS(paste0(dir_ENDO,"_Data/03_PCA_analysis/p_EndoAtlas_DEGsamples.rds"))

#Figure 1e
eigencorplot(p,
             metavars = c("sample","EndometriosisStatus","EndometriosisGrade","CycleDay","ProgesteroneSerum","MenstrualCyclePhase","MenstrualCyclePhase_main","sequencing.batch", "X10x.capturing.batch"),
             col = c('darkorange3','darkorange', 'white','cadetblue2','cadetblue'),
             corMultipleTestCorrection = "BH"
)


## ---------------------------------- ##
##  Figure 1.f
## ---------------------------------- ##

#Supplementary Data Figure 1f
#CODE OK
#"/home/common/data/output/projects/ENDO/E044/A082/SCTassay/" @Lea redefine path here 

#Define subsets of main annotation
subsets <- c("myeloid", "endothelial", "epithelial", "lymphocyte","stromal")

#Plot eigencorplots for the main annotations
plotting_PCA <- function(SUBSET) {
  p <- readRDS(paste0(dir_ENDO,"_Data/03_PCA_analysis/p_EndoAtlas_",SUBSET,".rds"))
  
  pdf(paste0(dir_out,"eigencorplot_final_BH_",SUBSET,".pdf"), width = 10, height = 5)
  eigencorplot(p,
               metavars = c("sample","EndometriosisStatus","EndometriosisGrade","CycleDay","ProgesteroneSerum","MenstrualCyclePhase","MenstrualCyclePhase_main","sequencing.batch", "X10x.capturing.batch"),
               col = c('darkorange3','darkorange', 'white','cadetblue2','cadetblue'),
               corMultipleTestCorrection = "BH"
  )
  dev.off()
}

lapply(subsets, plotting_PCA)
