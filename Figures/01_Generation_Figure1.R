# Aim : Code used for the generation of figure 1
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
## Figure 1.c
## ---------------------------------- ##

#@Lea check these paths are correct 
EndoAtlas <- readRDS("/endometriosis_endometrium_scRNA_atlas/_Data/EndoAtlas.rds")
EndoAtlas_meta <- EndoAtlas@meta.data

order <- EndoAtlas_meta %>%
  select(AnnotationMain, AnnotationRefined) %>%
  filter(!is.na(AnnotationMain)) %>%
  distinct() %>%
  arrange(AnnotationRefined) %>%
  mutate(AnnotationMain = factor(AnnotationMain, levels = c( "endothelial","stromal","lymphocyte","epithelial","myeloid"))) %>%
  arrange(AnnotationMain)
EndoAtlas$AnnotationRefined <- factor(x = EndoAtlas$AnnotationRefined, levels = order$AnnotationRefined)
Fig1c <- DimPlot(EndoAtlas, reduction = "umap", group.by = "AnnotationRefined",  raster = FALSE)

## ---------------------------------- ##
## Figure 1.d
## ---------------------------------- ##

Fig1d_1 <- EndoAtlas_meta %>%
  group_by(AnnotationMain) %>%
  summarise(row_count = n()) %>% 
  ggplot(aes(x = AnnotationMain, y = row_count)) +
  geom_col() +
  coord_flip()

Fig1d_2 <- EndoAtlas_meta %>%
  group_by(sample, AnnotationMain) %>%
  summarise(row_count = n()) %>%
  ggplot(aes(x = AnnotationMain, y = row_count)) +
  geom_boxplot() +
  scale_y_log10() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

Fig1d_3 <- EndoAtlas_meta%>%
  filter(!is.na(AnnotationMain)) %>%
  group_by(AnnotationMain, EndometriosisGrade) %>%
  summarise(row_count = n()) %>%
  ggplot(aes(x = AnnotationMain, y = row_count, fill = EndometriosisGrade)) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip()

## ---------------------------------- ##
## Figure 1.e
## ---------------------------------- ##

#Figure 1e, 2a
#CODE running, need to extend the PCA analysis folder
dir_out <- "/home/common/data/output/projects/ENDO/E044/A072/allcells_EndoAtlas_SCTassay/"

p <- readRDS(paste0(dir_ENDO,"_Data/03_PCA_analysis/p_EndoAtlas.rds"))
#Figure 1e
pdf(paste0(dir_out,"eigencorplot_final_BH.pdf"), width = 10, height = 5)
eigencorplot(p,
             metavars = c("sample","EndometriosisStatus","EndometriosisGrade","CycleDay","ProgesteroneSerum","MenstrualCyclePhase","MenstrualCyclePhase_main","sequencing.batch", "X10x.capturing.batch"),
             col = c('darkorange3','darkorange', 'white','cadetblue2','cadetblue'),
             corMultipleTestCorrection = "BH"
)
dev.off()
