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

dir_out <- "../_Data/Figures/"

## ---------------------------------- ##
## Figure 1.c
## ---------------------------------- ##

EndoAtlas <- readRDS("../_Data/EndoAtlas.rds")

order <- EndoAtlas@meta.data %>%
  select(AnnotationMain, AnnotationRefined) %>%
  filter(!is.na(AnnotationMain)) %>%
  distinct() %>%
  arrange(AnnotationRefined) %>%
  mutate(AnnotationMain = factor(AnnotationMain, levels = c( "endothelial","stromal","lymphocyte","epithelial","myeloid"))) %>%
  arrange(AnnotationMain)
EndoAtlas$AnnotationRefined <- factor(x = EndoAtlas$AnnotationRefined, levels = order$AnnotationRefined)
Fig1c <- DimPlot(EndoAtlas, reduction = "umap", group.by = "AnnotationRefined",  raster = FALSE)
ggsave(filename = paste0(dir_out, 'Fig1c.pdf'),
       plot = Fig1c,
       width =20, height=9)

## ---------------------------------- ##
## Figure 1.d
## ---------------------------------- ##

Fig1d_1 <- EndoAtlas@meta.data %>%
  group_by(AnnotationMain) %>%
  summarise(row_count = n()) %>% 
  ggplot(aes(x = AnnotationMain, y = row_count)) +
  geom_col() +
  coord_flip()
ggsave(paste0(dir_out, "Fig1d_1.pdf"),
       plot = last_plot(),
       width = 15, height = 6)

Fig1d_2 <- EndoAtlas@meta.data %>%
  group_by(sample, AnnotationMain) %>%
  summarise(row_count = n()) %>%
  ggplot(aes(x = AnnotationMain, y = row_count)) +
  geom_boxplot() +
  scale_y_log10() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(dir_out, "Fig1d_2.pdf"),
       plot = last_plot(),
       width = 15, height = 6)


Fig1d_3 <- EndoAtlas@meta.data %>%
  filter(!is.na(AnnotationMain)) %>%
  group_by(AnnotationMain, EndometriosisGrade) %>%
  summarise(row_count = n()) %>%
  ggplot(aes(x = AnnotationMain, y = row_count, fill = EndometriosisGrade)) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip()
ggsave(paste0(dir_out, "Fig1d_3.pdf"),
       plot = last_plot(),
       width = 6, height = 4)

## ---------------------------------- ##
## Figure 1.e
## ---------------------------------- ##


p <- readRDS("../_Data/03_PCA_analysis/p_EndoAtlas.rds")
pdf(paste0(dir_out,"eigencorplot_final_BH.pdf"), width = 10, height = 5)
eigencorplot(p,
             metavars = c("sample","EndometriosisStatus","EndometriosisGrade","CycleDay","ProgesteroneSerum","MenstrualCyclePhase","MenstrualCyclePhase_main","sequencing.batch", "X10x.capturing.batch"),
             col = c('darkorange3','darkorange', 'white','cadetblue2','cadetblue'),
             corMultipleTestCorrection = "BH"
)
dev.off()
