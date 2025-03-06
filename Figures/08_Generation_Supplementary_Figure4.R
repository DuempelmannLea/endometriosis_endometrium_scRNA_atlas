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
##  Figure 4.a
## ---------------------------------- ##

#Supplementary Data Figure 4a
#CODE OK

##Barplots of cell type frequency split by the main 5 cell types

#Set colour palette
P63 = createPalette(63, c("#ff0000", "#00ff00", "#0000ff")) 
swatch(P63)
names(P63) <- NULL


#Make df with number of AnnotationRefined cells in each sample
df <- EndoAtlas_meta%>%
  filter(!is.na(AnnotationMain)) %>%
  group_by(sample, AnnotationMain, AnnotationRefined, EndometriosisGrade, MenstrualCyclePhase, epithelial.pseudotime, EndometriosisStatus) %>%
  summarise(row_count = n()) %>%
  ungroup() 

#Make bar plot
p <- ggplot(df, aes(x = sample, y = row_count, fill = AnnotationRefined)) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = sample(P63)) +
  xlab('Sample ID') +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "lightgray", size = 0.5)
  ) +
  facet_grid(AnnotationMain ~ ., scales = "free_y")
p
#Save barplot
ggsave(paste0(dir_out, "SFig4a_barplot.pdf"),
       plot = p,
       width = 15, height = 6)

## ---------------------------------- ##
##  Figure 4.b
## ---------------------------------- ##

#Supplementary Data Figure 4b
#CODE OK

##########################################################################################
##Plot differences for samples without minor exclusion criteria

#Make df 
#Introduce 0 for the ones with 0 count!
SFig4b <- EndoAtlas_meta%>%
  filter(Minor.exclusion.criteria.logical == "FALSE") %>%
  group_by(sample, AnnotationRefined) %>%
  summarise(AnnotationRefined_counts = n()) %>%
  ungroup() %>%
  complete(sample, AnnotationRefined, fill = list(AnnotationRefined_counts = 0)) %>%
  filter(!is.na(AnnotationRefined)) %>%
  group_by(sample) %>%
  mutate(AnnotationRefined_perc = AnnotationRefined_counts / sum(AnnotationRefined_counts) * 100) %>%
  ungroup() %>% #60 samples * 63 cell types = 3780 rows, correct
  filter(!sample %in% c("sample81", "sample87", "sample108", "sample115", "sample119", "sample121"))  #54 samples * 63 cell types = 3402 rows, correct
SFig4b <- left_join(SFig4b, dplyr::distinct(EndoAtlas_meta[,c("sample","EndometriosisStatus","MenstrualCyclePhase_main")]), by = "sample")
write.csv(SFig4b, paste0(dir_out,"SFig4b_AnnotationRefined_frequency.csv"))

#BoxPlots percentage, EndometriosisStatus, percentage ###for paper Extended Data Fig. S4
ggplot(SFig4b, aes(x = MenstrualCyclePhase_main, y = AnnotationRefined_perc, fill = EndometriosisStatus)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(colour = EndometriosisStatus),
             colour="black",pch=21, 
             position = position_jitterdodge(jitter.width = 0.2)) +
  facet_wrap(~AnnotationRefined, scales = "free", ncol = 7) +
  theme_bw() +
  expand_limits(y = 0)

ggsave(paste0(dir_out, "SFig4b_EndometriosisStatus_perc.pdf"),
       plot = last_plot(),
       width = 20, height = 12)

## Perform t-test Prolif vs Sec for each cell type
SFig4b_ttest_ProlifSec <- SFig4b %>%
  dplyr::filter(MenstrualCyclePhase_main %in% c("proliferative","secretory")) %>%
  group_by(AnnotationRefined) %>%
  summarise(p.value = t.test(AnnotationRefined_perc ~ MenstrualCyclePhase_main)$p.value) %>%
  ungroup()
SFig4b_ttest_ProlifSec$p.value_adjBH <- p.adjust(SFig4b_ttest_ProlifSec$p.value, method = "BH", n = length(SFig4b_ttest_ProlifSec$p.value))
write.csv(SFig4b_ttest_ProlifSec, paste0(dir_out,"SFig4b_ttest_ProlifSec_adjpval.csv"))

# Perform t-test CTL vs Endo for each cell type grouped by AnnotationRefined, MenstrualCyclePhase_main
SFig4b_ttest_EndometriosisStatus <- SFig4b %>%
  group_by(AnnotationRefined, MenstrualCyclePhase_main) %>%
  summarise(p.value = t.test(AnnotationRefined_perc ~ EndometriosisStatus)$p.value) %>%
  ungroup()
#Prv_VSMC secretory pval 0.009075467
SFig4b_ttest_EndometriosisStatus$p.value_adjBH <- p.adjust(SFig4b_ttest_EndometriosisStatus$p.value, method = "BH", n = length(SFig4b_ttest_EndometriosisStatus$p.value))
write.csv(SFig4b_ttest_EndometriosisStatus, paste0(dir_out,"SFig4b_ttest_EndometriosisStatus_adjpval.csv"))
#After correction for multiple testing Prv_VSMC secretory is not significant anymore
