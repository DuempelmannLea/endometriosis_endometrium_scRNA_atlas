# Aim : Code used for the generation of figure 2
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
## Figure 2.a
## ---------------------------------- ##

p <- readRDS(paste0(dir_ENDO,"_Data/03_PCA_analysis/p_EndoAtlas.rds"))
#Figure 2a
biplot(p,
       colby = 'MenstrualCyclePhase',
       # ellipse config
       ellipse = TRUE,
       ellipseType = 't',
       ellipseLevel = 0.90,
       ellipseFill = TRUE,
       ellipseAlpha = 1/4,
       ellipseLineSize = 1.0,
       labSize = 0,drawConnectors = FALSE,
       legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0)

## ---------------------------------- ##
## Figure 2.b
## ---------------------------------- ##

EndoAtlas <- readRDS("/endometriosis_endometrium_scRNA_atlas/_Data/EndoAtlas.rds")
EndoAtlas_meta <- EndoAtlas@meta.data

#Figure 2b entire atlas
Fig2b_1 <- DimPlot(EndoAtlas, reduction = "umap", group.by = "MenstrualCyclePhase", raster = FALSE)
Fig2b_1
ggsave(filename = paste0(dir_out, 'Fig2b_1.pdf'),
       plot = Fig2b_1,
       width =20, height=9)

#Figure 2b epithelial cells
#@ Lea add how this object was generated EndoAtlas_epithelial_reintegrated
Fig2b_2 <- DimPlot(EndoAtlas_epithelial_reintegrated, reduction = "umap", group.by = "MenstrualCyclePhase", raster = FALSE)

## ---------------------------------- ##
## Figure 2.c
## ---------------------------------- ##

#Figure 2c
#@LeaLook in IBU folder for object with trajectory line or A063/A066

#Figure 2c epithelial cells
epithelial_cells_monocle3 <- readRDS(paste0(dir_data, "epithelial_cells_monocle3.rds"))

#Pseudotime
Fig2c <- plot_cells(epithelial_cells_monocle3, color_cells_by = "pseudotime",
                    label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE,
                    graph_label_size=3, cell_size = 1)

## ---------------------------------- ##
## Figure 2.d
## ---------------------------------- ##

#Remove samples 65 and 105, since these samples are just in between menstrual cycle phases and cannot be clearly assigned
Idents(EndoAtlas) <- EndoAtlas$sample
EndoAtlas_filtereds65s105 <- subset(EndoAtlas, subset = sample %in% setdiff(unique(EndoAtlas$sample), c("sample65", "sample105")))
Idents(EndoAtlas_filtereds65s105) <- EndoAtlas_filtereds65s105$AnnotationMain

#### Plot heatmaps for main cell type specific markers for the different menstrual cycle phases. Compare also Wang et al 2020 marker genes (Wang et al. 2020, DOI: 10.1038/s41591-020-1040-z)

###########################
### main epithelial marker genes for menstrual cycle phases

# Subset celltype of interest
EndoAtlas_filtereds65s105 <- subset(EndoAtlas_filtereds65s105, idents = epithelial)
#main epithelial marker genes for menstrual cycle phases
epithelial_main <- c("SFRP4", "LAMC2", "MMP11", "MMP7", "SPRY1", "TFPI2","STEAP4","ADAMTS8","CAPN6","SLC37A2","SCGB1D2", "MT1F", "MT1X","RIMKLB", "S100P", "DEFB1","G0S2","NNMT","GPX3", "PAEP")

# Do MultiBarHeatmap
DoMultiBarHeatmap(object = data, features = epithelial_main, group.by="MenstrualCyclePhase", additional.group.by="epithelial.pseudotime")
# Save MultiBarHeatmap
ggsave(
  filename = paste0(dir_out, epithelial, "Fig2d_epithelial_DoMultiBarHeatmap.pdf"),
  plot = last_plot(),
  width = 15,
  height = 8)

##Plot our marker genes and Wang et al. 2020 marker genes (Wang et al. 2020, DOI: 10.1038/s41591-020-1040-z)
#Wang et al. 2020 menstrual phase marker genes
Wang_epithelial_mens_markers <- c("PLAU", "MMP7", "THBS1", "CADM1", "NPAS3", "ATP1A1", "ANK3", "ALPL", "TRAK1", "SCGB1D2", "MT1F", "MT1X", "MT1E", "MT1G", "CXCL14", "MAOA", "DPP4", "NUPR1", "GPX3", "PAEP")
#Plot first our genes, then Wang genes
heatmap <- DoHeatmap(
  object = data,
  features = c(setdiff(epithelial_main ,Wang_epithelial_mens_markers), Wang_epithelial_mens_markers),
  cells = NULL,
  group.by = "MenstrualCyclePhase",
  group.bar = TRUE)
#heatmap
ggsave(
  filename = paste0(dir_out, epithelial, "DoMultiBarHeatmap_epithelial_OursAndWangMarkers.pdf"),
  plot = last_plot(),
  width = 15,
  height = 8)


###########################
### main stromal marker genes for menstrual cycle phases

# Subset celltype of interest
EndoAtlas_filtereds65s105 <- subset(EndoAtlas_filtereds65s105, idents = stromal) #@ Lea add the path here 
#main stromal marker genes for menstrual cycle phases
stromal_main <- c("SFRP1", "H1-1", "H3C2", "WNT5A","CILP","BRINP1", "CFD","PTGDS","SLIT3","RPRM", "TIMP3", "IGFBP1","LEFTY2","LRRC15","LUM")
# Do MultiBarHeatmap
DoMultiBarHeatmap(object = data, features = stromal_main, group.by="MenstrualCyclePhase", additional.group.by="epithelial.pseudotime")
# Save MultiBarHeatmap
ggsave(
  filename = paste0(dir_out, stromal, "Fig2d_stromal_DoMultiBarHeatmap.pdf"),
  plot = last_plot(),
  width = 15,
  height = 10)

##Plot our marker genes and Wang et al. 2020 marker genes (Wang et al. 2020, DOI: 10.1038/s41591-020-1040-z)
#Wang et al. 2020 menstrual phase marker genes
Wang_stromal_mens_markers <- c('STC1', 'NFATC2', 'BMP2' ,'PMAIP1' ,'MMP11' ,'SFRP1', 'WNT5A' ,'ZFYVE21' ,'CILP' ,'SLF2' ,'MATN2', 'S100A4' ,'DKK1' ,'CRYAB', 'FOXO1', 'IL15', 'FGF7', 'LMCD1')
#Plot first our genes, then Wang genes
heatmap <- DoHeatmap(
  object = data,
  features = c(setdiff(stromal_main ,Wang_stromal_mens_markers), Wang_stromal_mens_markers),
  cells = NULL,
  group.by = "MenstrualCyclePhase",
  group.bar = TRUE)
#heatmap
ggsave(
  filename = paste0(dir_out, stromal, "DoMultiBarHeatmap_stromal_OursAndWangMarkers.pdf"),
  plot = last_plot(),
  width = 15,
  height = 15)


###########################
### main endothelial, myeloid and lymphocyte marker genes for menstrual cycle phases
endothelial_main <- c("RGS16","DUOX1","CFD","MYC","TIMP3","LEFTY2","IL1B")
myeloid_main <- c("S100A10","TREM2","STXBP2","CCL2","GNLY","CXCL3","CXCL2")
lymphocyte_main <- c("IL7R","DUSP4","MMP26","MT1G","CXCL3","CXCL2","G0S2")

#Make heatmaps in loop
CELLTYPES <- c("endothelial","lymphocyte","myeloid") 
for (CELLTYPE in CELLTYPES) {
  # Subset celltype of interest
  EndoAtlas_filtereds65s105 <- subset(EndoAtlas_filtereds65s105, idents = CELLTYPE)
  # Do MultiBarHeatmap
  DoMultiBarHeatmap(object = EndoAtlas_filtereds65s105, 
                    features = if (CELLTYPE == "endothelial") endothelial_main
                    else if (CELLTYPE == "myeloid") myeloid_main
                    else if (CELLTYPE == "lymphocyte") lymphocyte_main, 
                    group.by="MenstrualCyclePhase", 
                    additional.group.by="epithelial.pseudotime")
  # Save MultiBarHeatmap
  ggsave(
    filename = paste0(dir_out, "Fig2d_", CELLTYPE, "_DoMultiBarHeatmap.pdf"),
    plot = last_plot(),
    width = 15,
    height = 3)
  
}

###########################
## Additional strong Menstrual cycle phase markers
stromal_extended <- c("H1-1","H1-2","H1-4","H1-5","H2AC13","H2Ac14","H3C2","H4C3","MKI67","TNC", #prolif/periov
                      "CALB2", "PENK","POSTN",#periovulatory
                      "CAPN6", "CYP26A1","HGD","HLA-DMB","MMP26","UPK1B","TFPI2","ENPP3","PKHD1L1","CILP","CHODL","GDF15","SCGB1D2","SCGB1D4","SCGB2A1","CERNA2","CD9",#secretory_early
                      "MT1E","MT1F","MT1G","MT1H","MT1M", "LINC01320",#secretory_early/_mid
                      "ACKR1", #secretory_mid
                      "BRINP1","C1QC","CCN3","CFD","GPX3","MPTL2","MUC16","MYC","PAEP","PDK4","PTGDS","RPRM","SLIT3","SLPI","SPP1","TGM2","TIMP3","TM1E","TM4SF1",
                      "IGFBP1","IGFBP6","KRT7","LEFTY2","LMCD1","LRRC15","LUM","NID2","PLAGL1","PLAT","PRLR","RBP4","BAMBI","CAB39L","COMP", #secretory_late
                      "APOD","CLDN4","CLU","CP","CXCL13","CXCL14","DEFB1","DEPP1","DKK1") #secretory_early/_mid/_late
epithelial_extended <- unique(c("SFRP4", "LAMC2", "MMP11", "TMEM107", "C12orf75", "UBE2C", "CENPW", "H2AC11", "H3C3", "H2AC13", "CSF3", "CENPE","SPRY1", "TFPI2", "STEAP4", "CKB", "LRRC26", "SNX30", "ADAMTS8", "CAPN6", "FXYD4", "XDH","CCN3", "SLC37A2","MMP26", "SCGB1D2", "SCGB1D4", "CYP26A1", "ENPP3", "PTGDS",  "PLEKHH2", "ZBTB16", "STC1", "IGFBP5", "PDGFRA", "REV3L", "HSPA1A",   "MT1F", "MT1G", "PAEP", "RIMKLB", "S100P", "DEFB1", "G0S2", "GDF15", "C12orf75", "DEFB1", "SPP1", "DKK1", "NNMT", "RBP4"))
endothelial_extended <- c("POSTN","RGS16","TFPI2", "SCGB1D4",  "CYP26A1","MT1G", "CFD", "FABP5","APLN","TIMP3","MYC","MUC16", "NAPSB", "LEFTY2","IL1B", "COMP")
myeloid_extended <- c("TREM2",  "PHYHIPL", "PAEP", "RIMKLB", "LEFTY2", "IL2RB", "APOD", "NNMT",  "STXBP2","GNLY", "S100A10", "OLR1", "MUC16", "CXCL13",  "INSIG1", "GZMA",  "PLAUR",  "MMP26",  "RGCC", "IGSF6",  "C15orf48", "UBD", "CCL19",  "FBP1", "CCL2", "FCN1", "CXCL3", "CORO1A","MT1G","INHBA", "CXCL5", "MMP10","SERPINB2","MMP1", "CCL18","SLC7A5","CCL20","THBS1","AREG", "APOC1","SERPINA1")
lymphocyte_extended <- c("CAPN6","MT1H",  "ENPP3",  "TFPI2", "CLDN3", "DUSP4", "PDK4", "MMP26", "CXCL13", "GIMAP4",  "RIMKLB", "MT1G", "DEFB1", "CFD","IL7R", "LDLRAD4", "LEFTY2", "DOCK10", "CYP26A1", "G0S2","MTND1P23")

## ---------------------------------- ##
## Figure 2.e
## ---------------------------------- ##

#Prepare data for boxplot
df2 <- EndoAtlas_meta%>%
  filter(Minor.exclusion.criteria.logical == "FALSE") %>%
  group_by(sample, AnnotationMain, EndometriosisGrade, MenstrualCyclePhase_main, EndometriosisStatus) %>%
  summarise(row_count = n()) %>%
  ungroup() %>%
  filter(!is.na(AnnotationMain)) %>%
  group_by(sample) %>%
  mutate(AnnotationMain_perc = row_count / sum(row_count) * 100) %>%
  ungroup()

#BoxPlots percentage, EndometriosisStatus 
ggplot(df2, aes(x = MenstrualCyclePhase_main, y = AnnotationMain_perc, fill = EndometriosisStatus)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(colour = EndometriosisStatus),
             colour="black",pch=21, 
             position = position_jitterdodge()) +
  RotatedAxis() +
  xlab('Sample ID') +
  facet_wrap(~AnnotationMain,scales = "free", ncol = 6)
ggsave(paste0(dir_out, "Fig2e_boxplot.pdf"),
       plot = last_plot(),
       width = 12, height = 4.4)


###Significance testing CTL vs Endo
# Perform t-test CTL vs Endo for each cell type grouped by AnnotationRefined, MenstrualCyclePhase_main
df2_ttest_EndometriosisStatus <- df2 %>%
  group_by(AnnotationMain, MenstrualCyclePhase_main) %>%
  summarise(p.value = t.test(AnnotationMain_perc ~ EndometriosisStatus)$p.value) %>%
  ungroup()
#epithelial, proliferative just not significant (p.value = 0.0898)

###Significance testing prolif vs periov vs sec, CTL + Endo
## Perform t-test Prolif vs Peri for each cell type
df2_ttest_ProlifPeri <- df2 %>%
  dplyr::filter(MenstrualCyclePhase_main %in% c("proliferative","periovulatory")) %>%
  group_by(AnnotationMain) %>%
  summarise(p.value = t.test(AnnotationMain_perc ~ MenstrualCyclePhase_main)$p.value) %>%
  ungroup()
#none significant
## Perform t-test Prolif vs Sec for each cell type
df2_ttest_ProlifSec <- df2 %>%
  dplyr::filter(MenstrualCyclePhase_main %in% c("proliferative","secretory")) %>%
  group_by(AnnotationMain) %>%
  summarise(p.value = t.test(AnnotationMain_perc ~ MenstrualCyclePhase_main)$p.value) %>%
  ungroup()
#epi (0.00977) and str (0.0255) significant

## Perform t-test periovulatory vs sec for each cell type
df2_ttest_PeriSec <- df2 %>%
  dplyr::filter(MenstrualCyclePhase_main %in% c("periovulatory","secretory")) %>%
  group_by(AnnotationMain) %>%
  summarise(p.value = t.test(AnnotationMain_perc ~ MenstrualCyclePhase_main)$p.value) %>%
  ungroup()
#epi (0.00794) and str (0.0137) significant


## ---------------------------------- ##
## Figure 2.f
## ---------------------------------- ##

#Early and Mid Secretory Prv_VSMCs 
Fig2f <- EndoAtlas_meta%>%
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
Fig2f <- left_join(Fig2f, dplyr::distinct(EndoAtlas_meta[,c("sample","EndometriosisStatus","MenstrualCyclePhase_main","MenstrualCyclePhase")]), by = "sample")

secEarylMid_Prv_VSMC_MenstrualCyclePhase_main <- Fig2f %>%
  filter(AnnotationRefined %in% c("Prv_VSMC secretory", "Prv_VSMC secretory_late")) %>% #same result with and without "Prv_VSMC secretory_late"
  dplyr::filter(MenstrualCyclePhase %in% c("proliferative","periovulatory","secretory_mid","secretory_early")) %>%
  group_by(sample, EndometriosisStatus, MenstrualCyclePhase_main) %>%
  summarize(
    total_counts_MenstrualCyclePhase_main = sum(AnnotationRefined_counts, na.rm = TRUE),
    total_perc_MenstrualCyclePhase_main = sum(AnnotationRefined_perc, na.rm = TRUE)
  ) %>%
  ungroup()

ggplot(secEarylMid_Prv_VSMC_MenstrualCyclePhase_main, aes(x = MenstrualCyclePhase_main, y = total_perc_MenstrualCyclePhase_main, fill = EndometriosisStatus)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(colour = EndometriosisStatus),
             colour="black",pch=21, 
             position = position_jitterdodge(jitter.width = 0.2))
ggsave(paste0(dir_out, "Fig2f_boxplot.pdf"),
       plot = last_plot(),
       width = 17, height = 12)

#Statistical testing
secEarylMid_Prv_VSMC_MenstrualCyclePhase_main_pval <- secEarylMid_Prv_VSMC_MenstrualCyclePhase_main %>%
  dplyr::filter(MenstrualCyclePhase_main %in% c("secretory")) %>% 
  summarise(p.value = t.test(total_perc_MenstrualCyclePhase_main ~ EndometriosisStatus)$p.value) %>%
  ungroup()
#significant, p.value 0.0261

secEarylMid_Prv_VSMC_MenstrualCyclePhase_main_pval_prolif <- secEarylMid_Prv_VSMC_MenstrualCyclePhase_main %>%
  dplyr::filter(MenstrualCyclePhase_main %in% c("proliferative")) %>% 
  summarise(p.value = t.test(total_perc_MenstrualCyclePhase_main ~ EndometriosisStatus)$p.value) %>%
  ungroup()
#not significant, p.value NA

secEarylMid_Prv_VSMC_MenstrualCyclePhase_main_pval_peri <- secEarylMid_Prv_VSMC_MenstrualCyclePhase_main %>%
  dplyr::filter(MenstrualCyclePhase_main %in% c("periovulatory")) %>% 
  summarise(p.value = t.test(total_perc_MenstrualCyclePhase_main ~ EndometriosisStatus)$p.value) %>%
  ungroup()
#not significant, p.value 0.7639853 

#Multiple testing correction
p.adjust(c(0.0261), method = "BH", n = length(3)) #Value does not change
p.adjust(c(0.7639853), method = "BH", n = length(3)) #Value does not change