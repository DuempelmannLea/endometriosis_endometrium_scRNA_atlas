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

dir_out <- "../_Data/Figures/"

## ---------------------------------- ##
##  Supplementary Data Figure 3.a
## ---------------------------------- ##

#Load data
p <- readRDS("../_Data/03_PCA_analysis/p_EndoAtlas.rds")

#Plot
biplot(p,
       colby = 'EndometriosisStatus',
       labSize = 0,drawConnectors = FALSE, 
       hline = 0, vline = 0,
       legendPosition = 'right')

## ---------------------------------- ##
##  Supplementary Data Figure 3.b
## ---------------------------------- ##

#Load data
epithelial_cells_monocle3 <- readRDS(paste0(dir_data, "epithelial_cells_monocle3.rds"))

#Plot Menstrual Cycle Day
SFig3b <- plot_cells(epithelial_cells_monocle3, color_cells_by = "CycleDay",
                     label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE,
                     graph_label_size=3, cell_size = 1)
SFig3b
ggsave(filename = paste0(dir_out, 'SFig3b.pdf'),
       plot = SFig3b,
       width =20, height=9)

## ---------------------------------- ##
##  Supplementary Data Figure 3.c
## ---------------------------------- ##

#Load data and clean
SupplementaryTable1 <- read.delim("../Tables/SupplementaryTable1.txt")
SupplementaryTable1 <- SupplementaryTable1 %>%
  filter(!(is.na(FinalPhaseRefined) | FinalPhaseRefined == "")) %>%
  mutate(Median_Wang_Phase_epithelial = as.factor(Median_Wang_Phase_epithelial))

#Plot and save
SFig3c <- ggplot(SupplementaryTable1, aes(x = Menstrual.cycle.day..Base.1., y = Progesterone.in.serum..nmol.L., color = FinalPhaseRefined, shape = endometriosis.status)) +
  geom_point() +
  geom_hline(yintercept=2, linetype="dashed", color = "dark grey") +
  coord_trans(y="log2")
ggsave(filename = paste0(dir_out, 'SFig3c.pdf'),
       plot = SFig3c,
       width =8 ,height=5)

## ---------------------------------- ##
##  Supplementary Data Figure 3.d
## ---------------------------------- ##

#Load data and clean
SupplementaryTable1 <- read.delim("../Tables/SupplementaryTable1.txt")
SupplementaryTable1 <- SupplementaryTable1 %>%
  filter(!(is.na(FinalPhaseRefined) | FinalPhaseRefined == ""))

#Plot and save
SFig3d <- ggplot(SupplementaryTable1, aes(x = Menstrual.cycle.day..Base.1., y = Progesterone.in.serum..nmol.L., color = Median_Wang_Phase_epithelial, shape = endometriosis.status)) +
  geom_point() +
  geom_hline(yintercept=2, linetype="dashed", color = "dark grey") +
  coord_trans(y="log2")
ggsave(filename = paste0(dir_out, 'SFig3d.pdf'),
       plot = SFig3d,
       width =8 ,height=5)

## ---------------------------------- ##
##  Supplementary Data Figure 3.e
## ---------------------------------- ##

#Load data
epithelial_cells_monocle3 <- readRDS(paste0(dir_data, "epithelial_cells_monocle3.rds"))

#Menstual Cycle Phase Wang
SFig3e <- plot_cells(epithelial_cells_monocle3 , color_cells_by = "WangMenstrualCyclePhase",
                     label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE,
                     graph_label_size=3, cell_size = 1)
SFig3e
ggsave(filename = paste0(dir_out, 'SFig3e.pdf'),
       plot = SFig3e,
       width =20, height=9)

## ---------------------------------- ##
##  Supplementary Data Figure 3.f
## ---------------------------------- ##

#Load data
EndoAtlas <- readRDS("../_Data/EndoAtlas.rds")

#Plot and save
SFig3f <- DimPlot(EndoAtlas, reduction = "umap", group.by = "WangMenstrualCyclePhase", raster = FALSE)
SFig3f 
ggsave(filename = paste0(dir_out, 'SFig3f.pdf'),
       plot = SFig3f,
       width =20, height=9)


## ---------------------------------- ##
##  Supplementary Data Figure 3.g
## ---------------------------------- ##

####### Download data from Tan et al. 2022 (DOI: 10.1038/s41556-022-00961-5) through the command line
#cd endometriosis_endometrium_scRNA_atlas/_Data
#wget --no-check-certificate https://singlecell.jax.org/datasets/endometriosis/h5ad/endo-2022_epithelial.h5ad

dir_data <- "../_Data/"

####### Convert .h5ad files to .h5seurat files
library(SeuratDisk)
Convert(paste0(dir_data,"endo-2022_epithelial.h5ad"), dest = paste0(dir_data,"endo-2022_epithelial.h5seurat"), assay = "RNA", overwrite = TRUE)
Tan_epithelial <- LoadH5Seurat(paste0(dir_data,"endo-2022_epithelial.h5seurat"), assays = "RNA")

#Create scale.data slot for the RNA assay required for heatmapping
Tan_epithelial_features <- Tan_epithelial@assays[["RNA"]]@data@Dimnames[[1]]
Tan_epithelial <- ScaleData(object = Tan_epithelial, features = Tan_epithelial_features)

#Subset cells from the endometrium (exclude cells from endometriosis lesions)
Idents(Tan_epithelial) <- Tan_epithelial$sample_type_rename
Tan_epithelial <- subset(Tan_epithelial, idents = c("Ctrl", "EuE"))

# factor levels in better order (E10 and E11 don't have endometrium cells)
Tan_epithelial$PID <- factor(Tan_epithelial$PID, levels = c("C01", "C02", "C03", "E01", "E02", "E03", "E04", "E05", "E06", "E07", "E08", "E09"))

# Generate heatmap of the endometrium cells from Tan et al. 2022 with the same epithelial menstrual marker genes as in our manuscript Figure 2d, for comparison
feature_list_figure2 <- c("SFRP4", "LAMC2", "MMP11", "MMP7", "SPRY1", "TFPI2","STEAP4","ADAMTS8","CAPN6","SLC37A2","SCGB1D2", "MT1F", "MT1X","RIMKLB", "S100P", "DEFB1","G0S2","NNMT","GPX3", "PAEP")

heatmap <- DoHeatmap(
  object = Tan_epithelial,
  features = intersect(feature_list_figure2, Tan_epithelial_features),
  cells = NULL,
  group.by = "PID",
  assay = "RNA",
  group.bar = TRUE) + 
  scale_fill_gradientn(limits = c(-2, 2), colors = c("#FF00FF", "#000000", "#FFFF00"))
heatmap
ggsave(
  filename = paste0(dir_out, "SFig3g_Tan_epithelial_Fig2EpiMarkers.pdf"),
  plot = last_plot(),
  width = 12.3,
  height = 5)


#As an alternative control heatmap of menstrual marker genes from Wang et al. 2020 (DOI: 10.1038/s41591-020-1040-z)
feature_list_Wang <- c("PLAU", "MMP7", "THBS1", "CADM1", "NPAS3", "ATP1A1", "ANK3", "ALPL", "TRAK1", "SCGB1D2", "MT1F", "MT1X", "MT1E", "MT1G", "CXCL14", "MAOA", "DPP4", "NUPR1", "GPX3", "PAEP")

heatmap <- DoHeatmap(
  object = Tan_epithelial,
  features = intersect(feature_list_Wang, Tan_epithelial_features),
  cells = NULL,
  group.by = "PID",
  assay = "RNA",
  group.bar = TRUE) + 
  scale_fill_gradientn(limits = c(-2, 2), colors = c("#FF00FF", "#000000", "#FFFF00"))
heatmap
ggsave(
  filename = paste0(dir_out, "SFig3g_Tan_epithelial_WangEpiMarkers.pdf"),
  plot = last_plot(),
  width = 15,
  height = 15)

#The menstrual marker gene expression varies largely between the endometrium samples in the Tan et al. 2022 dataset, strongly indicating not homogeneous menstrual cycle phase between the endometrium samples. 
#The expression distribution of menstrual marker genes from our manuscript figure 2d and Wang et al. 2020 look comparable.


## ---------------------------------- ##
##  Supplementary Data Figure 3.h
## ---------------------------------- ##

####### Download data from Tan et al. 2022 (DOI: 10.1038/s41556-022-00961-5) through the command line
#cd endometriosis_endometrium_scRNA_atlas/_Data
#wget --no-check-certificate https://singlecell.jax.org/datasets/endometriosis/h5ad/endo-2022_epithelial.h5ad

####### Convert .h5ad files to .h5seurat files
library(SeuratDisk)
Convert(paste0(dir_data,"endo-2022_stromal.h5ad"), dest = paste0(dir_data,"endo-2022_stromal.h5seurat"), assay = "RNA", overwrite = TRUE)
Tan_stromal <- LoadH5Seurat(paste0(dir_data,"endo-2022_stromal.h5seurat"), assays = "RNA")

#Create scale.data slot for the RNA assay required for heatmapping
Tan_stromal_features <- Tan_stromal@assays[["RNA"]]@data@Dimnames[[1]]
Tan_stromal <- ScaleData(object = Tan_stromal, features = Tan_stromal_features)

#Subset cells from the endometrium (exclude cells from endometriosis lesions)
Idents(Tan_stromal) <- Tan_stromal$sample_type_rename
Tan_stromal <- subset(Tan_stromal, idents = c("Ctrl", "EuE"))

# factor levels in better order (E10 and E11 don't have endometrium cells)
Tan_stromal$PID <- factor(Tan_stromal$PID, levels = c("C01", "C02", "C03", "E01", "E02", "E03", "E04", "E05", "E06", "E07", "E08", "E09"))

# Generate heatmap of the endometrium cells from Tan et al. 2022 with the same epithelial menstrual marker genes as in our manuscript Figure 2d, for comparison
feature_list_figure2_stromal <- c("SFRP1", "H1-1", "H3C2", "WNT5A","CILP","BRINP1", "CFD","PTGDS","SLIT3","RPRM", "TIMP3", "IGFBP1","LEFTY2","LRRC15","LUM")

heatmap <- DoHeatmap(
  object = Tan_stromal,
  features = intersect(feature_list_figure2_stromal, Tan_stromal_features),
  cells = NULL,
  group.by = "PID",
  assay = "RNA",
  group.bar = TRUE) + 
  scale_fill_gradientn(limits = c(-2, 2), colors = c("#FF00FF", "#000000", "#FFFF00"))
heatmap
ggsave(
  filename = paste0(dir_out, "SFig3h_Tan_stromal_Fig2EpiMarkers.pdf"),
  plot = last_plot(),
  width = 12.3,
  height = 5)


#As an alternative control heatmap of menstrual marker genes from Wang et al. 2020 (DOI: 10.1038/s41591-020-1040-z)
feature_list_Wang_stromal <- c('STC1', 'NFATC2', 'BMP2' ,'PMAIP1' ,'MMP11' ,'SFRP1', 'WNT5A' ,'ZFYVE21' ,'CILP' ,'SLF2' ,'MATN2', 'S100A4' ,'DKK1' ,'CRYAB', 'FOXO1', 'IL15', 'FGF7', 'LMCD1')

heatmap <- DoHeatmap(
  object = Tan_stromal,
  features = intersect(feature_list_Wang_stromal, Tan_stromal_features),
  cells = NULL,
  group.by = "PID",
  assay = "RNA",
  group.bar = TRUE) + 
  scale_fill_gradientn(limits = c(-2, 2), colors = c("#FF00FF", "#000000", "#FFFF00"))
heatmap
ggsave(
  filename = paste0(dir_out, "SFig3h_Tan_stromal_WangEpiMarkers.pdf"),
  plot = last_plot(),
  width = 15,
  height = 15)

#The menstrual marker gene expression varies largely between the endometrium samples in the Tan et al. 2022 dataset, strongly indicating not homogeneous menstrual cycle phases between the endometrium samples. 
#The expression distribution of menstrual marker genes from our manuscript figure 2d and Wang et al. 2020 look comparable.
