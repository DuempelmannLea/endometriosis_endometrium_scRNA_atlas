########################
#### 00. Load libraries and set paths
#######################
library(Seurat)
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggbreak) 
library(patchwork)
library(ComplexUpset)
library(ggrepel)
library(PCAtools)
library(scCustomize)
dir_out <- "ENDO/E044/A074/"
dir.create(dir_out, recursive = TRUE)
setwd(dir_out)

########################
#### 01. Load data and plot
#######################
###Load data
data <- readRDS("./ENDOatlas.rds")

###Plotting for overview
DimPlot_scCustom(data, reduction = "umap", group.by = "MenstrualCyclePhase")
DimPlot_scCustom(data, reduction = "umap", group.by = "Symphony_Global_main5")
DimPlot_scCustom(data, reduction = "umap", group.by = "Symphony_Refined")
DimPlot_scCustom(data, reduction = "umap", group.by = "sample", split.by = "MenstrualCyclePhase")
#Highlight samples 65 and 130 (they are in the transition between periovulatory and secretory menstual cycle phase)
Idents(data) <- data$sample
cellshighlight <- WhichCells(data, idents = c("sample65","sample130"))
DimPlot(data, label=T,reduction = "umap", group.by = "treatment", cells.highlight = cellshighlight, cols.highlight = c("darkred"), cols= "grey", raster = FALSE)

############################################################################
############################################################################
########################
#### 02. Clustering of entire atlas
#######################
## Make clustering with different resolutions
data <- FindNeighbors(data, dims = 1:30)
data <- Seurat::FindClusters(data, algorithm = 2 ,resolution = 0.1) #random.seed is set by default to 0 #algorithm = 2 gave better results than algorithm = 1
data <-Seurat::FindClusters(data, algorithm = 2 ,resolution = 0.2) 
data <-Seurat::FindClusters(data, algorithm = 2 ,resolution = 0.5) 
data <-Seurat::FindClusters(data, algorithm = 2 ,resolution = 1)

#Plot clustering with different resolutions and compare
DimPlot_scCustom(data, reduction = "umap", split.by = "integrated_snn_res.0.1",group.by = "MenstrualCyclePhase",  raster=FALSE)
DimPlot_scCustom(data, reduction = "umap", split.by = "integrated_snn_res.0.2",group.by = "MenstrualCyclePhase",  raster=FALSE)
DimPlot_scCustom(data, reduction = "umap", split.by = "integrated_snn_res.0.5",group.by = "MenstrualCyclePhase",  raster=FALSE)
DimPlot_scCustom(data, reduction = "umap", split.by = "integrated_snn_res.1",group.by = "MenstrualCyclePhase",  raster=FALSE)


############################################################################
############################################################################
########################
#### 03. Menstrual cycle specific mesenchymal cells
#######################
Idents(data) <- data$Symphony_Global_main5
data <- subset(data, idents = "stromal")
Idents(data) <- data$integrated_snn_res.0.1
#Cluster 2 corresponds to secretory samples and clusters 7 + 14 to secretory_late samples for mesenchymal cells
################################################
##dS3 and fib C7 secretory mesenchymal cells 
################################################
cluster2 <- subset(data, idents = "2" )
table(cluster2$Symphony_Refined)
table(cluster2$MenstrualCyclePhase,cluster2$sample)
#The numbers fit very nicely with the progression of cycle stage, including that sample 65 and 105 have higher numbers than the other periovulatory samples
DimPlot(cluster2, group.by = "Symphony_Refined")
DimPlot(cluster2, group.by = "Symphony_Refined", split.by = "MenstrualCyclePhase")
DimPlot(cluster2, group.by = "MenstrualCyclePhase", split.by = "Symphony_Refined")

#Cluster2 "dS1-myofibroblast","dS2","eF1","eF2","eF3" become dS3 (decidualized stromal 3)
meta_final <- mutate(meta_final,
                     Symphony_Refined_Final = case_when(
                       integrated_snn_res.0.1 == 2 & Symphony_Refined_Final %in% c("dS1-myofibroblast", "dS2", "eF1", "eF2", "eF3") ~ "dS3",
                       TRUE ~ Symphony_Refined_Final
                     ))
table(meta_final$Symphony_Refined_Final) #sanity check

#Cluster2 "fib C7","fib C7-SFRP2" become "fib C7 secretory"
meta_final <- mutate(meta_final,
                     Symphony_Refined_Final = case_when(
                       integrated_snn_res.0.1 == 2 & Symphony_Refined_Final %in% c("fib C7","fib C7-SFRP2") ~ "fib C7 secretory",
                       TRUE ~ Symphony_Refined_Final
                     ))
table(meta_final$Symphony_Refined_Final) #sanity check


################################################
##dS4, fib C7 secretory_late and Prv_VSMC secretory_late mesenchymal cells 
################################################
Idents(data) <- data$integrated_snn_res.0.1
cluster7_14 <- subset(data, idents = c("7","14"))
table(cluster7_14$MenstrualCyclePhase,cluster7_14$sample)
table(cluster7_14$Symphony_Refined) #corresponds mainly to eF3 (7067), dS2 (1550) and dS1 (1737)
#Fits nicely, samples 104,131 and 87 high numbers

DimPlot(cluster7_14, group.by = "Symphony_Refined")
DimPlot(cluster7_14, group.by = "Symphony_Refined", split.by = "MenstrualCyclePhase")
DimPlot(cluster7_14, group.by = "MenstrualCyclePhase", split.by = "Symphony_Refined")
table(data$Symphony_Refined)

table(meta_final$Symphony_Refined_Final)
#Cluster7_14 "dS1-myofibroblast","dS2","eF1","eF2","eF3" become dS4
meta_final <- mutate(meta_final,
                     Symphony_Refined_Final = case_when(
                       integrated_snn_res.0.1 %in% c(7, 14) & Symphony_Refined_Final %in% c("dS1-myofibroblast", "dS2", "eF1", "eF2", "eF3") ~ "dS4",
                       TRUE ~ Symphony_Refined_Final
                     ))
table(meta_final$Symphony_Refined_Final) #sanity check

#Cluster7_14 "fib C7","fib C7-SFRP2" become "fib C7 secretory_late"
meta_final <- mutate(meta_final,
                     Symphony_Refined_Final = case_when(
                       integrated_snn_res.0.1 %in% c(7, 14) & Symphony_Refined_Final %in% c("fib C7","fib C7-SFRP2") ~ "fib C7 secretory_late",
                       TRUE ~ Symphony_Refined_Final
                     ))
table(meta_final$Symphony_Refined_Final) #sanity check

#Cluster7_14 "Prv-CCL19","Prv-MYH11", "Prv-STEAP4" , "VSMC" become "Prv_VSMC secretory_late"
meta_final <- mutate(meta_final,
                     Symphony_Refined_Final = case_when(
                       integrated_snn_res.0.1 %in% c(7, 14) & Symphony_Refined_Final %in% c("Prv-CCL19","Prv-MYH11", "Prv-STEAP4" , "VSMC") ~ "Prv_VSMC secretory_late",
                       TRUE ~ Symphony_Refined_Final
                     ))
table(meta_final$Symphony_Refined_Final) #sanity check

table(meta_final$Symphony_Refined_Final,meta_final$integrated_snn_res.0.1)

################################################
##Prv and VSMC from secretory samples
################################################
Idents(data) <- data$integrated_snn_res.0.1
cluster8 <- subset(data, idents = "8")
table(cluster8$MenstrualCyclePhase,cluster8$sample)
table(cluster8$Symphony_Refined)
#The numbers fit very nicely with the progression of cycle stage, including that sample65, 105 have much more than the other periovulatory, and sample 99 and 103 have some more than other periovulatory samples
DimPlot(cluster8, group.by = "Symphony_Refined") 
DimPlot(cluster8, group.by = "MenstrualCyclePhase") +
  DimPlot(cluster8, group.by = "integrated_snn_res.1")
DimPlot(cluster8, group.by = "Symphony_Refined", split.by = "MenstrualCyclePhase")
DimPlot(cluster8, group.by = "MenstrualCyclePhase", split.by = "Symphony_Refined")

#From integrated_snn_res.1, cluster 8 corresponds best to the secretory samples, however not perfectly. After redoing the clustering there is great correspondence
cluster8 <- FindNeighbors(cluster8, dims = 1:30)
cluster8 <- Seurat::FindClusters(cluster8, algorithm = 2 ,resolution = 0.1) 
DimPlot_scCustom(cluster8, reduction = "umap", group.by = "integrated_snn_res.0.1") +
  DimPlot_scCustom(cluster8, reduction = "umap", group.by = "MenstrualCyclePhase")
DimPlot_scCustom(cluster8, reduction = "umap", group.by = "Symphony_Refined")
DimPlot_scCustom(cluster8, reduction = "umap", group.by = "MenstrualCyclePhase", split.by = "integrated_snn_res.0.1") #This clustering corresponds very well with the MenstrualCyclePhase
DimPlot_scCustom(cluster8, reduction = "umap", group.by = "sample", split.by = "MenstrualCyclePhase")

Idents(cluster8) <- cluster8$integrated_snn_res.0.1
cluster8_1 <- subset(cluster8, idents = "1")
table(cluster8_1$MenstrualCyclePhase,cluster8_1$sample)
table(cluster8_1$Symphony_Refined)

#Join 
meta_final <- left_join(meta_final, cluster8_1_integrated_snn_res.0.1, by = "cell_ident")

#Cluster8_1 "Prv-CCL19","Prv-MYH11", "Prv-STEAP4" , "VSMC" become "Prv_VSMC secretory"
meta_final <- mutate(meta_final,
                     Symphony_Refined_Final = case_when(
                       cluster8_1_integrated_snn_res.0.1 == 1 & Symphony_Refined_Final %in% c("Prv-CCL19","Prv-MYH11", "Prv-STEAP4" , "VSMC") ~ "Prv_VSMC secretory",
                       TRUE ~ Symphony_Refined_Final
                     ))
table(meta_final$Symphony_Refined_Final) #sanity check

############################################################################
############################################################################
########################
#### 04. Menstrual cycle specific epithelial cells
#######################
###Look into difference between glandular (proliferative), (glandular) mid-secretory, glandular early-secretory
data <- readRDS("./epithelial_reintegrated.rds")
data <- AddMetaData(data, meta_final)
DimPlot(data, split.by = "MenstrualCyclePhase", group.by = "Symphony_Refined_Final")
#there are quite some glandular early-secretory already in proliferative/periovulatory phase and a lot in secretory_early. lumenal 2 are mainly in lumenal 2. mid-secretory are in secretory_mid/late, this matches

########################
#### Prv_VSMC secretory
#######################
data <- FindNeighbors(data, dims = 1:30)
data <- Seurat::FindClusters(data, algorithm = 2 ,resolution = 0.05) 
DimPlot_scCustom(data, reduction = "umap", group.by = "integrated_snn_res.0.05") +
  DimPlot_scCustom(data, reduction = "umap", group.by = "MenstrualCyclePhase")

Idents(data) <- data$integrated_snn_res.0.05
epithelial_clusters_0_1 <- subset(data, idents = c("0","3"))
DimPlot_scCustom(epithelial_clusters_0_1, reduction = "umap", group.by = "integrated_snn_res.0.05")
table(epithelial_clusters_0_1$MenstrualCyclePhase,epithelial_clusters_0_1$sample)
table(epithelial_clusters_0_1$Symphony_Refined)

epithelial_clusters_0_1_res.0.05 <- as.data.frame(epithelial_clusters_0_1$integrated_snn_res.0.05)
epithelial_clusters_0_1_res.0.05 <- epithelial_clusters_0_1_res.0.05 %>%
  dplyr::rename(epithelial_clusters_0_1_res.0.05 = "epithelial_clusters_0_1$integrated_snn_res.0.05")
epithelial_clusters_0_1_res.0.05$cell_ident <- row.names(epithelial_clusters_0_1_res.0.05)

#Join with meta_final
meta_final <- left_join(meta_final, epithelial_clusters_0_1_res.0.05, by = "cell_ident")

#rename of cluster7_14 "Prv-CCL19","Prv-MYH11", "Prv-STEAP4" , "VSMC" to "Prv_VSMC secretory"
meta_final <- mutate(meta_final,
                     Symphony_Refined_Final = case_when(
                       epithelial_clusters_0_1_res.0.05 %in% c("0","3") & Symphony_Refined_Final %in% c("glandular early-secretory","mid-secretory") ~ "glandular",
                       TRUE ~ Symphony_Refined_Final
                     ))
table(meta_final$Symphony_Refined_Final) #sanity check

############################################################################
############################################################################
########################
#### 05. Save Final annotation
#######################
saveRDS(meta_final, paste0(dir_out,"Annotation_Final.rds"))

