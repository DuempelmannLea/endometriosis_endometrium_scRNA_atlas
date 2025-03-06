
####################
### 00. Load libs and set paths
####################
library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)


dir_out <- "/endometriosis_endometrium_scRNA_atlas/_Data/07_Ligand_receptor_analysis"
dir.create(paste0(dir_out,"output/"), recursive = TRUE)
setwd(dir_out)

path_seurat_object <- "/endometriosis_endometrium_scRNA_atlas/_Data/ENDO_global.rds"
dir_out <- "CellChat/"
setwd(dir_out)

####################
### 01. Prepare Seurat for CellChat
####################

##Load seurat object and add metadata
data  <- readRDS(path_seurat_object)

##Filter Seurat for proliferative phase samples with stringet exclusion criteria and annotated cells
DEGsamples <- unique(data@meta.data %>% filter(DEG.analysis..Figure.3a.b. == "TRUE") %>% pull(sample))
Idents(data) <- data$sample
data <- subset(data, idents = DEGsamples) 
data <- subset(data, subset = AnnotationRefined %in% levels(as.factor(data$AnnotationRefined)))

#Filter out secretory specific cell types (no/low cell counts stop script 02)
secretory_celltypes <- c("dS3","dS4","fib C7 secretory","glandular early-secretory","mid-secretory","Prv_VSMC secretory","Prv_VSMC secretory_late","fib C7 secretory_late") 
data <- subset(data, subset = AnnotationRefined %in% setdiff(levels(as.factor(data$AnnotationRefined)), secretory_celltypes))

# Define the conditions for Case and Control
Endo <- data@meta.data$sample %in% c(filter(data@meta.data, EndometriosisStatus == "Endo") %>% distinct(sample) %>% pull())
Control <- data@meta.data$sample %in% c(filter(data@meta.data, EndometriosisStatus == "CTL") %>% distinct(sample) %>% pull())

# Create the Condition column based on Endo/Control
data@meta.data$Condition <- NA
data@meta.data$Condition[Endo] <- "Endo"
data@meta.data$Condition[Control] <- "Control"
DimPlot(data, group.by = "Condition", raster = FALSE)


####################
### 02. Run CellChat
####################

runCellChat <- function(CONDITION){
      #Create cellchat 
      seurat_sens <- subset(x = data, subset = (Condition == CONDITION))
      DimPlot(seurat_sens, group.by = "Condition")
      data.input <- GetAssayData(seurat_sens, assay = "SCT", slot = "data") # normalized data matrix
      Idents(seurat_sens) <- seurat_sens@meta.data$AnnotationRefined
      labels <- Idents(seurat_sens)
      meta <- data.frame(labels = labels, row.names = names(labels))
      meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
      table(meta$labels)
      cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
      
      # load a scRNA-seq data matrix and its associated cell meta data
      CellChatDB <- CellChatDB.human 
      showDatabaseCategory(CellChatDB)
      CellChatDB.use <- CellChatDB # use the default CellChatDB
      
      # set the used database in the object
      cellchat@DB <- CellChatDB.use
      
      # subset the expression data of signaling genes for saving computation cost
      cellchat <- subsetData(cellchat) 
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat)
      
      # project gene expression data onto PPI
      cellchat <- projectData(cellchat, PPI.human)
      
      ptm = Sys.time()
      cellchat <- computeCommunProb(cellchat, type = "triMean")
      cellchat <- filterCommunication(cellchat, min.cells = 10)
      cellchat <- computeCommunProbPathway(cellchat)
      cellchat <- aggregateNet(cellchat)
      
      groupSize <- as.numeric(table(cellchat@idents))
      pdf(paste0(dir_out,"output/",CONDITION,"_NumberInteractions_InteractionWeights.pdf"))
      par(mfrow = c(1,2), xpd=TRUE)
      netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
      netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
      dev.off()
      
      saveRDS(cellchat, file = paste0(dir_out,CONDITION,"_CellChat_AnnotationRefined.rds"))
}

runCellChat("Endo")
runCellChat("Control")

rm(data)
gc()

####################
### 03. General visualization
####################

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
      dir.create(paste0(dir_out,"Endo_netVisual/"))
      pathways_list <- cellchat@netP$pathway
      for (PATHWAY in pathways_list) {
        pathways.show <- PATHWAY
        netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
        pdf(file = paste0(dir_out,"output_comparison_Refined/Endo_netVisual/",CONDITION,"_",PATHWAY,"_circle_plot.pdf"))
        netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
        dev.off() 
      }
}

CellChat_visualization("Endo")
CellChat_visualization("Control")
