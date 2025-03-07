####################################
## 00. Set path
####################################

##Set paths
dir_out <- paste0("../_Data/05a_DEG_analysis_muscat/")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)


DEGanalysis_muscat <- function(MenstrualCyclePhase, Annotation) {

####################################
## 01. Prepare Data
####################################

##load seurat object 
data <- readRDS("../_Data/EndoAtlas.rds")
print("loaded data")

##Subset samples by menstrual cycle phase
Idents(data) <- data$sample

# Proliferative phase
if (MenstrualCyclePhase == "Proliferative") {
  proliferativeDEG_samples <- data$DEG.analysis.Figure.3a.b
  # Alternatively, manually specify the samples
  # proliferativeDEG_samples <- c("sample76", "sample83", "sample84", "sample91", "sample94", 
  #                               "sample106", "sample111", "sample125", "sample126", "sample130", 
  #                               "sample134", "sample135", "sample70", "sample71", "sample73", 
  #                               "sample82", "sample95", "sample96", "sample110", "sample112", 
  #                               "sample114", "sample116", "sample123")
  
  data <- subset(data, idents = proliferativeDEG_samples)
}

# Secretory phase
if (MenstrualCyclePhase == "Secretory") {
  #Subset of balanced secretory phase samples in terms of EndometriosisStatus, CycleDay and MenstrualCyclePhase_Refined
  secretoryDEG_samples <- c("sample133", "sample113", "sample117", "sample88", 
                            "sample128", "sample132", "sample127", "sample75")
  data <- subset(data, idents = secretoryDEG_samples)
}

##Rename required Annotation to 'Annotation'
if (Annotation == "AnnotationMain") {
  data$Annotation <- data$AnnotationMain
}
if (Annotation == "AnnotationRefined") {
  data$Annotation <- data$AnnotationRefined
}
if (Annotation == "AnnotationUnited") {
  data$Annotation <- data$AnnotationUnited
}

##Reduce metadata
data@meta.data <- data@meta.data[,c("sample", "EndometriosisStatus", "Annotation","batch")]

##make SingleCellExperiment from RNA assay
sce <- as.SingleCellExperiment(data, assay = "RNA")

# remove undetected genes
sce <- sce[rowSums(counts(sce) > 0) > 0, ]

# remove lowly expressed genes (less than 10 counts)
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]

##Data preparation
sce$id <- paste0(sce$treatment, sce$sample)

(sce <- prepSCE(sce,
                kid = "Annotation", # subpopulation assignments
                gid = "EndometriosisStatus",  # group IDs (ENDO/non-ENDO)
                sid = "sample",   # sample IDs 
                drop = FALSE))  

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

####################################
## 02. Cell-level analysis: Mixed models
####################################

print("start DE testing:")
mm_dream <- mmDS(sce, method = "dream", coef = "group_idEndo",
                 n_cells = 50, n_samples = 3, covs = "batch")

saveRDS(mm_dream, paste0(dir_out, MenstrualCyclePhase,"_", Annotation,"/mm_dream.rds"))

tbl_fil <- plyr::rbind.fill(mm_dream) %>%
  dplyr::filter(p_adj.loc < 0.05) %>%
  dplyr::arrange(p_adj.loc)

write.csv(tbl_fil, paste0(dir_out, MenstrualCyclePhase,"_", Annotation,"/mm_dream.csv"))

sessionInfo()
date()

}






