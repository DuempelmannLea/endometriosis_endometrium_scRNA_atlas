####################################
## 00. Load library and environment
####################################


suppressPackageStartupMessages({
library(ggplot2)
library(dplyr)
library(limma)
library(muscat)
library(purrr)
library(Seurat)
library(scater)
library(SingleCellExperiment)
library(Matrix)
})

####################################
## 01. Prepare Data
####################################

#load seurat object 
data <- readRDS(snakemake@input[[1]])
print("loaded data")

#Filter for proliferative samples with stringent exclusion criteria
Idents(data) <- data$DEG.analysis.Figure.3a.b
data <- subset(data, idents = "TRUE")



# List of annotations to loop through
ANNOTATIONS <- c("AnnotationMain", "AnnotationRefined", "AnnotationUnited")

# Loop through each annotation
for (Annotation in ANNOTATIONS) {

	#make SingleCellExperiment from RNA assay
	sce <- as.SingleCellExperiment(data, assay = "RNA")

	# remove undetected genes
	sce <- sce[rowSums(counts(sce) > 0) > 0, ]

	# remove lowly expressed genes (less than 10 counts)
	sce <- sce[rowSums(counts(sce) > 1) >= 10, ]

	##Data preparation
	sce$id <- paste0(sce$treatment, sce$sample)

	(sce <- prepSCE(sce,
					kid = Annotation, # subpopulation assignments
					gid = "condition",  # group IDs (ENDO/non-ENDO)
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

	saveRDS(mm_dream, paste0("./mm_dream",Annotation,".csv"))
}

sessionInfo()
date()
