########################################################
##### 00. Load packages and utils
#######################################################

# load packages
suppressPackageStartupMessages({
  source('../02_Cell_annotation/01_libs.R')# imports
}) 

#source utils
source('../02_Cell_annotation/01_utils.R') # color definitions, symphony function and plotting functions

########################################################
##### 01. Download data from Tan et al. 2022 (DOI: 10.1038/s41556-022-00961-5) as annotation reference
#######################################################

####### Download data in the command line
#cd ../_Data/02_Cell_annotation/
#wget --no-check-certificate https://singlecell.jax.org/datasets/endometriosis/h5ad/endo-2022_global.h5ad
#wget --no-check-certificate https://singlecell.jax.org/datasets/endometriosis/h5ad/endo-2022_endothelial.h5ad
#wget --no-check-certificate https://singlecell.jax.org/datasets/endometriosis/h5ad/endo-2022_epithelial.h5ad
#wget --no-check-certificate https://singlecell.jax.org/datasets/endometriosis/h5ad/endo-2022_lymphocyte.h5ad
#wget --no-check-certificate https://singlecell.jax.org/datasets/endometriosis/h5ad/endo-2022_myeloid.h5ad
#wget --no-check-certificate https://singlecell.jax.org/datasets/endometriosis/h5ad/endo-2022_stromal.h5ad

####### Convert .h5ad files to .h5seurat files and save as .rds
library(SeuratDisk)

# List of cell types
cell_types <- c("global", "lymphocyte", "epithelial", "endothelial", "myeloid", "stromal")

# Convert each file to h5seurat format and save as RDS
for (cell_type in cell_types) {
  # Convert to h5seurat
  Convert(paste0("endo-2022_", cell_type, ".h5ad"), dest = "h5seurat", assay = "RNA", overwrite = TRUE)
  
  # Load data and save as RDS
  data <- LoadH5Seurat(paste0("endo-2022_", cell_type, ".h5seurat"), assays = "RNA")
  saveRDS(data, paste0("endo-2022_", cell_type, ".rds"))
}


########################################################
##### 02. Run Symphony on entire EndoAtlas
#######################################################

#create first following symbolic link in Linux environment:
#ln -s ../_Data/EndoAtlas.rds ../_Data/02_Cell_annotation/ENDO_global.rds

#00 global annotation
symphony(
  CELLTYPEref = "Tan_global",
  CELLTYPEquery = "ENDO_global"
)


########################################################
##### 03. Save EndoAtlas split into the 5 main cell types
#######################################################

#Load integrated EndoAtlas
data <- readRDS("/endometriosis_endometrium_scRNA_atlas/_Data/EndoAtlas.rds")
Idents(data) <- data$AnnotationMain

# Define the vector of annotations
ANNOTATIONMAIN <- c("stromal", "epithelial", "endothelial", "lymphocyte", "myeloid")
# Loop over each annotation type
for (annotation in ANNOTATIONMAIN) {
  # Subset data based on the current annotation
  AnnotationMain_subset <- subset(data, idents = annotation)
  # Save the subsetted data to an RDS file
  saveRDS(AnnotationMain_subset, paste0("../_Data/02_Cell_annotation/ENDO_",annotation,".rds"))
}


########################################################
##### 04. Run Symphony on 5 main cell types separately for refined annotation
#######################################################


#01 lymphocyte
symphony(
  CELLTYPEref = "Tan_lymphocyte",
  CELLTYPEquery = "ENDO_lymphocyte"
)

#02 myeloid
symphony(
  CELLTYPEref = "Tan_myeloid",
  CELLTYPEquery = "ENDO_myeloid"
)

#03 epithelial
symphony(
  CELLTYPEref = "Tan_epithelial",
  CELLTYPEquery = "ENDO_epithelial"
)

#04 endothelial
symphony(
  CELLTYPEref = "Tan_endothelial",
  CELLTYPEquery = "ENDO_endothelial"
)

#05 stromal fibroblasts
symphony(
  CELLTYPEref = "Tan_stromal",
  CELLTYPEquery = "ENDO_stromal"
)
