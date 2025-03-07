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
##### 01. Run Symphony on entire EndoAtlas
#######################################################

#create first following symbolic link in Linux environment:
#ln -s ../_Data/EndoAtlas.rds ../_Data/02_Cell_annotation/ENDO_global.rds

#00 global annotation
symphony(
  CELLTYPEref = "Tan_global",
  CELLTYPEquery = "ENDO_global"
)


########################################################
##### 02. Save EndoAtlas split into the 5 main cell types
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
##### 03. Run Symphony on 5 main cell types separately for refined annotation
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
