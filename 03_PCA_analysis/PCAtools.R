
###############################################
### Load Libraries and set path links
###############################################

#Load libraries
library(PCAtools)
library(Seurat)


##Set paths
dir_ENDO <- "/home/duempelmann/analysis-projects/ENDO/E044/A085/endometriosis_endometrium_scRNA_atlas/"
dir_out <- paste0(dir_ENDO,"_Data/03_PCA_analysis/")
dir.create(dir_out, recursive = TRUE)


###############################################
### Load Data
###############################################
##Load Data
EndoAtlas <- readRDS(paste0(dir_ENDO, "_Data/EndoAtlas.rds"))


###############################################
### Entire EndoAtlas: Data preparation and Conduct PCA
###############################################

##Data preparation
#Aggregate Expression per sample 
bulk <-  AggregateExpression(EndoAtlas, group.by = "sample", assay = "SCT", return.seurat = TRUE)
# Extract counts matrix
mat <- bulk[["SCT"]]$data


##Modify Metadata to match mat
EndoAtlas_meta <- EndoAtlas@meta.data  %>%
  distinct()
# Reorder EndoAtlas_meta to match the order of columns in mat
EndoAtlas_meta <- EndoAtlas_meta[match(colnames(mat), EndoAtlas_meta$sample), ]
# Set the row names of EndoAtlas_meta to be the reordered sample column
rownames(EndoAtlas_meta) <- EndoAtlas_meta$sample
## check that sample names match exactly between pdata and expression data 
all(colnames(mat) == rownames(EndoAtlas_meta))
## [1] TRUE

##Conduct principal component analysis (PCA):
p <- pca(mat, metadata = EndoAtlas_meta, removeVar = 0.1)

##save pca file
saveRDS(p, paste0(dir_out,"p_EndoAtlas.rds"))

##Remove variables
rm(bulk, mat, EndoAtlas_meta, p)



###############################################
### DEG samples: Data preparation and Conduct PCA
###############################################

##Data preparation
#Subset DEG samples
Idents(EndoAtlas) <- EndoAtlas$DEG.analysis..Figure.3a.b.
EndoAtlas <- subset(EndoAtlas, idents = "TRUE")

#Aggregate Expression per sample 
bulk <-  AggregateExpression(EndoAtlas, group.by = "sample", assay = "SCT", return.seurat = TRUE)
# Extract counts matrix
mat <- bulk[["SCT"]]$data


##Modify Metadata to match mat
EndoAtlas_meta <- EndoAtlas@meta.data  %>%
  distinct()
# Reorder EndoAtlas_meta to match the order of columns in mat
EndoAtlas_meta <- EndoAtlas_meta[match(colnames(mat), EndoAtlas_meta$sample), ]
# Set the row names of EndoAtlas_meta to be the reordered sample column
rownames(EndoAtlas_meta) <- EndoAtlas_meta$sample
## check that sample names match exactly between pdata and expression data 
all(colnames(mat) == rownames(EndoAtlas_meta))
## [1] TRUE

##Conduct principal component analysis (PCA):
p <- pca(mat, metadata = EndoAtlas_meta, removeVar = 0.1)

##save pca file
saveRDS(p, paste0(dir_out,"p_EndoAtlas_DEGsamples.rds"))

##Remove variables
rm(bulk, mat, EndoAtlas_meta, p)



###############################################
### Main cell types: Data preparation and Conduct PCA
###############################################

# Set Idents of EndoAtlas to AnnotationMain
Idents(EndoAtlas) <- EndoAtlas$AnnotationMain

#SUBSET <- "epithelial"
subsets <- c("myeloid", "endothelial", "epithelial", "lymphocyte","stromal")

#Function
process_PCA <- function(SUBSET) {
  print(SUBSET)
  # subset with the main 5 cell types
  subset <- subset(EndoAtlas,idents = SUBSET)
  
  # pseudobulk cells only by sample
  bulk <- AggregateExpression(subset, group.by = "sample", assay = "SCT", return.seurat = TRUE)
    # Extract counts matrix
  mat <- bulk[["SCT"]]$data 

  ##Modify Metadata to match mat
  subset_meta <- subset@meta.data  %>%
    distinct()
  # Reorder EndoAtlas_meta to match the order of columns in mat
  subset_meta <- subset_meta[match(colnames(mat), subset_meta$sample), ]
  # Set the row names of EndoAtlas_meta to be the reordered sample column
  rownames(subset_meta) <- subset_meta$sample
  

  ## check that sample names match exactly between pdata and expression data 
  all(colnames(mat) == rownames(subset_meta))
  ## [1] TRUE
  
  ##Conduct principal component analysis (PCA):
  p <- pca(mat, metadata = subset_meta, removeVar = 0.1)

  ##save pca file
  saveRDS(p, paste0(dir_out,"p_EndoAtlas_",SUBSET,".rds"))
  
  gc()
  
}

lapply(subsets, process_PCA)

rm(data)
gc()


###############################################
### Save sessionInfo
###############################################
writeLines(capture.output(sessionInfo()), paste0(dir_out,"sessionInfo.txt"))
