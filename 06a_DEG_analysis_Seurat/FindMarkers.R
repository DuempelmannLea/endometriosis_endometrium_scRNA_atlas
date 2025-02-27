## ------------------------------------------ ##
## Setup
## ------------------------------------------ ##

# load packages 
library(dplyr)
library(tibble)
library(SeuratWrappers)
library(Seurat)

# Define paths 
output_path <- ''
input_seurat_proliferative_object_path <- ''


## ------------------------------------------ ##
## Process Seurat object 
## ------------------------------------------ ##

# Load seurat object 
data <- readRDS(input_seurat_proliferative_object_path)

# Identify lowly expressed genes and remove them from the seurat object 
genes_keep <- rowSums(data@assays$RNA@counts > 1) >= 10
seurat_object_filtered <- data[names(which(genes_keep)),]

## ------------------------------------------ ##
## Perform find markers 
## ------------------------------------------ ##

##0. Functions
clean_ct_naming <- function(cell_types_vector){
  cell_types_vector <- gsub('M$\\Phi$','Macrophages_PHI_',cell_types_vector,fixed=T)
  cell_types_vector <- gsub('[\\$\\{\\}\\/]','',cell_types_vector)
}
perform_FindMarkers_per_cell_type_wbatch_correction <- function(cell_type,seurat_object,output_path,cell_type_annotation, gene_filtering){

  if (gene_filtering){
    output_file <- file.path(output_path,paste0('DEG_CTL_vs_ENDO_',cell_type,'_',cell_type_annotation,'_with_batch_correction_with_gene_filtering.csv'))
  } else{
    output_file <- file.path(output_path,paste0('DEG_CTL_vs_ENDO_',cell_type,'_',cell_type_annotation,'_with_batch_correction.csv'))
  }
  
  if (!file.exists(output_file)) {
    tryCatch({
      print(paste0('Running DEG analysis for the following cell-type: ', cell_type))
      
      if (gene_filtering){
        Idents(seurat_object) <- seurat_object$cell_type_annotation 
        seurat_object_subset <- subset(seurat_object,cell_type_annotation == cell_type)
        # Added this additional filtering to fit with the default  muscat function
        genes_keep <- rowSums(seurat_object_subset@assays$RNA@counts > 1) >= 20
        seurat_object_filtered <- seurat_object_subset[names(which(genes_keep)),]
        DefaultAssay(seurat_object_filtered) <- 'RNA'
        deg_results <- FindMarkers(seurat_object_filtered,
                                   ident.1 = 'CTL',
                                   ident.2 = 'Endo',
                                   group.by = 'treatment',
                                   subset.ident = cell_type, 
                                   min.pct = 0, 
                                   test.use = 'LR',
                                   latent.vars=c('batch'),
                                   assay='RNA')
      } else{
        Idents(seurat_object) <- seurat_object$cell_type_annotation
        seurat_object_subset <- subset(seurat_object,cell_type_annotation == cell_type)
        DefaultAssay(seurat_object_subset) <- 'RNA'
        deg_results <- FindMarkers(seurat_object_subset,
                                   ident.1 = 'CTL',
                                   ident.2 = 'Endo',
                                   group.by = 'treatment',
                                   subset.ident = cell_type,
                                   test.use = 'LR',
                                   latent.vars=c('batch'),
                                   min.pct = 0)
      }
      
      

      deg_results <- deg_results  %>% rownames_to_column('gene')
      write.csv(deg_results,file = output_file,
                quote=F,row.names = F)
      print("DEG analysis completed successfully.")
    }, error = function(e) {
      print(paste("Error in DEG analysis:", e$message))
    })
  } else {
    print("File already exists. Skipping analysis.")
  }
}
perform_FindMarkers_per_cyclephase <- function(
                                         seurat_object,
                                         output_path,cell_type_column){
  
  # clean up cell-types names 
  seurat_object$cell_type_annotation <- clean_ct_naming(as.vector(unlist(seurat_object[[cell_type_column]])))
  cell_types <- unique(seurat_object$cell_type_annotation)
  
  # Run find markers DEGS with correction and with filtering
  lapply(cell_types,
         perform_FindMarkers_per_cell_type_wbatch_correction,
         seurat_object,
         output_path,cell_type_annotation=cell_type_column,gene_filtering=TRUE)
}
##1. Run find markers analysis for proliferative samples 
cell_annotation_variables <- c("AnnotationMain","AnnotationRefined","AnnotationUnited" )
for (var in cell_annotation_variables) {
  perform_FindMarkers_per_cyclephase(seurat_object = seurat_object_filtered,
                                     output_path = output_path,
                                     cell_type_column = var)
  
}
