####################################
## 00. Load library and environment
####################################

# log
log <- file(snakemake@log$log_file, open="wt")
sink(log)
sink(log, type="message")

# image
save.image(file = paste0(snakemake@log$log_object) )

# load packages
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(future)
  library(RColorBrewer)
  library(remotes)
  library(ggplot2)
  library(SingleCellExperiment)
})

# dir_out
dir_out <- snakemake@config$integration_output_dir


# load seurat object
seurat_obj_split_input <- readRDS(snakemake@input[[1]])
print("input loaded")



####################################
## 01. Run integration
####################################

rpca_integration = function(seurat_object_list, snakemake, features_integration, verbose = TRUE){

  # RPCA integration on SCT transformed counts
  features_integration <- SelectIntegrationFeatures(object.list = seurat_obj_split_input,
                                                    nfeatures = 3000)
  seurat_object_list <- Seurat::PrepSCTIntegration(object.list = seurat_object_list,
                                                   assay = "SCT",
                                                   anchor.features = features_integration,
                                                   verbose = verbose
  )
  seurat_object_list <- lapply(seurat_object_list,
                               FUN = Seurat::RunPCA,
                               anchor.features = features_integration
  )


  reference_dataset <- which(names(seurat_object_list) == snakemake@config$reference_dataset)


  integration_anchors <- FindIntegrationAnchors(object.list = seurat_object_list,
                                                normalization.method = "SCT",
                                                anchor.features = features_integration,
                                                dims = seq_len(snakemake@config$n_dimred_integration),
                                                reduction = "rpca",
                                                k.anchor = 20,
                                                reference = reference_dataset,
                                                verbose = verbose
  )
  print(paste0('Length of features used for integration',
               length(integration_anchors@anchor.features)))
  print(paste0('Number of cells used for integration ',
               dim(integration_anchors@anchors)[1]))
  seurat_object_list <- IntegrateData(anchorset = integration_anchors,
                                      normalization.method = "SCT",
                                      dims = seq_len(snakemake@config$n_dimred_integration),
                                      verbose = verbose
  )
  return(seurat_object_list )
}

# Run RPCA
integrated_data <- rpca_integration(seurat_obj_split_input, snakemake, features_integration)
DefaultAssay(integrated_data) <- "integrated"


seed <- 122334

# Perform linear dimensional reduction by PCA on HVGs
integrated_data <- RunPCA(object = integrated_data,
                     features = VariableFeatures(integrated_data),
                     seed.use = seed,
                     npcs = 50)
# Run UMAP
integrated_data <- RunUMAP(integrated_data,  dims = 1:30, seed.use = seed)

# Save the filtered, normalized and integrated data
saveRDS(integrated_data, snakemake@output[[1]])
print("saved integrated_data.rds")


sessionInfo()
date()
