####################################
## 00. Load library and environment
####################################

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

# load seurat object
seurat_obj_split_input <- readRDS("./seurat_obj_split.rds")
print("input loaded")

# dir_out
dir_out <- snakemake@config$integration_output_dir

####################################
## 01. Run integration
####################################

# RPCA integration on SCT transformed counts
features_integration <- SelectIntegrationFeatures(object.list = seurat_object_list,
                                                  nfeatures = 3000)
seurat_object_list <- Seurat::PrepSCTIntegration(object.list = seurat_object_list,
                                                 assay = "SCT",
                                                 anchor.features = features_integration,
                                                 verbose = TRUE
)
seurat_object_list <- lapply(seurat_object_list,
                             FUN = Seurat::RunPCA,
                             anchor.features = features_integration
)


reference_dataset <- which(names(seurat_object_list) == "batch17")


integration_anchors <- FindIntegrationAnchors(object.list = seurat_object_list,
                                              normalization.method = "SCT",
                                              anchor.features = features_integration,
                                              dims = seq_len(20),
                                              reduction = "rpca",
                                              k.anchor = 20,
                                              reference = reference_dataset,
                                              verbose = TRUE
)

seurat_object_list <- IntegrateData(anchorset = integration_anchors,
                                    normalization.method = "SCT",
                                    dims = seq_len(snakemake@config$n_dimred_integration),
                                    verbose = TRUE
)

DefaultAssay(integrated_data) <- "integrated"

seed <- 122334

# Perform linear dimensional reduction by PCA on HVGs
integrated_data <- RunPCA(object = integrated_data,
                     features = VariableFeatures(integrated_data),
                     seed.use = seed,
                     npcs = 50)
# Run UMAP
integrated_data <- RunUMAP(integrated_data,  dims = 1:30, seed.use = seed)

#Save seurat
#integrated_data_diet <- Seurat::DietSeurat(integrated_data, assays = c("SCT","integrated"),
#                                      dimreducs = c("pca", "umap"))

# Save the filtered, normalized and integrated data
saveRDS(integrated_data, "./seurat_integrated.rds")
print("saved seurat_integrated.rds")

sessionInfo()
date()
