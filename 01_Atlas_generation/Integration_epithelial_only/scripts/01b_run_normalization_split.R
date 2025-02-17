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
})


# load seurat object
seurat_obj <- readRDS(snakemake@input[[1]])

# Subset epithelial cells
Idents(seurat_obj) <- seurat_obj@meta.data$AnnotationMain
seurat_obj <- subset(seurat_obj, idents = "epithelial")

#Remove samples with less than 100 cells, otherwise the subsequent integration will fail
samples_cells <- as.data.frame(table(seurat_obj@meta.data[["sample"]])) 
print(samples_cells)
samples_integration <- samples_cells %>% filter(Freq > 99)
samples_integration <- samples_integration$Var1
print("samples used for integration")
print(samples_integration)
seurat_obj <- subset(seurat_obj, idents = samples_integration)

####################################
## 01. Normalization split per batch
####################################

regr.vars <- strsplit(snakemake@params[[1]], split='_')[[1]]

# split by batch
Idents(seurat_obj) <- seurat_obj@meta.data$sample
seurat_obj <- SCTransform(seurat_obj,
                          vars.to.regress = regr.vars,
                          conserve.memory=FALSE,
                          seed.use = 123,
                          vst.flavor = "v2",
                          return.only.var.genes=TRUE)

seurat_obj_split <- SplitObject(seurat_obj, split.by = snakemake@params[[2]])

saveRDS(seurat_obj_split, snakemake@output[[1]])
sessionInfo()
date()
