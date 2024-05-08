####################################
## 00. Load library and environment
####################################

# load packages
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(future)
  library(RColorBrewer)
})


# load seurat object
seurat_obj <- readRDS(./seurat.rds)
Idents(seurat_obj) <- seurat_obj@meta.data$sample

####################################
## 01. Normalization split per batch
####################################

regr.vars <- c("percent.mt","sample")

# split by batch
seurat_obj <- SCTransform(seurat_obj,
                          vars.to.regress = regr.vars,
                          conserve.memory=FALSE,
                          seed.use = 123,
                          vst.flavor = "v2",
                          return.only.var.genes=TRUE)

seurat_obj_split <- SplitObject(seurat_obj, split.by = "batch")

saveRDS(seurat_obj_split, "./seurat_obj_split.rds")
sessionInfo()
date()
