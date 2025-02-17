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
Idents(seurat_obj) <- seurat_obj@meta.data$sample


####################################
## 01. Normalization split per batch
####################################

regr.vars <- strsplit(snakemake@params[[1]], split='_')[[1]]

# split by batch
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
