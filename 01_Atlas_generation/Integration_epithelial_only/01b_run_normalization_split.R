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
seurat_obj <- readRDS("../_Data/EndoAtlas.rds")

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

# SCTransform and regress out percent.mt and sample
regr.vars <- c("percent.mt","sample")
seurat_obj <- SCTransform(seurat_obj,
                          vars.to.regress = regr.vars,
                          conserve.memory=FALSE,
                          seed.use = 123,
                          vst.flavor = "v2",
                          return.only.var.genes=TRUE)

# split by 10x capture batch
seurat_obj_split <- SplitObject(seurat_obj, split.by = "batch")

# save split seurat object
saveRDS(seurat_obj_split, "../_Data/01_Atlas_generation/Integration_epithelial_only/seurat_obj_split.rds")
sessionInfo()
date()
