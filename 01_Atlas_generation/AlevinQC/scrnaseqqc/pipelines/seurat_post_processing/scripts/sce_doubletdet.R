#log
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")
#image
save.image(file = paste0(snakemake@log[[2]]) )

# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(scater)
  library(SingleCellExperiment)
  library(scran)
  library(Seurat)
  library(Matrix)
  library(scds)
  library(scDblFinder)
  library(BiocParallel)
  
})

##################
# load paths
##################
print(snakemake@input[[1]])
print(snakemake@output[[1]])
path_data <- paste0( snakemake@input[[1]] )
set.seed(100)

###########################################################
# Normalistaion, scaling, HVG detection, rename genes
###########################################################

# Load the dataset
raw_data <- readRDS( file = path_data )
# as sce object
sce_raw <- as.SingleCellExperiment(raw_data)
sce_raw <- scater::logNormCounts(sce_raw)
# Use the scDblFinder package 
sce_raw <- scDblFinder(sce_raw, 
                       samples="sample",
                       BPPARAM=MulticoreParam(length(unique(sce_raw$sample)))
                       )
table(sce_raw$scDblFinder.class, sce_raw$sample)

##  Using scds fromkostaklab, for more info see https://bioconductor.org/packages/release/bioc/vignettes/scds/inst/doc/scds.html 
##  Computational doublet annotation, cxds: by co-expression of gene pairs and works with absence/presence calls only, bcds: xgbosst classfier
sce_raw = tryCatch({cxds(sce_raw,retRes = FALSE)}, error=function(e) {
  message(e)
  sce_raw$cxds_score = NA
  return(sce_raw)
} )
sce_raw = tryCatch(bcds(sce_raw,retRes = FALSE,verb=TRUE), error = function(e){
  message(e)
  sce_raw$bcds_score = NA
  return(sce_raw)
})
sce_raw = tryCatch(cxds_bcds_hybrid(sce_raw), error = function(e) {
  message(e)
  sce_raw$hybrid_score = NA
  return(sce_raw)
})
# Combine Scores
df_doublet_scores <- cbind(cell_ident=colnames(sce_raw), 
                       scDblFinder.score = sce_raw$scDblFinder.score,
                       scDblFinder.class = sce_raw$scDblFinder.class,
                       scds.cxds.score = as.numeric(sce_raw$cxds_score),
                       scds.bcds.score = as.numeric(sce_raw$bcds_score), 
                       scds.hybrid.score = as.numeric(sce_raw$hybrid_score) )
rownames(df_doublet_scores) <- df_doublet_scores[,1]
# add to metadata to seurat object
raw_data = Seurat::AddMetaData(raw_data, df_doublet_scores, col.name = colnames(df_doublet_scores))
raw_data$scDblFinder.score <- as.numeric(raw_data$scDblFinder.score)
raw_data$scds.hybrid.score <- as.numeric(raw_data$scds.hybrid.score)

# Save the filtered, normalized and scaled data with annotaed doublets
saveRDS(raw_data, snakemake@output[[1]])

sessionInfo()
date()


