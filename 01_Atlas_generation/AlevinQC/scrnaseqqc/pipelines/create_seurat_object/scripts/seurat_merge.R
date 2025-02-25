options(error = traceback)
# Merging the sample wise Seurat objects into one Seurat object and adding the meta data
# log
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")
# image
save.image(file = paste0( snakemake@log[[2]] ) )


# Load libraries
suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(Seurat)
  library(biomaRt)
  library(testthat)
})

print(snakemake@input)
print(snakemake@output)

###########################
# load helper functions
###########################
source(paste(snakemake@scriptdir, "helperfunctions", "general.R", sep = "/"))

###########################
# load metadata
###########################
group_info <- read.csv( paste0(snakemake@config$group_info), comment.char="#")
row.names(group_info) = group_info$names
###################################################################
# Merging, metadata annotation and remaval of zero expressed genes
###################################################################

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells. Keep all cells with at
# least 200 detected genes
# add metadata for every read sample and merge the seurat objects


# read first seurat object
data_combined <- readRDS(file = snakemake@input[[1]])

# add metadata 
##  create metadata table 
df_meta_data <- data.frame(
  cell.names = as.character( colnames(data_combined) ),
  celline = rep(group_info[snakemake@params[[1]][1],"celline"]),
  lane = rep(group_info[snakemake@params[[1]][1],"lane"]),
  batch = rep(group_info[snakemake@params[[1]][1],"batch"]),
  sample = rep(group_info[snakemake@params[[1]][1],"names"]),
  treatment = rep(group_info[snakemake@params[[1]][1],"treatment"])
)

##  set cell names
stopifnot( length(unique(rownames(df_meta_data))) == length(unique(colnames(data_combined) )))
rownames(df_meta_data) <-  colnames(data_combined)

# add metadata
data_combined  <- AddMetaData(object = data_combined, metadata = df_meta_data  ) # use new function object$name <- vector
# merge by merge function, use for loop for not loading all data

for (i in 2:length(snakemake@input)) {
  print(paste0("merging seurat object ",snakemake@input[[i]], " with metadata from ",snakemake@params[[1]][i]   ))
  stopifnot( grepl(snakemake@params[[1]][i] , snakemake@input[[i]] ) )
  obj2 <- readRDS(file = snakemake@input[[i]])
  # create metadata
  df_meta_data <- data.frame(
    cell.names = as.character(colnames(obj2) ),
    celline = rep(group_info[snakemake@params[[1]][i],"celline"]),
    lane = rep(group_info[snakemake@params[[1]][i],"lane"]),
    batch = rep(group_info[snakemake@params[[1]][i],"batch"]),
    sample = rep(group_info[snakemake@params[[1]][i],"names"]),
    treatment = rep(group_info[snakemake@params[[1]][i],"treatment"])
  )
  rownames(df_meta_data) <-  colnames(obj2)
  # add metadata
  stopifnot(anyNA(df_meta_data) == FALSE )
  obj2   <- AddMetaData(object = obj2, metadata = df_meta_data  )
  # merge, append the read sample name to the cellnames (starting from i+1) # add charcter vector of cells
  cell.ids <- as.character(c(1, paste0(i)))
  data_combined  <- merge(x = data_combined, 
                          y = obj2, 
                          merge.data = TRUE,
                          add.cell.ids = cell.ids  ,
                          project = snakemake@config$project)

}

# add cell-wise metadata on perent mt expression, number of UMIs and genes per cells
data_combined <- PercentageFeatureSet(object = data_combined, 
                             pattern = paste0(snakemake@config$regex.mito), 
                             col.name = "percent.mt")

# sanity checks: Seurat object data_combined is not empty and contains more than one sample wise grouping factor
stopifnot({
  expect_is(data_combined, "Seurat")
  expect_gt(dim(data_combined)[1], 1)
  expect_gt(dim(data_combined)[2], 1)
  expect_gt(length(unique(data_combined$treatment)), 1)
})
################################
# save the combined dataset
################################
saveRDS(data_combined, file = snakemake@output$datacombined)

# Save the data summary for testing
## export csv with sum of counts, data, scale and variable features 

n_cells_total = length(data_combined@meta.data$cell.names)
n_cells_sample = table(data_combined@meta.data$sample)
n_cells_treatment = table(data_combined@meta.data$treatment)
data.merge.summary <- cbind( n_cells_sample, n_cells_treatment, n_cells_total )
write.csv(data.merge.summary, snakemake@output$testdata)

sessionInfo()
date()
