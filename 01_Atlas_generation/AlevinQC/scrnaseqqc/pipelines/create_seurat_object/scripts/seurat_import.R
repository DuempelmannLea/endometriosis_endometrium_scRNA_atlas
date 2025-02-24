options(error = traceback)
# Samplewise importing the count matrices into a Seurat object 
# log
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")
# image
save.image(file = paste0(snakemake@log[[2]]) )

# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(Matrix)
  library(biomaRt)
  library(testthat)
})

##################
# load paths
##################
print(snakemake@input)
print(snakemake@output)

###########################
# load helper functions
###########################
source(paste(snakemake@scriptdir, "helperfunctions", "general.R", sep = "/"))

######################################################################
# import alevin quants and create Seurat object
# parameters min.cells and min genes are filtering thresholds
# project is defined as the read sample name given from snakemake 
######################################################################

path_data <- paste0(snakemake@input[[1]])
path_lookup <- paste0(snakemake@input[[2]])

data <- readRDS(path_data)
gene_lookup_table <- readRDS(path_lookup)

ifelse(snakemake@config$counting == TRUE, 
       count_matrix <- data$counts,
       count_matrix <- data)
#########################################################################################################
# rename the gene names to external gene names and create lookup table for gene names
#########################################################################################################
count_matrix  <- rename_features(matrix=count_matrix, 
                                 lookup=gene_lookup_table, 
                                 config = snakemake@config )

# create a seurat object from count matrix
data <- CreateSeuratObject(counts = count_matrix, 
                           min.cells = snakemake@config$min.cells, 
                           min.features = snakemake@config$min.genes, 
                           project = paste0(snakemake@wildcards$sample))

# sanity check: Object is a Seurat object,  not empty, contains features and cells in the raw data slot
stopifnot({
  expect_is(data, "Seurat")
  expect_gt(dim(data@assays$RNA)[1], 1)
  expect_gt(dim(data@assays$RNA)[2], 1)
})
# save results
saveRDS(data, file = snakemake@output[[1]])

sessionInfo()
date()

