options(error = traceback)
# Creation of sparse Matrix from counts from Alevin or Cellrangers mtx format
# log
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")
# mage
save.image(file = paste0(snakemake@log[[2]]) )

# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(tximeta)
  library(fishpond)
  library(testthat)
  library(SummarizedExperiment)
})

##################
# load paths
##################
print(snakemake@input)

###########################
# load helper functions
###########################
source(paste(snakemake@scriptdir, "helperfunctions", "general.R", sep = "/"))
############################################
# import counts from cellranger, mmx format
############################################
# import counts

import_counts <- function(snakemake){
  if (snakemake@params$count_format_convention == "alevin"){
    stopifnot(file.exists(snakemake@input[[1]]))
    summarised_experiment <- tximeta(snakemake@input[[1]], type="alevin",alevinArgs=list(filterBarcodes=TRUE))
    count_matrix = assay(summarised_experiment)
  } else if (snakemake@params$count_format_convention == "mtx" | snakemake@params$count_format_convention == "10X_new" | snakemake@params$count_format_convention == "10X_old") {
    count_matrix<- as.matrix(Seurat::Read10X(dirname(snakemake@input[[1]])))
  } else {
    print("Wrong file format")
    stop()
  }
}

count_matrix <- import_counts(snakemake)

# save results
saveRDS( count_matrix, file = snakemake@output[[1]] )

sessionInfo()
date()

