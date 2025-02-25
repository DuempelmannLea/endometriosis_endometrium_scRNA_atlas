options(error = traceback)
# Importing of count matrices from the Salmon Alevin Quantification tool using the tximport package
# log
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")
# mage
save.image(file = paste0(snakemake@log[[2]]) )

# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tximport)
  library(jsonlite)
  library(testthat)
})

##################
# load paths
##################
print(snakemake@input)
print(snakemake@output)

################################
# import counts data from Alevin
################################

# path to the given output directory of Alevin
file_path <- file.path(snakemake@input[[1]])
alevin_data <- tximport(file_path, type="alevin")

# sanity cheks
stopifnot({
    expect_is(alevin_data, "lisst")
    expect_is(alevin_data[['counts']], "dgTMatrix")
    expect_gt(ncol(alevin_data[['counts']]), 1)
    })

# save results
saveRDS(alevin_data , file = snakemake@output[[1]])

sessionInfo()
date()

