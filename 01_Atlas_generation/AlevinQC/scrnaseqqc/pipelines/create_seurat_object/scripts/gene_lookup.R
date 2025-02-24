options(error = traceback)
# Create a lookup table containg the gene names in the ensemble and HUGO format
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
# helper functions
###########################
source(paste(snakemake@scriptdir, "helperfunctions", "general.R", sep = "/"))

######################################################################
# load data
######################################################################

path <- paste0(snakemake@input[[1]])
data <- readRDS(path)

ifelse(snakemake@config$counting == TRUE, 
       count_matrix <- data$counts,
       count_matrix <-  data)
#########################################################################################################
# create gene lookuptable using the gene annotation
#########################################################################################################
gene_look_table <- entrez_lookup(matrix = count_matrix,
                                 config = snakemake@config)

# what are the ensembl_gene_id returning NAs in entrez
ensembles_genes_na <-  gene_look_table[is.na(gene_look_table$entrezgene),]
table(ensembles_genes_na$gene_biotype)

# sanity check: Gene look table is not a empty data frame, columns entrezgene_id  andexternal_gene_name are numeric and character"
stopifnot({
  expect_is(gene_look_table, "data.frame")
  expect_gt(dim(gene_look_table)[1], 1)
  expect_is(gene_look_table$entrezgene_id, "integer")
  expect_is(gene_look_table$external_gene_name, "character")
})

# save results
saveRDS(gene_look_table, file = snakemake@output[[1]])

sessionInfo()
date()

