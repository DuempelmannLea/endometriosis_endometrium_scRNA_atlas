# Creation of QC reports of the Salmon quantification step using C.Soneson AlevinQC package
#log
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")
#image
save.image(file = paste0(snakemake@log[[2]]) )

# Load libraries
suppressPackageStartupMessages({
  library(alevinQC)
})

##################
# load paths
##################
print(snakemake@input[[1]])
print(snakemake@output[[1]])
path <- paste0( snakemake@input[[1]] )
baseDir <- dirname(path)
checkAlevinInputFiles(baseDir = dirname(baseDir) )

#######################################################################
# export summary metrics from salmon Alevin Quantification Step
######################################################################

alevin <- readAlevinQC(baseDir = dirname(baseDir))
median_mapping_rate <- median(alevin$cbTable$mappingRate, na.rm = TRUE)
alevin$summaryTables$mappingRate <- median_mapping_rate
final_whitelist <- alevin$summaryTables$finalWhitelist

#######################################################################
# Compute alevinQC results, create the alevin QC , assumes fixed paths
#######################################################################

alevinQCReport(baseDir = dirname(baseDir), 
               sampleId =  snakemake@wildcards[[1]], 
               forceOverwrite = TRUE, ignorePandoc=TRUE,
               outputFile = paste0( as.character( snakemake@wildcards[[1]] ),"_alevinReport.html"), 
               outputFormat = "html_document",outputDir = paste0(dirname(snakemake@output[[1]]))
               ) 

saveRDS(alevin, snakemake@output[[2]])

## Session info

sessionInfo()
date()
