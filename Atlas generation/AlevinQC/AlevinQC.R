## ----v1, include = FALSE------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = TRUE,
  crop = NULL
)

## ----install, eval=FALSE------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("alevinQC")

## ----v2, message=FALSE--------------------------------------------------------
library(alevinQC)

## ----v3a, eval=FALSE, results='asis'------------------------------------------
#  baseDir
#    |--alevin
#    |    |--featureDump.txt
#    |    |--filtered_cb_frequency.txt
#    |    |--MappedUmi.txt
#    |    |--quants_mat_cols.txt
#    |    |--quants_mat_rows.txt
#    |    |--quants_mat.gz
#    |    |--raw_cb_frequency.txt
#    |    |--whitelist.txt
#    |--aux_info
#    |    |--meta_info.json
#    |--cmd_info.json

## ----v3b, eval=FALSE, results='asis'------------------------------------------
#  baseDir
#    |--alevin
#    |    |--featureDump.txt
#    |    |--raw_cb_frequency.txt
#    |    |--whitelist.txt  (depending on how alevin was run)
#    |--aux_info
#    |    |--meta_info.json
#    |    |--alevin_meta_info.json
#    |--cmd_info.json

## ----v4-----------------------------------------------------------------------
baseDir <- system.file("extdata/alevin_example_v0.14", package = "alevinQC")
checkAlevinInputFiles(baseDir = baseDir)

## ----v5, message=FALSE, warning=FALSE-----------------------------------------
outputDir <- tempdir()
alevinQCReport(baseDir = baseDir, sampleId = "testSample", 
               outputFile = "alevinReport.html", 
               outputFormat = "html_document",
               outputDir = outputDir, forceOverwrite = TRUE)

## ----v7, message=FALSE--------------------------------------------------------
app <- alevinQCShiny(baseDir = baseDir, sampleId = "testSample")

## ----v8, eval=FALSE-----------------------------------------------------------
#  shiny::runApp(app)

## -----------------------------------------------------------------------------
if (interactive()) {
  out <- shiny::runApp(app)
}

## ----v9-----------------------------------------------------------------------
alevin <- readAlevinQC(baseDir = baseDir)

## ----v10----------------------------------------------------------------------
head(alevin$cbTable)

## ----v11----------------------------------------------------------------------
knitr::kable(alevin$summaryTables$fullDataset)
knitr::kable(alevin$summaryTables$initialWhitelist)
knitr::kable(alevin$summaryTables$finalWhitelist)

## ----v12----------------------------------------------------------------------
knitr::kable(alevin$versionTable)

## ----v13, warning=FALSE-------------------------------------------------------
plotAlevinKneeRaw(alevin$cbTable)
plotAlevinBarcodeCollapse(alevin$cbTable)
plotAlevinQuant(alevin$cbTable)
plotAlevinKneeNbrGenes(alevin$cbTable)

## ----v14----------------------------------------------------------------------
sessionInfo()
