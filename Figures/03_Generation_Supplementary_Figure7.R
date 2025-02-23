## ---------------------------------- ##
##  Setup
## ---------------------------------- ##

# load librairies
library(ggplot2)
library(rlang)
library(stringr)
library(dplyr)
library(tidyr)
library(cowplot)
library(gridExtra)
library(patchwork)
library(tibble)
library(tinter)
library(Seurat)

# input path 
input_path <- '../Tables'

# source common functions 
source('./common_functions.R')

## ---------------------------------- ##
##  Figure 6.a
## ---------------------------------- ##

## TO do lea 

## ---------------------------------- ##
##  Figure 6.b
## ---------------------------------- ##

# this figure represent a scheme and was done manually 

## ---------------------------------- ##
##  Figure 6.c
## ---------------------------------- ##

#1. Get plot for evaluation
#load performance summary table evaluation  
evaluation_metrics <- read.csv(file.path(input_path,'Performance_cycling_fibroblasts_evaluation.csv'))
evaluation_metrics$dimred <- factor(evaluation_metrics$dimred,levels = c('pca','corrgpca','hvg','diffgen'))
evaluation_metrics$fold <- gsub(".*(cvsplit\\d+).*", "\\1", evaluation_metrics$learner_uid)
evaluation_metrics_list <- split(evaluation_metrics, evaluation_metrics$dimred)
list_plot_evaluation_auc <- lapply(evaluation_metrics_list,plot_scatterplot,y_variable='auc')
dimreds <- names(list_plot_evaluation_auc)
for (i in 1:length(dimreds)){
  list_plot_evaluation_auc[[i]]$title <- dimreds[i]
}
list_plot_evaluation_auc <- lapply(list_plot_evaluation_auc,function(x){x + ggtitle(x$title)})
plot_grid(plotlist = c(list_plot_evaluation_auc),nrow = 1)

## ---------------------------------- ##
##  Figure 6.d
## ---------------------------------- ##

validation_auc_cf <- read.csv(file.path(input_path,'Performance_cycling_fibroblasts_signature_confirmation.csv')) 
validation_auc_cf$dimred <- factor(validation_auc_cf$dimred,levels = c("PCA","corrGPCA" ,"corrGPCA+PCA" ))
plot_boxplot_scailyte_custom(validation_auc_cf)

## ---------------------------------- ##
##  Figure 6.e
## ---------------------------------- ##

pca_signature <- read.csv(file.path(input_path,'PCA_cycling_fibroblast_signature.csv')) %>%
  dplyr::pull(gene)
corrgpca_signature <- read.csv(file.path(input_path,'Corrgpca_cycling_fibroblast_signature.csv')) %>%
  dplyr::pull(gene)

go_results_cf <- gprofiler2::gost(intersect(pca_signature,corrgpca_signature), 
                                         domain_scope = "annotated", 
                                         significant = FALSE)
go_enrichment_plot(go_results_cf,top_gene_threshold=20)
