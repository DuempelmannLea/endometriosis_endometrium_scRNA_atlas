# Aim : Code used for the generation of figure 4 
# Author: Shaoline Sheppard

## ---------------------------------- ##
## Setup
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

# input path 
input_path <- '../Tables' 

# source common functions 
source('./common_functions.R')

## ---------------------------------- ##
##  Figure 4.b
## ---------------------------------- ##

#1. Get plot for validation 
#load performance summary table validation 
validation_metrics <- read.csv(file.path(input_path,'Performance_all_cells_validation.csv'))
validation_metrics$dimred <- factor(validation_metrics$dimred,levels = c('pca','corrgpca','hvg','diffgen'))
validation_metrics_list <- split(validation_metrics, validation_metrics$dimred)
list_plot_validation_auc <- lapply(validation_metrics_list,plot_boxplot,y_variable='auc')


#2. Get plot for evaluation
#load performance summary table evaluation  
evaluation_metrics <- read.csv(file.path(input_path,'Performance_all_cells_evaluation.csv'))
evaluation_metrics$dimred <- factor(evaluation_metrics$dimred,levels = c('pca','corrgpca','hvg','diffgen'))
evaluation_metrics_list <- split(evaluation_metrics, evaluation_metrics$dimred)
list_plot_evaluation_auc <- lapply(evaluation_metrics_list,plot_scatterplot,y_variable='auc')

#3. Combine the 2 plots
dimreds <- names(list_plot_validation_auc)
for (i in 1:length(dimreds)){
  list_plot_validation_auc[[i]]$title <- dimreds[i]
}
list_plot_validation_auc <- lapply(list_plot_validation_auc,function(x){x + ggtitle(x$title)})
plot_grid(plotlist = c(list_plot_validation_auc,
                       list_plot_evaluation_auc),nrow = 2)
validation_metrics %>%
  group_by(dimred) %>%
  summarize(median(auc))
## ---------------------------------- ##
##  Figure 4.c
## ---------------------------------- ##

validation_auc <- read.csv(file.path(input_path,'Performance_all_cells_signature_confirmation.csv')) 
validation_auc$dimred <- factor(validation_auc$dimred,levels = c("PCA","corrGPCA" ,"corrGPCA + PCA" ))
plot_boxplot_scailyte_custom(validation_auc)

## ---------------------------------- ##
##  Figure 4.d
## ---------------------------------- ##

pca_signature <- read.csv(file.path(input_path,'SupplementaryTable6_PCAsignature.csv')) %>%
  dplyr::pull(gene)
corrgpca_signature <- read.csv(file.path(input_path,'SupplementaryTable6_CorrgPCA_signature.csv')) %>%
  dplyr::pull(gene)

go_results_all_cells <- gprofiler2::gost(intersect(pca_signature,corrgpca_signature), 
                                         domain_scope = "annotated", 
                                         significant = FALSE)

go_enrichment_plot(go_results_all_cells)

## ---------------------------------- ##
##  Figure 4.e
## ---------------------------------- ##

## predictive cells enrichment 


## ---------------------------------- ##
##  Figure 4.f
## ---------------------------------- ##

validation_auc <- read.csv(file.path(input_path,'Performance_11_gene_signature.csv')) 
plot_boxplot_scailyte_custom(validation_auc,accross_dimreds = F)

